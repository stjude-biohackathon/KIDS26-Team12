#!/usr/bin/env Rscript
# =============================================================================
# scripts/loco_one_fold.R - run ONE outer leave-one-cancer-out fold.
# =============================================================================
#
# WHY THIS EXISTS
#   scripts/train_baseline.R runs all 30 outer folds in a serial loop. Measured
#   on the real matrix that is ~45 minutes per fold, or roughly 22 hours end to
#   end (see results/benchmark/preprocess_benchmark.tsv).
#
#   The folds are completely independent of one another: fold "BRCA" never reads
#   anything produced by fold "LUAD". That makes this an embarrassingly parallel
#   problem, so each fold can run as its own LSF array task and the wall time
#   collapses to roughly the cost of a single fold.
#
#   This script is deliberately a THIN WRAPPER around the same fit_en() /
#   predict_en() calls that train_baseline.R makes. It does not reimplement the
#   modelling. Keeping one implementation means the parallel path cannot quietly
#   drift away from the serial one.
#
# USAGE
#   Rscript scripts/loco_one_fold.R <beta.tsv> <master.tsv> <out_dir> <fold_index> [transform] [feature_rank]
#
#   fold_index is 1-based and indexes into the sorted vector of development
#   cancer types, which is exactly the order train_baseline.R iterates in. LSF
#   array indices are also 1-based, so LSB_JOBINDEX maps across directly.
#
#   transform is optional and defaults to "identity" (the run 01 behaviour).
#   Set "log1p" to fit on log(1+HRDsum) and back-transform - see C2 in docs/21
#   and R/calibration.R for why this is a real trade-off and not a free win.
#   "clip" is NOT offered here because clipping is a post-hoc transform of the
#   predictions and needs no refit; the merge step applies it.
#
#   feature_rank is optional and defaults to "pooled" (the run 01 behaviour:
#   the 5,000-probe unsupervised filter ranks probes by TOTAL variance across
#   the training fold). Set "within_tissue" for V3-abs, which ranks the same
#   5,000 probes by POOLED WITHIN-TISSUE variance instead, so probes that only
#   separate tissues from one another cannot win a slot. Set "lineage_penalized"
#   for V4 (docs/27 section 6): a SUPERVISED within-tissue HRD meta-association
#   score minus a lineage-discriminability penalty whose weight lambda_pen is
#   selected inside the inner folds by the same macro-averaged MAE that picks
#   alpha and lambda. Everything else about the fold - missingness filter,
#   medians, centre/scale, tuning grid, target - is unchanged. Whichever is used
#   is recorded in metrics_<TYPE>.tsv, and on the V4 path the selected
#   lambda_pen is recorded there too.
#
# OUTPUT
#   One .rds bundle and two .tsv files per fold, written into <out_dir>/folds/.
#   scripts/loco_merge.R then combines them into the same files train_baseline.R
#   would have produced.
#
# WHAT THIS DOES NOT DO
#   No Phase 2 freezing, no calibration reservation, no conformal interval.
#   Those need all folds complete and belong in the merge step.
# =============================================================================

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4 || length(args) > 7) {
  stop("Usage: Rscript scripts/loco_one_fold.R <beta.tsv> <master.tsv> <out_dir> <fold_index> [transform] [feature_rank] [lambda_pen_grid]")
}
beta_path <- args[1]; meta_path <- args[2]; out_dir <- args[3]
fold_index <- as.integer(args[4])
transform_name <- if (length(args) >= 5) args[5] else "identity"
feature_rank <- if (length(args) >= 6) args[6] else "pooled"
# lambda_pen_grid is a comma-separated list of candidate lineage penalties and
# is meaningful ONLY on the lineage_penalized path. It defaults to the grid
# named in docs/27 section 6. Whatever is supplied is SEARCHED IN THE INNER
# FOLDS - passing a grid does not pick a value, it only says which values the
# inner-fold macro MAE is allowed to choose between.
lambda_pen_grid <- if (length(args) >= 7) {
  g <- as.numeric(strsplit(args[7], ",", fixed=TRUE)[[1]])
  if (!length(g) || any(!is.finite(g))) stop("lambda_pen_grid must be finite comma-separated numbers")
  g
} else c(0, 0.25, 0.5, 1, 2)
if (is.na(fold_index) || fold_index < 1L) stop("fold_index must be a positive integer")
if (!feature_rank %in% c("pooled", "within_tissue", "lineage_penalized")) {
  stop("feature_rank must be 'pooled', 'within_tissue' or 'lineage_penalized'")
}
# A grid on an unsupervised rank would be silently ignored, which is exactly the
# kind of quiet no-op that makes a run unauditable. Refuse it.
if (length(args) >= 7 && feature_rank != "lineage_penalized") {
  stop("lambda_pen_grid is only meaningful with feature_rank='lineage_penalized'")
}

fold_dir <- file.path(out_dir, "folds")
dir.create(fold_dir, recursive=TRUE, showWarnings=FALSE)

cat("host       :", Sys.info()[["nodename"]], "\n")
cat("fold index :", fold_index, "\n")
cat("transform  :", transform_name, "\n")
cat("feature_rank:", feature_rank, "\n")
if (feature_rank == "lineage_penalized")
  cat("lambda_pen grid (searched in the INNER folds):",
      paste(lambda_pen_grid, collapse=", "), "\n")
cat("started    :", format(Sys.time()), "\n\n")

source("R/model.R")
source("R/calibration.R")
tf <- get_transform(transform_name)
if (!requireNamespace("data.table", quietly=TRUE)) stop("Install dependencies with scripts/setup.R")

# --- Load matrix and metadata ----------------------------------------------
# This repeats the load in every array task. That is a deliberate trade: ~100 s
# and ~75 GB per task, against the complexity and fragility of a shared-memory
# scheme. With 30 tasks the total CPU cost of reloading is minutes, while the
# wall-clock saving from parallelism is many hours.
cat("reading beta matrix ...\n")
beta <- data.table::fread(beta_path, data.table=FALSE, check.names=FALSE)
if (anyDuplicated(beta[[1]])) stop("Duplicate probe IDs")
x <- t(as.matrix(beta[,-1,drop=FALSE])); storage.mode(x) <- "double"; colnames(x) <- beta[[1]]
rm(beta); invisible(gc(verbose=FALSE))
if (any(is.finite(x) & (x < 0 | x > 1))) stop("Expected beta values in [0,1]")
cat(sprintf("  %d samples x %d probes\n", nrow(x), ncol(x)))

meta <- read.delim(meta_path, check.names=FALSE, stringsAsFactors=FALSE)
if (anyDuplicated(meta$sample_id) || anyDuplicated(meta$patient_id)) stop("Duplicate patient/specimen")
if (any(!rownames(x) %in% meta$sample_id)) stop("Missing sample metadata")
meta <- meta[match(rownames(x), meta$sample_id),]
if (any(!is.finite(meta$HRDsum)) || any(meta$HRDsum < 0)) stop("Invalid target")
if (!all(abs(meta$HRDsum - meta$HRD_LOH - meta$LST - meta$TAI) < 1e-6)) stop("Label sum mismatch")

# --- Resolve which cancer type this task owns ------------------------------
# sort() must match train_baseline.R exactly, or task N here would not be the
# same fold as iteration N there.
cns <- meta$cancer_type %in% c("GBM","LGG")
# DEFECT 3: the hardcoded list above and the `partition` column that
# build_master.py writes are two sources of truth for the SAME lock. Compare
# them before either is used. This changes nothing about which samples are
# locked; it refuses to run when the two definitions disagree, which is the only
# situation in which the hardcoded list could silently be wrong.
assert_partition_matches_cns(meta$cancer_type, meta$partition,
                             context = basename(meta_path))
types <- sort(unique(meta$cancer_type[!cns]))
if (fold_index > length(types)) {
  stop(sprintf("fold_index %d exceeds the %d development cancer types", fold_index, length(types)))
}
type <- types[fold_index]
cat(sprintf("\nthis task holds out: %s (%d of %d)\n", type, fold_index, length(types)))

# --- Fit the fold -----------------------------------------------------------
# Both masks carry !cns, so locked CNS samples appear in neither train nor test.
tr <- !cns & meta$cancer_type != type
te <- !cns & meta$cancer_type == type
cat(sprintf("train n=%d  test n=%d\n\n", sum(tr), sum(te)))

t0 <- Sys.time()
# C2: fit on the TRANSFORMED target. tf$forward is identity by default, so this
# reproduces run 01 exactly unless a transform was requested.
# feature_rank is threaded straight through: fit_en() re-ranks inside EVERY
# inner fold and again for the outer-training refit, always on that partition's
# own rows. meta$cancer_type[tr] carries the same `tr` mask as the matrix, so the
# held-out type is absent from the ranking statistic by construction.
b <- fit_en(x[tr,,drop=FALSE], tf$forward(meta$HRDsum[tr]), meta$patient_id[tr], meta$cancer_type[tr],
            feature_rank = feature_rank, lambda_pen_grid = lambda_pen_grid)
b$target_transform <- transform_name
cat(sprintf("fit_en completed in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units="mins"))))
cat(sprintf("  selected alpha=%.2g lambda=%.5g  (boundary: %s)\n",
            b$alpha, b$lambda, b$lambda_boundary_side))
if (!is.null(b$lambda_pen)) {
  cat(sprintf("  selected lambda_pen=%.4g from grid {%s} (inner-fold macro MAE)\n",
              b$lambda_pen, paste(b$lambda_pen_grid, collapse=", ")))
}

p <- predict_en(b, x[te,,drop=FALSE])
# Back-transform IMMEDIATELY so every downstream metric, null and artefact is on
# the original HRDsum scale. Mixing scales across folds would be silently wrong.
p$predicted_reference_HRDsum <- tf$inverse(p$predicted_reference_HRDsum)
p$actual <- meta$HRDsum[te]
p$cancer_type <- type
p$patient_id <- meta$patient_id[te]
# training_mean is stored on the transformed scale inside the bundle, so it must
# be back-transformed too before it can serve as a null on the original scale.
p$null_prediction <- tf$inverse(b$training_mean)
p$split <- "development_LOCO"

mm <- metrics(p$actual, p$predicted_reference_HRDsum)
mm$cancer_type <- type
mm$train_n <- sum(tr)
mm$target_transform <- transform_name
# Provenance: every metrics row states which probe-ranking strategy produced it,
# so a pooled run and a V3-abs run are never confusable after the fact.
mm$feature_rank <- feature_rank
# V4 only: the inner-fold-selected lineage penalty, and whether it landed on an
# endpoint of its grid (which would say the grid was too narrow).
mm$lambda_pen <- if (is.null(b$lambda_pen)) NA_real_ else b$lambda_pen
mm$lambda_pen_at_boundary <- if (is.null(b$lambda_pen_at_boundary)) NA else b$lambda_pen_at_boundary
mm$selected_features <- length(b$preprocess$features)
mm$alpha <- b$alpha
mm$lambda <- b$lambda
mm$lambda_at_boundary <- b$lambda_at_boundary
mm$lambda_boundary_side <- b$lambda_boundary_side
mm$null_MAE <- mean(abs(p$actual - b$training_mean))
mm$abstention_rate <- mean(!p$reportable)

# Per-fold null panel. Within one fold every sample is the same cancer type, so
# the tissue-mean null collapses to this fold's own mean. The between-tissue
# picture only appears once folds are pooled in the merge step.
np <- null_panel(p$actual, p$predicted_reference_HRDsum, p$cancer_type, b$training_mean)
np$cancer_type <- type

# --- Write per-fold artefacts ----------------------------------------------
# One file set per fold, named by type. The merge step globs these, so a failed
# task is detectable simply by its outputs being absent.
saveRDS(b, file.path(fold_dir, paste0("loco_", type, ".rds")))
write.table(b$inner_folds, file.path(fold_dir, paste0("inner_folds_", type, ".tsv")),
            sep="\t", row.names=FALSE, quote=FALSE)
write.table(p, file.path(fold_dir, paste0("predictions_", type, ".tsv")),
            sep="\t", row.names=FALSE, quote=FALSE)
write.table(mm, file.path(fold_dir, paste0("metrics_", type, ".tsv")),
            sep="\t", row.names=FALSE, quote=FALSE)
write.table(np, file.path(fold_dir, paste0("nullpanel_", type, ".tsv")),
            sep="\t", row.names=FALSE, quote=FALSE)

cat("\nwrote artefacts for fold:", type, "\n")
cat("finished   :", format(Sys.time()), "\n")

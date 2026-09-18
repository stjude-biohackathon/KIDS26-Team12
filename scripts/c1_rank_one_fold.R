#!/usr/bin/env Rscript
# =============================================================================
# scripts/c1_rank_one_fold.R - run ONE leave-one-cancer-out fold of the
#                              TISSUE-RELATIVE (rank) model.
# =============================================================================
#
# WHY THIS EXISTS
#   scripts/loco_one_fold.R runs one outer fold of the ABSOLUTE HRDsum model.
#   This is its counterpart for the relative-target model in R/rank_model.R: same
#   fold indexing, same positional CLI, same per-fold artefact layout, so the two
#   probes can be compared fold for fold. It is a THIN WRAPPER - all modelling
#   lives in R/rank_model.R.
#
# USAGE
#   Rscript scripts/c1_rank_one_fold.R <beta.tsv> <master.tsv> <out_dir> <fold_index> \
#                                      [target_type] [weight_mode] [feature_rank]
#
#   fold_index  1-based, indexes sort(unique(cancer_type[!cns])) - the same
#               ordering loco_one_fold.R uses, so LSB_JOBINDEX maps across.
#   target_type "normal_score" (default) or "centered"
#   weight_mode "tissue" (default) or "sample"
#   feature_rank "pooled" (default) or "within_tissue"
#
#   Variants:  V1 = normal_score sample  pooled
#              V2 = normal_score tissue  pooled
#              V3 = normal_score tissue  within_tissue
#
#   Example (ONE sentinel fold, relative-model probe directory):
#     Rscript scripts/c1_rank_one_fold.R data/processed/beta.tsv \
#             data/processed/master.tsv results/c1_rank_probe 1 \
#             normal_score tissue pooled
#
# OUTPUT (all under <out_dir>/folds/, default out_dir results/c1_rank_probe)
#   rank_<type>.rds           frozen bundle (relative score only - see below)
#   predictions_<type>.tsv    sample_id, cancer_type, relative_score, percentile,
#                             true relative target (EVAL-ONLY)
#   metrics_<type>.tsv        within-tissue Spearman/Pearson, macro relative MAE
#   nullpanel_<type>.tsv      within-tissue permutation null for the ordering
#   inner_folds_<type>.tsv    inner fold assignment, for audit
#
# WHAT THE OUTPUT IS NOT
#   relative_score is NOT an absolute HRDsum and cannot be converted into one.
#   The bundle stores no tissue mean/SD/quantile of y by construction.
#
# LOCKED PARTITION
#   GBM and LGG are locked. The mask is built and asserted BEFORE any target
#   transform, fit, calibration, feature selection or summary happens.
# =============================================================================

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4 || length(args) > 7) {
  stop("Usage: Rscript scripts/c1_rank_one_fold.R <beta.tsv> <master.tsv> <out_dir> <fold_index> [target_type] [weight_mode] [feature_rank]")
}
beta_path <- args[1]; meta_path <- args[2]; out_dir <- args[3]
fold_index <- as.integer(args[4])
target_type  <- if (length(args) >= 5) args[5] else "normal_score"
weight_mode  <- if (length(args) >= 6) args[6] else "tissue"
feature_rank <- if (length(args) >= 7) args[7] else "pooled"
if (is.na(fold_index) || fold_index < 1L) stop("fold_index must be a positive integer")
stopifnot(target_type %in% c("normal_score","centered"),
          weight_mode %in% c("sample","tissue"),
          feature_rank %in% c("pooled","within_tissue"))

fold_dir <- file.path(out_dir, "folds")
dir.create(fold_dir, recursive=TRUE, showWarnings=FALSE)

cat("host        :", Sys.info()[["nodename"]], "\n")
cat("fold index  :", fold_index, "\n")
cat("target_type :", target_type, "\n")
cat("weight_mode :", weight_mode, "\n")
cat("feature_rank:", feature_rank, "\n")
cat("started     :", format(Sys.time()), "\n\n")

source("R/rank_model.R")   # which itself sources R/model.R
if (!requireNamespace("data.table", quietly=TRUE)) stop("Install dependencies with scripts/setup.R")

# --- Load matrix and metadata ----------------------------------------------
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

# =============================================================================
# LOCKED PARTITION - ASSERTED BEFORE ANY TRANSFORM, FIT, OR SUMMARY
# =============================================================================
cns <- meta$cancer_type %in% c("GBM","LGG")
types <- sort(unique(meta$cancer_type[!cns]))
if (fold_index > length(types)) {
  stop(sprintf("fold_index %d exceeds the %d development cancer types", fold_index, length(types)))
}
type <- types[fold_index]

# Both masks carry !cns, so locked CNS samples appear in neither train nor test.
tr <- !cns & meta$cancer_type != type
te <- !cns & meta$cancer_type == type

# Executable enforcement. Nothing below this line may run if it fails.
stopifnot(!any(meta$cancer_type[tr] %in% c("GBM","LGG")),
          !any(meta$cancer_type[te] %in% c("GBM","LGG")))
stopifnot(!type %in% c("GBM","LGG"))
assert_cns_locked(meta$cancer_type[tr], meta$cancer_type[te], type)
stopifnot(sum(tr) > 0, sum(te) > 0)

cat(sprintf("\nthis task holds out: %s (%d of %d)\n", type, fold_index, length(types)))
cat(sprintf("train n=%d  test n=%d  (CNS locked out: %d)\n\n", sum(tr), sum(te), sum(cns)))

# --- Fit the fold -----------------------------------------------------------
t0 <- Sys.time()
b <- fit_rank_en(x[tr,,drop=FALSE], meta$HRDsum[tr], meta$patient_id[tr],
                 meta$cancer_type[tr], target_type=target_type,
                 weight_mode=weight_mode, feature_rank=feature_rank)
cat(sprintf("fit_rank_en completed in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units="mins"))))
cat(sprintf("  selected alpha=%.2g lambda=%.5g  (boundary: %s)  objective=%s\n",
            b$alpha, b$lambda, b$lambda_boundary_side, b$tuning_objective))

p <- predict_rank(b, x[te,,drop=FALSE])
p$cancer_type <- type
p$patient_id <- meta$patient_id[te]
# EVAL-ONLY: the held-out tissue's own relative target, computed after the fact
# from labels that are unavailable at deployment. It is a grading key and never
# touched the fit - the bundle `b` was frozen before this line ran.
p$true_relative_target <- as.numeric(relative_target(meta$HRDsum[te],
                                                     meta$cancer_type[te], target_type))
p$split <- "development_LOCO_relative"

mm <- rank_metrics(p$true_relative_target, p$relative_score, p$cancer_type)
mm$cancer_type <- type
mm$train_n <- sum(tr)
mm$target_type <- target_type
mm$weight_mode <- weight_mode
mm$feature_rank <- feature_rank
mm$selected_features <- length(b$preprocess$features)
mm$alpha <- b$alpha
mm$lambda <- b$lambda
mm$lambda_at_boundary <- b$lambda_at_boundary
mm$tuning_objective <- b$tuning_objective
mm$abstention_rate <- mean(!p$reportable)

# --- Null comparison --------------------------------------------------------
# Within-tissue permutation of the true relative target against fixed scores.
np <- rank_null_panel(p$true_relative_target, p$relative_score, p$cancer_type)
np$cancer_type <- type
# A shuffled-score null for reference: the same scores in random order.
set.seed(b$seed)
shuffled <- p$relative_score[sample.int(nrow(p))]
np$shuffled_score_spearman <- if (nrow(p) > 2 && sd(p$true_relative_target) > 0 && sd(shuffled) > 0)
  suppressWarnings(cor(p$true_relative_target, shuffled, method="spearman")) else NA_real_

# --- Write per-fold artefacts ----------------------------------------------
saveRDS(b, file.path(fold_dir, paste0("rank_", type, ".rds")))
write.table(b$inner_folds, file.path(fold_dir, paste0("inner_folds_", type, ".tsv")),
            sep="\t", row.names=FALSE, quote=FALSE)
write.table(p, file.path(fold_dir, paste0("predictions_", type, ".tsv")),
            sep="\t", row.names=FALSE, quote=FALSE)
write.table(mm, file.path(fold_dir, paste0("metrics_", type, ".tsv")),
            sep="\t", row.names=FALSE, quote=FALSE)
write.table(np, file.path(fold_dir, paste0("nullpanel_", type, ".tsv")),
            sep="\t", row.names=FALSE, quote=FALSE)

cat("\nwrote relative-model artefacts for fold:", type, "\n")
cat("finished    :", format(Sys.time()), "\n")

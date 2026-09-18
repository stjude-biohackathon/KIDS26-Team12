#!/usr/bin/env Rscript
# =============================================================================
# scripts/freeze_final_model.R - fit and FREEZE the one shipping model.
# =============================================================================
#
# WHAT THIS IS
#   The Phase 2 half of scripts/train_baseline.R, extracted and run on its own.
#   It fits ONE elastic net on ALL 7,065 development samples (30 non-CNS cancer
#   types), derives a split-conformal interval on a reserved calibration subset,
#   and writes the inference bundle that scripts/predict_frozen.R consumes to
#   score the locked GBM/LGG cohort exactly once.
#
# WHY IT IS NOT train_baseline.R
#   train_baseline.R runs 30 serial LOCO folds (~22 h) BEFORE reaching Phase 2.
#   We already have those folds: results/loco_run01 (array 323078995), merged by
#   scripts/loco_merge.R. Re-running them to reach the freeze would cost most of
#   a day and change nothing. This script therefore runs Phase 2 ONLY. It is a
#   thin wrapper around the same fit_en()/predict_en()/conformal_q() calls, so
#   the frozen artefact cannot drift from the architecture run 01 evaluated.
#
#   Consequence, stated plainly: this script produces NO development metrics.
#   The performance numbers that describe this architecture are run 01's, and
#   run 01's folds are different fits (they each drop one tissue; this one drops
#   the calibration reservation). macro_metrics.txt does not describe this .rds.
#
# ARCHITECTURE - frozen, matches run 01 exactly
#   feature_rank    = "pooled"     (V3-abs / within_tissue is NOT being frozen)
#   target at fit   = identity     (no log1p; REJECTED on the full 30-fold run)
#   target at score = clip         (C2 ADOPTED: pmax(pred, 0))
#   inner folds     = one per TRAINING cancer type (cross-tissue tuning rule)
#   standardize     = FALSE        (data is pre-standardised by apply_preprocess)
#   seed            = 260910       (the project model seed)
#
#   The "identity at fit / clip at score" pair is expressed as the single
#   bundle field `target_transform = "clip"`. See the block at THE C2 FIELD
#   below for why that exact string, and not "identity", is required.
#
# THE LOCK
#   GBM + LGG (642 samples) are excluded from EVERY step here: preprocessing,
#   feature selection, hyperparameter selection, the fit, and the conformal
#   calibration set. Their HRDsum is never read. Three independent guards:
#     1. assert_partition_matches_cns() - the hardcoded GBM/LGG list must agree
#        row-for-row with build_master.py's `partition` column, both directions.
#     2. An explicit post-split assertion that no locked row appears in `tr` or
#        `cal`, and that tr + cal + locked partitions the cohort exactly.
#     3. final_partitions.tsv records the role of every one of the 7,707 rows,
#        so the claim is auditable from disk afterwards rather than from this
#        comment.
#
# USAGE
#   Rscript scripts/freeze_final_model.R <beta.tsv> <master.tsv> <out_dir> [--allow-fixture]
#
#   Run from the repository root: source("R/model.R") is a relative path.
#   --allow-fixture is for the synthetic dry run ONLY. It downgrades the
#   provenance gate to a recorded admission and stamps the manifest
#   IS_FIXTURE=TRUE so a fixture bundle can never be mistaken for the artefact.
#
# OUTPUT (all into <out_dir>)
#   frozen_nonCNS.rds        the shipping bundle
#   frozen_coefficients.tsv  non-zero coefficients, sorted by |coefficient|
#   final_partitions.tsv     role of every sample: locked_CNS/calibration/training
#   freeze_manifest.tsv      provenance: commit, seed, hyperparameters, sha256
#   sessionInfo.txt          the R environment that produced it
# =============================================================================

t_start <- Sys.time()

args <- commandArgs(trailingOnly = TRUE)
allow_fixture <- "--allow-fixture" %in% args
args <- args[args != "--allow-fixture"]
if (length(args) != 3) {
  stop("Usage: Rscript scripts/freeze_final_model.R <beta.tsv> <master.tsv> <out_dir> [--allow-fixture]")
}
beta_path <- args[1]; meta_path <- args[2]; out <- args[3]

cat("=============================================================\n")
cat("KIDS26 Team 12 - FREEZE THE SHIPPING MODEL\n")
cat("=============================================================\n")
cat("host       :", Sys.info()[["nodename"]], "\n")
cat("started    :", format(t_start), "\n")
cat("beta       :", beta_path, "\n")
cat("master     :", meta_path, "\n")
cat("out_dir    :", out, "\n")
cat("fixture ok :", allow_fixture, "\n\n")

source("R/model.R")
source("R/provenance.R")
# C2 clipping lives in R/calibration.R as a matched (forward, inverse) pair and
# is reached through get_transform(), the SAME object loco_one_fold.R and
# predict_frozen.R use. Re-spelling pmax(p, 0) here would recreate the
# two-sources-of-truth problem this repo already has with the CNS lock.
source("R/calibration.R")
if (!requireNamespace("data.table", quietly = TRUE)) stop("Install dependencies with scripts/setup.R")
if (!requireNamespace("glmnet", quietly = TRUE)) stop("Install glmnet via scripts/setup.R")

SEED         <- 260910L
FEATURE_RANK <- "pooled"
# THE C2 FIELD ---------------------------------------------------------------
# scripts/predict_frozen.R line 299 reads:
#     bundle_transform<-if(is.null(b$target_transform))"clip" else b$target_transform
# and immediately does tf <- get_transform(bundle_transform);
#     p$predicted_reference_HRDsum <- tf$inverse(p$predicted_reference_HRDsum_raw)
#
# R/calibration.R defines the pair as
#     clip = list(forward = function(y) y, inverse = function(p) pmax(p, 0))
#
# So "clip" IS "identity at fit time, non-negativity clipping at prediction
# time" - exactly the adopted C2 rule - expressed in one field. Writing
# "identity" here would still be clipped by predict_frozen.R's belt-and-braces
# pmax() on the following line, but the recorded clipping_rule string would then
# read inverse_identity, i.e. the artefact would misdescribe its own contract.
# The bundle must say what it means.
TARGET_TRANSFORM <- "clip"
CAL_FRACTION     <- 0.2   # matches train_baseline.R's Phase 2 exactly
COVERAGE         <- 0.95  # conformal_q() default, stated explicitly for the manifest

tf <- get_transform(TARGET_TRANSFORM)
stopifnot(identical(tf$forward(c(-1, 0, 3)), c(-1, 0, 3)))   # fit-time is identity
stopifnot(identical(tf$inverse(c(-1, 0, 3)), c(0, 0, 3)))    # score-time clips

dir.create(out, recursive = TRUE, showWarnings = FALSE)

# --- Provenance gate --------------------------------------------------------
# Same gate the fold scripts run, via the shared R/provenance.R implementation.
# An unlabelled or engineering-only matrix is a hard stop: freezing a fixture
# would produce a .rds whose schema is indistinguishable from the real artefact,
# and this one is meant to be loaded months from now by someone who was not here.
matrix_pr <- read_matrix_provenance(beta_path)
gate <- gate_matrix_provenance(matrix_pr, allow_fixture = allow_fixture)
cat("provenance :", gate$class, "-", gate$reason, "\n")
if (isTRUE(gate$fatal)) {
  cat("*** THIS IS A FIXTURE RUN. The bundle it writes is NOT the shipping artefact. ***\n")
}

# --- Load the beta matrix (memory-safe pattern from loco_one_fold.R) --------
# Probes in ROWS, samples in COLUMNS, first column probe_id. check.names=FALSE
# preserves TCGA barcodes verbatim (R would mangle the hyphens into dots).
#
# NEVER add select= here. `select=` filters COLUMNS by name; probes are ROWS, so
# every probe ID is skipped WITH A WARNING and the transpose below silently
# yields 0 samples x n_probes. That shipped once (commit 4ec3493) and was
# reverted; tests/test_load_path.R guards it. The genuine saving from that
# commit was the rm(beta)+gc() below, which is kept.
cat("\nreading beta matrix ...\n")
beta <- data.table::fread(beta_path, data.table = FALSE, check.names = FALSE)
if (anyDuplicated(beta[[1]])) stop("Duplicate probe IDs")
x <- t(as.matrix(beta[, -1, drop = FALSE])); storage.mode(x) <- "double"; colnames(x) <- beta[[1]]
rm(beta); invisible(gc(verbose = FALSE))
if (any(is.finite(x) & (x < 0 | x > 1))) stop("Expected beta values in [0,1]")
if (!nrow(x) || !ncol(x)) stop("Matrix loaded with a zero dimension - check the load orientation")
cat(sprintf("  %d samples x %d probes\n", nrow(x), ncol(x)))

# --- Load and validate the sample metadata ---------------------------------
meta <- read.delim(meta_path, check.names = FALSE, stringsAsFactors = FALSE)
if (anyDuplicated(meta$sample_id) || anyDuplicated(meta$patient_id)) stop("Duplicate patient/specimen")
if (any(!rownames(x) %in% meta$sample_id)) stop("Missing sample metadata")
# EVERY subsequent meta$... lookup depends on this alignment holding.
meta <- meta[match(rownames(x), meta$sample_id), ]
if (any(!is.finite(meta$HRDsum)) || any(meta$HRDsum < 0)) stop("Invalid target")
# Label integrity: catches HRD_LOH/LST/TAI column mis-mapping, which would
# otherwise train the model against a scrambled target invisibly.
if (!all(abs(meta$HRDsum - meta$HRD_LOH - meta$LST - meta$TAI) < 1e-6)) stop("Label sum mismatch")
if (!"quality_annotation" %in% names(meta)) stop("Missing quality_annotation column")

# =============================================================================
# THE LOCK - guard 1 of 3
# =============================================================================
# The hardcoded GBM/LGG list and build_master.py's `partition` column are two
# sources of truth for the SAME lock. Compare them, in both directions, before
# either is used. This changes nothing about which samples are locked; it
# refuses to run when the two definitions disagree, which is the only situation
# in which the hardcoded list could silently be wrong.
assert_partition_matches_cns(meta$cancer_type, meta$partition,
                             context = basename(meta_path))
cns <- meta$cancer_type %in% CNS_LOCKED_TYPES
cat(sprintf("\nlock       : %d locked_CNS rows (%s), %d development rows\n",
            sum(cns), paste(sort(unique(meta$cancer_type[cns])), collapse = "+"), sum(!cns)))

# Anti-fixture guards: refuse to make transfer claims from a toy cohort.
dev_types <- sort(unique(meta$cancer_type[!cns]))
if (length(dev_types) < 3) stop("Need >=3 development cancers")
if (sum(!cns) < 50) stop("Need >=50 development patients for this protocol")

# =============================================================================
# CALIBRATION RESERVATION
# =============================================================================
# Stratified 20% per cancer type, so the calibration set mirrors the development
# cohort's tissue composition rather than being dominated by BRCA (n=743). The
# seed makes the split reproducible; it is the same seed, the same fraction and
# the same function body as train_baseline.R's Phase 2, so this reservation is
# byte-identical to the one that script would have produced.
#
# sample.int(length(ii), k) NOT sample(ii, k): with length(ii) == 1 the latter
# is interpreted as sample(1:ii, k) and returns a random integer in 1..ii rather
# than ii itself - an arbitrary row, possibly a locked_CNS one. It does not fire
# on the current cohort (smallest development group is OV, n=10) but it would on
# any re-partition producing a singleton group.
pick_calibration <- function(ii) {
  k <- max(1L, floor(length(ii) * CAL_FRACTION))
  ii[sample.int(length(ii), k)]
}
set.seed(SEED)
dev <- which(!cns)
cal <- unlist(lapply(split(dev, meta$cancer_type[dev]), pick_calibration), use.names = FALSE)
cal <- sort(cal)
stopifnot(!anyDuplicated(cal), all(cal %in% dev))
tr <- setdiff(dev, cal)
if (!length(tr)) stop("Calibration reservation consumed every development sample")

# =============================================================================
# THE LOCK - guard 2 of 3: hard-stop if any CNS row would enter the fit
# =============================================================================
# assert_partition_matches_cns() above proves the lock DEFINITION is unambiguous.
# This proves the lock was actually APPLIED to the row sets about to be used.
# Both are needed: the first can pass while a downstream indexing bug still
# leaks a locked row into `tr`.
locked_rows <- which(cns)
leak_fit <- intersect(tr, locked_rows)
leak_cal <- intersect(cal, locked_rows)
if (length(leak_fit) || length(leak_cal)) {
  stop(sprintf(paste0("LOCK VIOLATION: %d locked_CNS row(s) in the training set and %d in the ",
                      "calibration set (types: %s). Refusing to fit. The locked partition must ",
                      "never enter fitting, preprocessing, feature selection, hyperparameter ",
                      "selection or conformal calibration."),
               length(leak_fit), length(leak_cal),
               paste(sort(unique(meta$cancer_type[c(leak_fit, leak_cal)])), collapse = ", ")),
       call. = FALSE)
}
# The three roles must partition the cohort exactly - no row unaccounted for,
# no row counted twice. This is the statement final_partitions.tsv encodes.
stopifnot(length(tr) + length(cal) + length(locked_rows) == nrow(meta),
          !length(intersect(tr, cal)),
          setequal(c(tr, cal, locked_rows), seq_len(nrow(meta))),
          !any(cns[tr]), !any(cns[cal]),
          setequal(sort(unique(meta$cancer_type[tr])), dev_types))
cat(sprintf("lock guard : PASS - 0 locked rows in %d training + %d calibration rows\n",
            length(tr), length(cal)))
cat(sprintf("split      : train n=%d  calibration n=%d  across %d cancer types\n",
            length(tr), length(cal), length(dev_types)))

# =============================================================================
# THE FIT
# =============================================================================
# fit_en() runs its own inner CV. inner_folds() assigns one fold per TRAINING
# cancer type, so hyperparameters are selected for CROSS-TISSUE generalisation -
# the repo's rule, not a random-sample CV. Preprocessing (missingness filter,
# medians, centre/scale, and the 5,000-probe pooled-variance filter) is
# re-learned inside every inner fold from that fold's own training rows.
#
# tf$forward is the IDENTITY, so this fits on raw HRDsum exactly as run 01 did.
# The clip lives on the inverse and is applied at prediction time only.
cat("\nfitting final model (fit_en, feature_rank='", FEATURE_RANK, "') ...\n", sep = "")
t_fit <- Sys.time()
b <- fit_en(x[tr, , drop = FALSE],
            tf$forward(meta$HRDsum[tr]),
            meta$patient_id[tr],
            meta$cancer_type[tr],
            seed = SEED,
            feature_rank = FEATURE_RANK)
cat(sprintf("fit_en completed in %.1f min\n", as.numeric(difftime(Sys.time(), t_fit, units = "mins"))))
cat(sprintf("  alpha=%.2g lambda=%.5g  (lambda at path boundary: %s)\n",
            b$alpha, b$lambda, b$lambda_boundary_side))
stopifnot(identical(b$feature_rank, FEATURE_RANK), identical(b$seed, SEED))

# =============================================================================
# SPLIT-CONFORMAL INTERVAL on the reserved calibration set
# =============================================================================
# The calibration rows were never seen by the fit, the preprocessing or the
# hyperparameter search, which is what makes the interval's marginal coverage
# guarantee real. No refit happens here.
#
# Residuals are taken against the CLIPPED prediction, because predict_frozen.R
# centres [pred-q, pred+q] on the clipped column. Deriving q from the raw column
# and applying it to the clipped one would be a half-width for a different point
# estimate.
cp <- predict_en(b, x[cal, , drop = FALSE])
cal_pred <- tf$inverse(cp$predicted_reference_HRDsum)
b$interval_q <- conformal_q(meta$HRDsum[cal], cal_pred, coverage = COVERAGE)
cat(sprintf("\nconformal  : q=%.4f at %.0f%% coverage from n=%d calibration samples\n",
            b$interval_q, 100 * COVERAGE, length(cal)))
if (!is.finite(b$interval_q)) {
  warning("Conformal q is Inf: too few calibration points for the requested coverage. ",
          "predict_frozen.R will report interval_status='unavailable_insufficient_calibration'.",
          call. = FALSE)
}
cal_mae <- mean(abs(meta$HRDsum[cal] - cal_pred))

# =============================================================================
# BUNDLE THE ARTEFACT
# =============================================================================
# Self-contained by design: anyone loading this file must be able to answer
# "what was it fitted on, under what rule, and what must be done to its output"
# without this repository.
b$target_transform  <- TARGET_TRANSFORM       # <- the field predict_frozen.R:299 reads
b$calibration_ids   <- meta$patient_id[cal]
b$training_ids      <- meta$patient_id[tr]
b$training_cancer_types <- dev_types
b$probe_provenance  <- matrix_pr              # predict_frozen.R:223 compares allowlists against this
b$conformal_coverage <- COVERAGE
b$interval_note     <- "Split-conformal residual interval; no pediatric/domain-shift coverage guarantee"
b$frozen_at         <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
b$locked_partition  <- list(types = CNS_LOCKED_TYPES, n = sum(cns),
                            note = "Never fitted, preprocessed, selected on, tuned on or calibrated on.")
b$is_fixture        <- isTRUE(gate$fatal)
b$provenance_class  <- gate$class

rds_path <- file.path(out, "frozen_nonCNS.rds")
saveRDS(b, rds_path)
cat("\nwrote", rds_path, "\n")

# --- Coefficients -----------------------------------------------------------
# Non-zero only, sorted by |coefficient| descending. This count, not
# `selected_features`, is the real sparsity measure: selected_features counts
# probes ENTERING glmnet, and the elastic net then shrinks most of them to zero.
# The intercept is kept (it is the model's offset, and dropping it would make
# the file unable to reproduce a prediction) but tagged so it is never mistaken
# for a probe when the table is ranked.
cc <- as.matrix(coef(b$model, s = b$lambda))
coef_df <- data.frame(probe = rownames(cc),
                      coefficient = as.numeric(cc[, 1]),
                      stringsAsFactors = FALSE)
coef_df$term_type <- ifelse(coef_df$probe == "(Intercept)", "intercept", "probe")
n_features_in <- length(b$preprocess$features)
n_nonzero <- sum(coef_df$coefficient != 0 & coef_df$term_type == "probe")
coef_df <- coef_df[coef_df$coefficient != 0, , drop = FALSE]
coef_df <- coef_df[order(coef_df$term_type != "intercept",
                         -abs(coef_df$coefficient)), , drop = FALSE]
write.table(coef_df, file.path(out, "frozen_coefficients.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Full role assignment for every sample ---------------------------------
# The audit trail proving the locked set was never fitted. One row per specimen
# in the master table, including the 642 locked ones.
roles <- rep("training", nrow(meta))
roles[cal] <- "calibration"
roles[cns] <- "locked_CNS"
write.table(data.frame(sample_id = meta$sample_id,
                       patient_id = meta$patient_id,
                       cancer_type = meta$cancer_type,
                       partition = meta$partition,
                       role = roles,
                       stringsAsFactors = FALSE),
            file.path(out, "final_partitions.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Manifest ---------------------------------------------------------------
# sha256 of the .rds is computed AFTER writing, so it certifies the bytes on
# disk rather than an in-memory object that could still have been mutated.
file_sha256 <- function(path) {
  if (!file.exists(path)) return("file_absent")
  if (requireNamespace("digest", quietly = TRUE))
    return(digest::digest(file = path, algo = "sha256"))
  o <- suppressWarnings(try(system2("sha256sum", shQuote(path), stdout = TRUE, stderr = FALSE), silent = TRUE))
  if (inherits(o, "try-error") || !length(o)) return("sha256_unavailable")
  sub("\\s.*$", "", o[1])
}
git_field <- function(a) {
  o <- suppressWarnings(try(system2("git", a, stdout = TRUE, stderr = FALSE), silent = TRUE))
  if (inherits(o, "try-error") || !length(o)) return(NA_character_)
  paste(o, collapse = " ")
}
rds_sha <- file_sha256(rds_path)

man <- data.frame(
  field = c("artifact", "artifact_sha256", "is_fixture", "matrix_provenance_class",
            "git_commit", "git_dirty", "frozen_at_utc", "elapsed_min",
            "seed", "feature_rank", "target_transform", "target_transform_forward",
            "target_transform_inverse", "alpha", "lambda", "lambda_at_boundary",
            "lambda_boundary_side", "n_features_entering_glmnet", "n_nonzero_coefficients",
            "n_training_samples", "n_calibration_samples", "n_locked_CNS_excluded",
            "n_training_cancer_types", "training_cancer_types", "calibration_fraction",
            "conformal_coverage", "conformal_q", "calibration_MAE", "ood_cut",
            "training_mean_HRDsum", "n_inner_folds", "beta_matrix", "beta_matrix_sha256",
            "master_table", "master_table_sha256", "probe_allowlist_sha256",
            "R_version", "glmnet_version", "host", "user"),
  value = c(rds_path, rds_sha, as.character(b$is_fixture), gate$class,
            git_field(c("rev-parse", "HEAD")),
            as.character(nzchar(git_field(c("status", "--porcelain")))),
            format(Sys.time(), "%Y-%m-%dT%H:%M:%S", tz = "UTC"),
            sprintf("%.1f", as.numeric(difftime(Sys.time(), t_start, units = "mins"))),
            as.character(SEED), FEATURE_RANK, TARGET_TRANSFORM, "identity", "pmax(p, 0)",
            format(b$alpha), format(b$lambda, digits = 10),
            as.character(b$lambda_at_boundary), as.character(b$lambda_boundary_side),
            as.character(n_features_in), as.character(n_nonzero),
            as.character(length(tr)), as.character(length(cal)), as.character(sum(cns)),
            as.character(length(dev_types)), paste(dev_types, collapse = ","),
            format(CAL_FRACTION), format(COVERAGE), format(b$interval_q, digits = 10),
            sprintf("%.4f", cal_mae), format(b$ood_cut, digits = 10),
            format(b$training_mean, digits = 10),
            as.character(length(unique(b$inner_folds$fold))),
            beta_path,
            if (!is.null(matrix_pr$sha256)) as.character(matrix_pr$sha256) else "not_recorded",
            meta_path, file_sha256(meta_path),
            if (!is.null(matrix_pr$probe_allowlist_sha256)) as.character(matrix_pr$probe_allowlist_sha256) else "not_recorded",
            R.version.string, as.character(utils::packageVersion("glmnet")),
            Sys.info()[["nodename"]], Sys.info()[["user"]]),
  stringsAsFactors = FALSE)
write.table(man, file.path(out, "freeze_manifest.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n=============================================================\n")
cat("FROZEN MODEL SUMMARY\n")
cat("=============================================================\n")
cat(sprintf("  artifact             : %s\n", rds_path))
cat(sprintf("  sha256               : %s\n", rds_sha))
cat(sprintf("  git commit           : %s\n", man$value[man$field == "git_commit"]))
cat(sprintf("  seed                 : %d\n", SEED))
cat(sprintf("  feature_rank         : %s\n", FEATURE_RANK))
cat(sprintf("  target_transform     : %s   (fit: identity, score: pmax(p,0) = C2 ADOPTED)\n", TARGET_TRANSFORM))
cat(sprintf("  alpha / lambda       : %.3g / %.6g  (boundary: %s)\n", b$alpha, b$lambda, b$lambda_boundary_side))
cat(sprintf("  features -> glmnet   : %d\n", n_features_in))
cat(sprintf("  non-zero coefficients: %d  (%.1f%% of features entering)\n",
            n_nonzero, 100 * n_nonzero / max(1L, n_features_in)))
cat(sprintf("  training samples     : %d  across %d cancer types\n", length(tr), length(dev_types)))
cat(sprintf("  calibration samples  : %d  (stratified %.0f%% per cancer type)\n", length(cal), 100 * CAL_FRACTION))
cat(sprintf("  inner folds          : %d  (one per training cancer type)\n", length(unique(b$inner_folds$fold))))
cat(sprintf("  conformal q (%.0f%%)    : %.4f   calibration MAE %.4f\n", 100 * COVERAGE, b$interval_q, cal_mae))
cat(sprintf("  OOD cutoff           : %.6g\n", b$ood_cut))
cat(sprintf("  locked_CNS EXCLUDED  : %d rows (%s) - never read\n",
            sum(cns), paste(CNS_LOCKED_TYPES, collapse = "+")))
cat(sprintf("  elapsed              : %.1f min\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
if (isTRUE(b$is_fixture)) {
  cat("\n  *** FIXTURE RUN - NOT THE SHIPPING ARTEFACT ***\n")
}
cat("\nNOTE: this script produces no development metrics. The performance\n")
cat("numbers for this architecture are run 01's (results/loco_run01), from\n")
cat("DIFFERENT fits. Do not present macro_metrics.txt as describing this .rds.\n")
cat("\nNext: score the locked cohort ONCE via scripts/predict_frozen.R.\n")
cat("finished   :", format(Sys.time()), "\n")

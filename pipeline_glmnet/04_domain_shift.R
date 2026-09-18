#!/usr/bin/env Rscript
# =============================================================================
# 04_domain_shift.R  —  KIDS26 Team 12
# Quantifies the cross-tissue OFFSET problem (issue C1) using LOCO predictions
# that already exist. No refitting, runs in seconds.
#
#   Rscript 04_domain_shift.R --run_dir=results/full_none
#
# Produces:
#   domain_shift_per_cancer.tsv  per held-out cancer: rank quality (tissue-free),
#                                mean bias, and MAE before/after offset correction
#   calibration_curve.tsv        MAE vs number of labelled samples k used to fix
#                                the intercept in the new tissue (k = 0,3,5,10,20,
#                                repeated `--reps` times with random draws)
#   Console summary of both.
#
# Why: LOCO shows within-cancer ranking transfers but the level does not. This
# asks how many labelled samples from a NEW tissue are needed to anchor it.
# =============================================================================
suppressPackageStartupMessages(library(data.table))
arg <- function(name, default = NULL) {
  a <- commandArgs(trailingOnly = TRUE)
  hit <- grep(paste0("^--", name, "="), a, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
run_dir <- arg("run_dir", "results/full_none")
ks      <- as.integer(strsplit(arg("k", "0,3,5,10,20"), ",")[[1]])
reps    <- as.integer(arg("reps", "200"))
min_n   <- as.integer(arg("min_n", "20"))
seed    <- as.integer(arg("seed", "260910"))
hrd_cut <- as.numeric(arg("hrd_threshold", "42"))
set.seed(seed)

p <- fread(file.path(run_dir, "loco_predictions.tsv"))
setnames(p, "predicted_reference_HRDsum", "pred", skip_absent = TRUE)
stopifnot(all(c("cancer_type", "actual", "pred") %in% names(p)))

auc_rank <- function(label, score) {
  label <- as.logical(label); n1 <- sum(label); n0 <- sum(!label)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score); (sum(r[label]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# ---- per-cancer summary: ranking (offset-free) vs level -------------------------
per <- p[, {
  rho  <- if (.N >= 3 && sd(actual) > 0 && sd(pred) > 0) cor(actual, pred, method = "spearman") else NA_real_
  bias <- mean(pred - actual)
  .(n = .N, spearman = round(rho, 3),
    auc_ge_cut = round(auc_rank(actual >= hrd_cut, pred), 3),
    mean_actual = round(mean(actual), 1), mean_pred = round(mean(pred), 1),
    bias = round(bias, 1),
    MAE_raw = round(mean(abs(pred - actual)), 2),
    MAE_oracle_offset = round(mean(abs(pred - bias - actual)), 2),   # best possible shift
    MAE_null = round(mean(abs(actual - mean(actual))), 2))
}, by = cancer_type][order(-n)]

# ---- calibration curve: fix the intercept with k labelled samples ---------------
cal <- rbindlist(lapply(split(p, p$cancer_type), function(d) {
  if (nrow(d) < min_n) return(NULL)
  rbindlist(lapply(ks, function(k) {
    if (k == 0) {
      return(data.table(cancer_type = d$cancer_type[1], k = 0L,
                        MAE = mean(abs(d$pred - d$actual)),
                        MAE_sd = NA_real_, n_eval = nrow(d)))
    }
    v <- vapply(seq_len(reps), function(i) {
      idx <- sample.int(nrow(d), k)
      off <- mean(d$pred[idx] - d$actual[idx])       # intercept learned on k labels
      te  <- setdiff(seq_len(nrow(d)), idx)
      mean(abs(d$pred[te] - off - d$actual[te]))
    }, 0.0)
    data.table(cancer_type = d$cancer_type[1], k = k, MAE = mean(v),
               MAE_sd = sd(v), n_eval = nrow(d) - k)
  }))
}))

summary_k <- cal[, .(median_MAE = round(median(MAE), 2),
                     mean_MAE = round(mean(MAE), 2),
                     n_cancers = .N), by = k][order(k)]
null_macro <- per[n >= min_n, round(mean(MAE_null), 2)]

fwrite(per, file.path(run_dir, "domain_shift_per_cancer.tsv"), sep = "\t")
fwrite(cal, file.path(run_dir, "calibration_curve.tsv"), sep = "\t")

cat("\n== Per held-out cancer: ranking transfers, level does not ==\n")
print(per)
cat("\n== Calibration curve: MAE vs labelled samples used for the intercept ==\n")
print(summary_k)
cat("\nMacro null MAE (predict that cancer's own mean):", null_macro, "\n")
cat("Wrote domain_shift_per_cancer.tsv and calibration_curve.tsv to", run_dir, "\n")

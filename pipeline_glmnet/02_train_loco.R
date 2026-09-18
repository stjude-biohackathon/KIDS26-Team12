#!/usr/bin/env Rscript
# =============================================================================
# 02_train_loco.R  —  KIDS26 Team 12
# Memory-safe, parallel, resumable replacement for scripts/train_baseline.R.
# Same design and outputs as train_baseline.R:
#   * leave-one-cancer-out (LOCO) over non-CNS development cancers
#   * per fold: <=5000 features selected on training rows, cv.glmnet alpha x lambda tuning
#   * frozen non-CNS model: 20% per-cancer calibration split (seed 260910),
#     split-conformal 95% interval, NO refit on calibration
#   * CNS (GBM/LGG) is never loaded or evaluated here
#
# Usage:
#   Rscript 02_train_loco.R --features=features --out=results/run01 --cores=8
# Options: --inner=cancer|random  --control=none|permute|permute_within_cancer
#          --alphas=0.1,0.25,0.5,0.75,1 --inner_folds=5 --max_features=5000
#          --target=raw|within_z (C1c relative-target; score with Spearman/AUC)
#          --lambda_rule=1se|min --transform=beta|m --only=ACC,KICH (smoke test)
# Re-running the same command skips folds already saved (safe after a crash).
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(parallel) })
arg <- function(name, default = NULL) {
  a <- commandArgs(trailingOnly = TRUE)
  hit <- grep(paste0("^--", name, "="), a, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(arg("model_r", file.path(script_dir, "R", "model.R")))

feat_dir     <- arg("features", "features")
out          <- arg("out", "results/run01")
alphas       <- as.numeric(strsplit(arg("alphas", "0.1,0.25,0.5,0.75,1"), ",")[[1]])
n_inner      <- as.integer(arg("inner_folds", "5"))
max_features <- as.integer(arg("max_features", "5000"))
max_missing  <- as.numeric(arg("max_missing", "0.05"))
lambda_rule  <- arg("lambda_rule", "1se")
transform    <- arg("transform", "beta")
seed         <- as.integer(arg("seed", "260910"))
coverage     <- as.numeric(arg("coverage", "0.95"))
cal_frac     <- as.numeric(arg("calibration_fraction", "0.2"))
cores        <- as.integer(arg("cores", as.character(max(1L, min(8L, detectCores() - 1L)))))
only         <- arg("only", NULL)
inner        <- arg("inner", "cancer")      # cancer (nested LOCO, team design) | random
control      <- arg("control", "none")      # none | permute | permute_within_cancer
hrd_cut      <- as.numeric(arg("hrd_threshold", "42"))
target       <- arg("target", "raw")        # raw HRDsum | within_z (C1c relative-target)
stopifnot(target %in% c("raw", "within_z"), inner %in% c("cancer", "random"), control %in% c("none", "permute", "permute_within_cancer"),
          transform %in% c("beta", "m"), lambda_rule %in% c("1se", "min"), all(alphas >= 0 & alphas <= 1))
dir.create(file.path(out, "folds"), recursive = TRUE, showWarnings = FALSE)
msg <- function(...) { cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = ""); flush.console() }
t_start <- Sys.time()

# ---- load prepared features ------------------------------------------------------
X    <- readRDS(file.path(feat_dir, "X_dev.rds"))
meta <- as.data.frame(fread(file.path(feat_dir, "meta_dev.tsv"), sep = "\t"))
stopifnot(identical(rownames(X), as.character(meta$sample_id)))
if (any(meta$cancer_type %in% c("GBM", "LGG"))) stop("CNS samples present in development features")
if (transform == "m") X <- to_m(X)
y <- meta$HRDsum
# Negative controls: models are TRAINED on shuffled labels and EVALUATED on true
# labels. permute = pure null; permute_within_cancer keeps each cancer's HRDsum
# distribution, so it is the tissue-identity control: any score it achieves is
# what lineage alone buys. The real model must clearly beat both.
y_train <- y
if (control != "none") {
  set.seed(seed + 1L)
  y_train <- if (control == "permute") sample(y) else
    ave(y, meta$cancer_type, FUN = function(v) if (length(v) > 1) sample(v) else v)
  msg("NEGATIVE CONTROL: training labels shuffled (", control, ")")
}
# C1c relative-target: train on HRDsum standardised WITHIN each training cancer, so
# the between-tissue level is removed from the target and only within-tissue ordering
# is learnable. Predictions are z-scores: judge with Spearman/AUC, NOT MAE/R2.
if (target == "within_z") {
  y_train <- ave(y_train, meta$cancer_type,
                 FUN = function(v) if (length(v) > 1 && sd(v) > 0) (v - mean(v)) / sd(v) else rep(0, length(v)))
  msg("RELATIVE TARGET (C1c): within-cancer standardised HRDsum; MAE/R2 are not on the HRDsum scale")
}
msg("Loaded ", nrow(X), " samples x ", ncol(X), " prefiltered probes; transform=", transform, "; target=", target,
    "; cores=", cores, "; alphas=", paste(alphas, collapse = ","), "; inner=", inner, "; control=", control)

# Global column sums: per-fold training variance = global minus held-out rows,
# so we never copy the full training matrix just to rank features.
S1 <- colSums(X, na.rm = TRUE); S2 <- colSums(X^2, na.rm = TRUE); NN <- colSums(!is.na(X))
select_features <- function(test_rows) {
  xt <- X[test_rows, , drop = FALSE]
  n  <- NN - colSums(!is.na(xt))
  s1 <- S1 - colSums(xt, na.rm = TRUE); s2 <- S2 - colSums(xt^2, na.rm = TRUE)
  v  <- s2 / n - (s1 / n)^2
  miss <- 1 - n / (nrow(X) - length(test_rows))
  ok <- which(miss <= max_missing & is.finite(v) & v > 1e-12)
  ok[order(v[ok], decreasing = TRUE)][seq_len(min(max_features, length(ok)))]
}

fit_on <- function(train_rows) {
  keep <- select_features(setdiff(seq_len(nrow(X)), train_rows))
  b <- fit_en(X[train_rows, keep, drop = FALSE], y_train[train_rows], meta$patient_id[train_rows],
              meta$cancer_type[train_rows], alphas = alphas, n_inner = n_inner, inner = inner,
              max_features = max_features, max_missing = max_missing, lambda_rule = lambda_rule,
              transform_label = "none", seed = seed)
  b$preprocess$transform <- transform      # predict_en() receives raw betas
  b
}

cancers <- sort(unique(meta$cancer_type))
if (!is.null(only)) cancers <- intersect(cancers, strsplit(only, ",")[[1]])
if (length(unique(meta$cancer_type)) < 3 || nrow(meta) < 50) stop("Need >=3 cancers and >=50 patients")
tasks <- c(cancers, "__frozen__")

run_task <- function(task) {
  f_out <- file.path(out, "folds", paste0(if (task == "__frozen__") "frozen" else paste0("loco_", task), ".rds"))
  if (file.exists(f_out)) return(paste("skipped (exists):", task))
  t0 <- Sys.time()
  if (task == "__frozen__") {
    set.seed(seed)
    idx <- seq_len(nrow(meta))
    cal <- unlist(lapply(split(idx, meta$cancer_type), function(ii)
      if (length(ii) == 1) ii else sample(ii, max(1L, floor(length(ii) * cal_frac)))))
    tr <- setdiff(idx, cal)
    b <- fit_on(tr)
    cp <- predict_en(b, X_raw_rows(cal))
    b$interval_q <- conformal_q(y[cal], cp$predicted_reference_HRDsum, coverage)
    b$coverage <- coverage
    b$calibration_ids <- meta$patient_id[cal]; b$training_ids <- meta$patient_id[tr]
    b$interval_note <- "Split-conformal residual interval; no pediatric/domain-shift coverage guarantee"
    b$calibration_metrics <- metrics(y[cal], cp$predicted_reference_HRDsum)
    b$role <- ifelse(idx %in% cal, "calibration", "training")
  } else {
    te <- which(meta$cancer_type == task); tr <- which(meta$cancer_type != task)
    b <- fit_on(tr)
    p <- predict_en(b, X_raw_rows(te))
    p$actual <- y[te]; p$cancer_type <- task; p$patient_id <- meta$patient_id[te]
    p$null_prediction <- b$training_mean; p$split <- "development_LOCO"
    b$predictions <- p
    mm <- metrics(p$actual, p$predicted_reference_HRDsum)
    mm$cancer_type <- task; mm$train_n <- length(tr); mm$selected_features <- length(b$preprocess$features)
    mm$nonzero_coefficients <- sum(as.matrix(coef(b$model, s = b$lambda))[-1, 1] != 0)
    mm$alpha <- b$alpha; mm$lambda <- b$lambda
    mm$null_MAE <- mean(abs(p$actual - b$training_mean)); mm$abstention_rate <- mean(!p$reportable)
    b$fold_metrics <- mm
  }
  b$minutes <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  saveRDS(b, paste0(f_out, ".tmp")); file.rename(paste0(f_out, ".tmp"), f_out)   # atomic save
  paste0("done: ", task, " (", b$minutes, " min, alpha=", b$alpha, ")")
}
# predict_en expects raw betas; undo the M transform only for prediction input
X_raw_rows <- function(rows) {
  x <- X[rows, , drop = FALSE]
  if (transform == "m") { x <- 2^x / (1 + 2^x) }
  x
}

msg("Running ", length(tasks), " tasks (", length(cancers), " LOCO folds + frozen model)")
res <- mclapply(tasks, function(t) tryCatch(run_task(t), error = function(e) paste("FAILED:", t, "-", conditionMessage(e))),
                mc.cores = cores, mc.preschedule = FALSE)
for (r in res) msg(r)
failed <- grep("^FAILED", unlist(res), value = TRUE)

# ---- collect outputs (same file names as train_baseline.R) ------------------------
fold_files <- file.path(out, "folds", paste0("loco_", cancers, ".rds"))
have <- file.exists(fold_files)
folds <- lapply(fold_files[have], readRDS)
if (length(folds)) {
  pred <- rbindlist(lapply(folds, `[[`, "predictions"))
  mt   <- rbindlist(lapply(folds, `[[`, "fold_metrics"))
  tune <- rbindlist(lapply(folds, function(b) cbind(cancer_type = b$fold_metrics$cancer_type,
                                                    selected = b$tuning$alpha == b$alpha, b$tuning)))
  fwrite(pred, file.path(out, "loco_predictions.tsv"), sep = "\t")
  fwrite(mt,   file.path(out, "loco_metrics.tsv"), sep = "\t")
  fwrite(tune, file.path(out, "tuning_by_fold.tsv"), sep = "\t")
  for (b in folds) fwrite(b$inner_folds, file.path(out, paste0("inner_folds_", b$fold_metrics$cancer_type, ".tsv")), sep = "\t")
  pooled <- metrics(pred$actual, pred$predicted_reference_HRDsum)
  auc42 <- auc_rank(pred$actual >= hrd_cut, pred$predicted_reference_HRDsum)
  within <- pred[, .(rho = if (.N >= 20 && sd(actual) > 0 && sd(predicted_reference_HRDsum) > 0) cor(actual, predicted_reference_HRDsum, method = "spearman") else NA_real_), by = cancer_type]
  writeLines(c(paste("macro_MAE", mean(mt$MAE)),
               paste("macro_null_MAE", mean(mt$null_MAE)),
               paste("pooled_MAE", pooled$MAE), paste("pooled_spearman", pooled$spearman),
               paste("median_within_cancer_spearman_n20", median(within$rho, na.rm = TRUE)),
               paste0("pooled_AUC_HRDsum_ge", hrd_cut, " ", auc42),
               paste("control", control), paste("inner_folds", inner), paste("transform", transform),
               paste("target", target),
               paste("folds_completed", sum(have), "of", length(cancers))),
             file.path(out, "macro_metrics.txt"))
  msg("macro MAE = ", round(mean(mt$MAE), 2), " vs null ", round(mean(mt$null_MAE), 2),
      "; pooled Spearman = ", round(pooled$spearman, 3), "; AUC(>=", hrd_cut, ") = ", round(auc42, 3))
}
frozen_file <- file.path(out, "folds", "frozen.rds")
if (file.exists(frozen_file)) {
  b <- readRDS(frozen_file)
  b$probe_provenance <- paste(readLines(file.path(feat_dir, "beta.provenance.json"), warn = FALSE), collapse = "\n")
  saveRDS(b, file.path(out, "frozen_nonCNS.rds"))
  cc <- as.matrix(coef(b$model, s = b$lambda))
  fwrite(data.table(feature = rownames(cc), coefficient = cc[, 1]), file.path(out, "frozen_coefficients.tsv"), sep = "\t")
  fwrite(data.table(sample_id = meta$sample_id, patient_id = meta$patient_id,
                    cancer_type = meta$cancer_type, role = b$role),
         file.path(out, "final_partitions.tsv"), sep = "\t")
  msg("Frozen model: alpha=", b$alpha, " lambda=", signif(b$lambda, 4), " interval +/-", round(b$interval_q, 2))
}
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
msg("Finished in ", round(as.numeric(difftime(Sys.time(), t_start, units = "hours")), 2), " h")
if (length(failed)) { cat(failed, sep = "\n"); quit(status = 1) }

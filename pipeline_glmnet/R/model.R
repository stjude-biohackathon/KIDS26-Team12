# =============================================================================
# R/model.R  —  KIDS26 Team 12 elastic-net core (replacement for the missing file)
#
# Defines the four functions train_baseline.R expects:
#   fit_en(), predict_en(), metrics(), conformal_q()
# Follows config/analysis_protocol.json: continuous reference_HRDsum, elastic net,
# <= 5000 features, <= 5% missing per probe, seed 260910, 95% conformal coverage.
# Dependencies: glmnet (and base R). Works on R >= 4.1 (tested 4.3, 4.6).
# =============================================================================
suppressPackageStartupMessages(library(glmnet))

# Beta -> M-value (optional transform; bounded betas are heteroscedastic near 0/1)
to_m <- function(b, eps = 1e-3) { b <- pmin(pmax(b, eps), 1 - eps); log2(b / (1 - b)) }

# Fold ids stratified by a grouping vector (here: cancer type)
strat_folds <- function(strata, k, seed) {
  set.seed(seed)
  f <- integer(length(strata))
  for (s in unique(strata)) {
    i <- which(strata == s)
    f[i] <- sample(rep_len(sample.int(k), length(i)))
  }
  f
}

# Column variance with NAs, without copying beyond one squared matrix
col_var <- function(x) {
  mu <- colMeans(x, na.rm = TRUE)
  colMeans(x^2, na.rm = TRUE) - mu^2
}

# -----------------------------------------------------------------------------
# fit_en: feature selection + imputation learned on TRAINING rows only, then
# automated alpha x lambda tuning with cv.glmnet on one fixed inner foldid.
# glmnet standardize=TRUE is fold-safe: each internal fit standardises using only
# its own training rows, and predictions use coefficients on the original scale.
# The medians/variance ranking below use all outer-training rows (label-free);
# this can only affect hyperparameter choice, never the outer held-out cancer.
# -----------------------------------------------------------------------------
fit_en <- function(x, y, patient_id, cancer_type,
                   alphas = c(0.1, 0.25, 0.5, 0.75, 1), n_inner = 5L, inner = c("cancer", "random"),
                   max_features = 5000L, max_missing = 0.05,
                   lambda_rule = c("1se", "min"), type_measure = "mse",
                   transform_label = "beta", seed = 260910L, verbose = TRUE) {
  lambda_rule <- match.arg(lambda_rule); inner <- match.arg(inner)
  stopifnot(nrow(x) == length(y), length(y) == length(patient_id),
            length(y) == length(cancer_type), !is.null(colnames(x)), all(is.finite(y)))
  if (anyDuplicated(patient_id)) stop("fit_en: duplicate patients in training data")

  # 1) feature filter (label-free): missingness, then top variance
  na_frac <- colMeans(is.na(x))
  v <- col_var(x)
  ok <- which(na_frac <= max_missing & is.finite(v) & v > 0)
  if (!length(ok)) stop("fit_en: no features pass the missingness/variance filter")
  ok <- ok[order(v[ok], decreasing = TRUE)][seq_len(min(max_features, length(ok)))]
  features <- colnames(x)[ok]
  xs <- x[, ok, drop = FALSE]

  # 2) impute with training medians; record training ranges for OOD flags
  med <- apply(xs, 2, median, na.rm = TRUE)
  if (anyNA(xs)) { idx <- which(is.na(xs), arr.ind = TRUE); xs[idx] <- med[idx[, 2]] }
  q_lo <- apply(xs, 2, quantile, probs = 0.005, names = FALSE)
  q_hi <- apply(xs, 2, quantile, probs = 0.995, names = FALSE)

  # 3) automated tuning: same inner folds for every alpha
  # Inner folds. Default (team design): leave-one-cancer-out INSIDE the training
  # set, so alpha/lambda are chosen for cross-tissue generalisation. Deterministic,
  # so no repeated-CV needed for lambda stability. Fallback: stratified random.
  if (inner == "cancer" && length(unique(cancer_type)) >= 3) {
    foldid <- match(cancer_type, sort(unique(cancer_type)))
  } else {
    foldid <- strat_folds(cancer_type, n_inner, seed)
  }
  fits <- vector("list", length(alphas)); tune <- vector("list", length(alphas))
  for (j in seq_along(alphas)) {
    t0 <- Sys.time()
    cv <- cv.glmnet(xs, y, family = "gaussian", alpha = alphas[j], foldid = foldid,
                    type.measure = type_measure, standardize = TRUE)
    i_min <- which(cv$lambda == cv$lambda.min); i_1se <- which(cv$lambda == cv$lambda.1se)
    tune[[j]] <- data.frame(alpha = alphas[j], lambda_min = cv$lambda.min, lambda_1se = cv$lambda.1se,
                            cv_error_min = cv$cvm[i_min], cv_error_1se = cv$cvm[i_1se],
                            nonzero_min = cv$nzero[i_min], nonzero_1se = cv$nzero[i_1se],
                            lambda_min_at_path_end = cv$lambda.min == min(cv$lambda),
                            minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2))
    fits[[j]] <- cv
    if (verbose) message(sprintf("    alpha=%.2f lambda.min=%.4g (cv=%.4g, nz=%d) lambda.1se=%.4g (nz=%d)",
                                 alphas[j], cv$lambda.min, cv$cvm[i_min], cv$nzero[i_min],
                                 cv$lambda.1se, cv$nzero[i_1se]))
  }
  tune <- do.call(rbind, tune)
  best <- which.min(tune$cv_error_min)
  cv <- fits[[best]]
  lam <- if (lambda_rule == "1se") cv$lambda.1se else cv$lambda.min

  list(model = cv$glmnet.fit, cv = cv, alpha = alphas[best], lambda = lam,
       lambda_rule = lambda_rule, tuning = tune, training_mean = mean(y), n_train = length(y),
       preprocess = list(features = features, medians = med, q_lo = q_lo, q_hi = q_hi,
                         max_missing = max_missing, transform = transform_label),
       inner_folds = data.frame(patient_id = patient_id, cancer_type = cancer_type,
                                inner_fold = foldid, stringsAsFactors = FALSE),
       seed = seed, glmnet_version = as.character(packageVersion("glmnet")))
}

# -----------------------------------------------------------------------------
# predict_en: x = samples x probes on the ORIGINAL beta scale, colnames = probes.
# -----------------------------------------------------------------------------
predict_en <- function(b, x, ood_threshold = 0.25) {
  pp <- b$preprocess; f <- pp$features
  xs <- matrix(NA_real_, nrow(x), length(f), dimnames = list(rownames(x), f))
  common <- intersect(f, colnames(x))
  if (length(common)) xs[, common] <- x[, common, drop = FALSE]
  if (identical(pp$transform, "m")) xs <- to_m(xs)
  miss <- rowMeans(is.na(xs))
  if (anyNA(xs)) { idx <- which(is.na(xs), arr.ind = TRUE); xs[idx] <- pp$medians[idx[, 2]] }

  cf <- as.matrix(coef(b$model, s = b$lambda))[-1, 1]
  use <- which(cf != 0); if (!length(use)) use <- seq_along(f)
  outside <- sweep(xs[, use, drop = FALSE], 2, pp$q_lo[use], "<") |
             sweep(xs[, use, drop = FALSE], 2, pp$q_hi[use], ">")
  out_frac <- rowMeans(outside)

  pred <- as.numeric(predict(b$model, newx = xs, s = b$lambda))
  data.frame(sample_id = rownames(x),
             predicted_reference_HRDsum = pmax(0, pred),   # HRDsum cannot be negative
             missing_feature_fraction = miss, outside_range_fraction = out_frac,
             qc_fail = miss > pp$max_missing, ood = out_frac > ood_threshold,
             reportable = !(miss > pp$max_missing) & !(out_frac > ood_threshold),
             stringsAsFactors = FALSE)
}

# -----------------------------------------------------------------------------
auc_rank <- function(label, score) {           # rank-based AUC, no packages
  label <- as.logical(label); n1 <- sum(label); n0 <- sum(!label)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score); (sum(r[label]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

metrics <- function(actual, predicted) {
  n <- length(actual)
  cr <- function(m) if (n >= 3 && sd(actual) > 0 && sd(predicted) > 0)
    suppressWarnings(cor(actual, predicted, method = m)) else NA_real_
  data.frame(n = n, MAE = mean(abs(actual - predicted)),
             RMSE = sqrt(mean((actual - predicted)^2)),
             R2 = if (n >= 2 && var(actual) > 0) 1 - sum((actual - predicted)^2) / sum((actual - mean(actual))^2) else NA_real_,
             pearson = cr("pearson"), spearman = cr("spearman"))
}

# Split-conformal absolute-residual quantile with finite-sample correction
conformal_q <- function(actual, predicted, coverage = 0.95) {
  r <- sort(abs(actual - predicted)); n <- length(r)
  k <- ceiling((n + 1) * coverage)
  if (n == 0 || k > n) Inf else r[k]
}

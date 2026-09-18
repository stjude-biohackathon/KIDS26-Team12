# =============================================================================
# R/rank_model.R - Tissue-RELATIVE target models (Phase C) and a within-tissue
#                  pairwise ranker prototype (Phase D).
# =============================================================================
#
# WHAT THIS FILE IS
# -----------------
# A SEPARATE, self-contained modelling library that sits alongside R/model.R and
# does not modify it. R/model.R predicts ABSOLUTE reference HRDsum. This file
# predicts a TISSUE-RELATIVE quantity: "where does this tumour sit within its own
# tissue's HRD distribution", expressed as a unitless score and, after
# calibration, as a percentile.
#
# WHAT THE OUTPUT IS NOT
# ----------------------
# `relative_score` IS NOT AN ABSOLUTE HRDsum AND CANNOT BE CONVERTED INTO ONE.
# The bundles returned here deliberately store NO per-tissue mean, SD, or
# quantile of y, so back-transformation to the HRDsum scale is impossible by
# construction, not merely discouraged. Anyone who wants absolute HRDsum must use
# R/model.R. The n-of-1 claim this file supports is ordinal ("this patient looks
# high/low relative to a same-tissue reference distribution"), never a scar count.
#
# THE TWO RELATIVE TARGETS
# ------------------------
#   centered      y - mean(y within that sample's TRAINING tissue)
#   normal_score  z = qnorm(q), q = within-tissue midrank percentile
#                 (rank(y, ties="average") - 0.5) / n_tissue
# Both are pointwise functions of the rows handed to relative_target() and of
# nothing else.
#
# THE LEAKAGE BOUNDARY, WHICH IS STRICTER HERE THAN IN R/model.R
# --------------------------------------------------------------
# A relative target is itself a LEARNED quantity (it subtracts/ranks against a
# tissue statistic), so it must obey the same discipline as preprocessing:
#
#   * The final refit uses a target computed on the FULL training set.
#   * Every inner validation fold RECOMPUTES the target from its own inner
#     TRAINING partition only. The transform used to FIT never sees validation y.
#   * Inner-validation rows still need a relative target to be SCORED against.
#     That one is computed from the validation rows' own tissue statistics and is
#     used FOR EVALUATION ONLY - it is a label, never an input to any fit. This
#     mirrors deployment exactly: at deployment the held-out tissue's reference
#     distribution is unknown to the model, and the only way to *evaluate* the
#     ordering is to use the held-out tissue's own labels after the fact.
#     Every place this happens is marked EVAL-ONLY below.
#
# HYPERPARAMETER TUNING IS ALWAYS IN RELATIVE SPACE
# -------------------------------------------------
# The tuning loss is a macro (per-tissue) mean absolute error IN THE RELATIVE
# TARGET SPACE, never absolute HRDsum MAE. Macro 1 - Spearman is recorded as a
# secondary diagnostic. fit_rank_en() asserts this with stopifnot(); see
# RELATIVE_LOSS_COLUMN below.
#
# GBM/LGG
# -------
# GBM and LGG are a LOCKED partition. This library never sees them because the
# caller removes them first; assert_cns_locked() is provided so that removal is
# an executable assertion at the script level rather than a comment.
#
# V3 / feature_rank="within_tissue" - COMPUTATIONAL NOTE
# ------------------------------------------------------
# Ranking probes by POOLED WITHIN-TISSUE variance is done with a blocked
# sufficient-statistic sweep and never materialises a residualised copy of the
# matrix. For each tissue t and each block of columns we accumulate
#   S_t = colSums(x_t), Q_t = colSums(x_t^2)
# and then
#   within_SS = sum_t ( Q_t - S_t^2 / n_t ),  var_within = within_SS / (N - T).
# Peak extra memory is one block (nrow x block_size doubles) plus its square,
# i.e. ~2 x 8 x N x block_size bytes. At N = 7,000 training rows and
# block_size = 20,000 that is ~1.1 GB per temporary, ~2.5 GB peak including the
# non-finite mask - bounded and independent of the 336k total probe count.
# There is NO silent fallback to pooled/global variance: if this path cannot run,
# it errors.
#
# PHASE D (bottom of file) IS A CONTINGENCY PROTOTYPE
# ---------------------------------------------------
# sample_within_tissue_pairs() / fit_pairwise_ranker() / predict_pairwise()
# implement a within-tissue pairwise (RankNet-style logistic) ranker. It has been
# exercised on SYNTHETIC DATA ONLY and is not production code. It is included so
# that, if the pointwise relative target fails at scale, the fallback already
# exists and has a test.
# =============================================================================

source("R/model.R")   # fit_preprocess / apply_preprocess / inner_folds /
                      # lambda_path / check_lambda_boundary are reused verbatim.

# The single permitted name for the hyperparameter-selection loss. Anything else
# (in particular anything mentioning absolute HRDsum) is a bug, and fit_rank_en()
# asserts on this constant.
RELATIVE_LOSS_COLUMN <- "inner_macro_relative_MAE"
CNS_LOCKED <- c("GBM", "LGG")


# -----------------------------------------------------------------------------
# assert_cns_locked(): executable enforcement of the locked CNS partition.
# -----------------------------------------------------------------------------
# Call this immediately after building the train/test masks and BEFORE any target
# transform, fit, calibration, feature selection or summary. It errors - it does
# not warn - if a GBM or LGG sample has reached either side of the split, or if
# the fold's held-out type is itself locked.
assert_cns_locked <- function(cancer_train, cancer_test, held_out_type = NULL) {
  stopifnot(!any(cancer_train %in% CNS_LOCKED),
            !any(cancer_test  %in% CNS_LOCKED))
  if (!is.null(held_out_type)) stopifnot(!any(held_out_type %in% CNS_LOCKED))
  invisible(TRUE)
}


# -----------------------------------------------------------------------------
# relative_target(): the tissue-relative transform. Pure function of its inputs.
# -----------------------------------------------------------------------------
# Input : y       numeric outcome (HRDsum) for the rows handed in
#         cancer  tissue label, same length
#         type    "centered" or "normal_score"
# Output: numeric vector, same length and order as y, carrying attributes
#         "type" and "tissue_stats" (per-tissue n, mean, sd and a constant flag)
#         recording exactly which statistics were used.
#
# This function is deliberately PURE: it uses only the rows it is handed. The
# caller decides whether those rows are an inner-training partition (fitting) or
# an inner-validation partition (EVAL-ONLY scoring). Because the function cannot
# tell the difference, it cannot itself leak - the discipline lives at the call
# sites, all of which are labelled.
#
# normal_score detail: midranks (ties.method="average") mean tied y share a q, so
# tied samples receive an identical z. q = (rank - 0.5)/n is bounded in
# [0.5/n, 1 - 0.5/n] by construction, so qnorm(q) is always finite; the clip below
# is a defensive no-op that would catch a future change to the rank convention.
# A tissue whose y is constant maps to q = 0.5 and z = 0 for ALL its samples -
# handled by an explicit branch so it can never produce NaN.
relative_target <- function(y, cancer, type = c("centered", "normal_score")) {
  type <- match.arg(type)
  stopifnot(length(y) == length(cancer), all(is.finite(y)))
  cancer <- as.character(cancer)

  out <- numeric(length(y))
  groups <- sort(unique(cancer))
  stats <- data.frame(cancer_type = groups, n = NA_integer_, mean = NA_real_,
                      sd = NA_real_, constant = NA, stringsAsFactors = FALSE)

  for (k in seq_along(groups)) {
    i  <- which(cancer == groups[k])
    yt <- y[i]
    nt <- length(yt)
    const <- (nt == 1L) || (max(yt) - min(yt) == 0)
    stats$n[k] <- nt
    stats$mean[k] <- mean(yt)
    stats$sd[k] <- if (nt > 1L) sd(yt) else 0
    stats$constant[k] <- const

    if (type == "centered") {
      # A constant tissue centres to exactly 0, which is the correct answer and
      # needs no special case.
      out[i] <- yt - mean(yt)
    } else {
      if (const) {
        q <- rep(0.5, nt)               # explicit: never NaN, z becomes 0
      } else {
        q <- (rank(yt, ties.method = "average") - 0.5) / nt
        # Defensive clip to the half-rank limits. With midranks this is a no-op.
        q <- pmin(pmax(q, 0.5 / nt), 1 - 0.5 / nt)
      }
      out[i] <- qnorm(q)
    }
  }
  stopifnot(all(is.finite(out)))
  attr(out, "type") <- type
  attr(out, "tissue_stats") <- stats
  out
}


# -----------------------------------------------------------------------------
# tissue_weights(): observation weights giving each tissue equal TOTAL weight.
# -----------------------------------------------------------------------------
# weight_mode "sample": every row weight 1 (the R/model.R behaviour).
# weight_mode "tissue": w_i = 1 / n_{tissue(i)}, rescaled so sum(w) == n. Then
#   tapply(w, cancer, sum) == n / n_tissues for every tissue, i.e. a tissue with
#   800 samples and a tissue with 40 pull equally hard on the fit. Passed to
#   glmnet via weights= (glmnet renormalises to sum to nobs internally, which our
#   rescaling already matches, so the values are interpretable in the log).
tissue_weights <- function(cancer, weight_mode = c("sample", "tissue")) {
  weight_mode <- match.arg(weight_mode)
  cancer <- as.character(cancer)
  n <- length(cancer)
  if (weight_mode == "sample") return(rep(1, n))
  cnt <- table(cancer)
  w <- 1 / as.numeric(cnt[cancer])
  w * (n / sum(w))
}


# -----------------------------------------------------------------------------
# within_tissue_var(): pooled within-tissue variance per probe, blocked.
# -----------------------------------------------------------------------------
# Training rows only - the caller subsets before calling. Non-finite cells are
# filled with the supplied per-probe medians (the same training medians the
# preprocessing uses), so this never invents values and never looks at y.
#
# Identity used (see file header): within_SS = sum_t (Q_t - S_t^2 / n_t), and
# var_within = within_SS / (N - T). Columns are processed in blocks so peak
# memory is O(N * block_size), independent of the total probe count. No
# residualised copy of the matrix is ever created.
within_tissue_var <- function(x, cancer, med = NULL, block_size = 20000L) {
  stopifnot(is.matrix(x), !is.null(colnames(x)), nrow(x) == length(cancer))
  cancer <- as.character(cancer)
  groups <- sort(unique(cancer))
  idx <- lapply(groups, function(g) which(cancer == g))
  nt <- vapply(idx, length, integer(1))
  N <- nrow(x); Tn <- length(groups)
  if (N - Tn < 1L) stop("Need more training rows than tissues for within-tissue variance")

  p <- ncol(x)
  out <- numeric(p); names(out) <- colnames(x)
  starts <- seq.int(1L, p, by = block_size)

  for (st in starts) {
    cols <- st:min(st + block_size - 1L, p)
    xb <- x[, cols, drop = FALSE]
    storage.mode(xb) <- "double"
    nf <- !is.finite(xb)
    if (any(nf)) {
      if (is.null(med)) stop("within_tissue_var(): non-finite values need `med`")
      mb <- med[colnames(xb)]
      if (any(!is.finite(mb))) stop("within_tissue_var(): non-finite median supplied")
      # Fill by column, touching only the affected columns.
      for (j in which(matrixStatsOrBase_colAnys(nf))) xb[nf[, j], j] <- mb[j]
    }
    ss <- numeric(length(cols))
    for (k in seq_along(idx)) {
      i <- idx[[k]]
      if (!length(i)) next
      xt <- xb[i, , drop = FALSE]
      S <- colSums(xt)
      Q <- colSums(xt * xt)
      ss <- ss + (Q - (S * S) / nt[k])
      rm(xt, S, Q)
    }
    out[cols] <- ss / (N - Tn)
    rm(xb, nf, ss)
  }
  # Numerical guard: the identity can return tiny negatives for constant columns.
  out[out < 0 & out > -1e-8] <- 0
  out
}

# Tiny helper so within_tissue_var() can use the compiled reduction when
# matrixStats is available, exactly as fit_preprocess() does, without depending
# on it.
matrixStatsOrBase_colAnys <- function(m) {
  if (requireNamespace("matrixStats", quietly = TRUE)) matrixStats::colAnys(m)
  else apply(m, 2, any)
}


# -----------------------------------------------------------------------------
# fit_preprocess_ranked(): preprocessing with a selectable feature ranking.
# -----------------------------------------------------------------------------
# feature_rank "pooled"        : exactly fit_preprocess() from R/model.R - the
#                                pooled-variance top-max_features filter.
# feature_rank "within_tissue" : same missingness filter, same training medians,
#                                same centre/scale, but probes are ranked by
#                                POOLED WITHIN-TISSUE variance (training rows
#                                only). Implemented by asking fit_preprocess()
#                                for ALL surviving probes, then re-ranking with
#                                within_tissue_var() and subsetting the returned
#                                transform. Nothing is recomputed on held-out
#                                rows and the returned object has the identical
#                                shape, so apply_preprocess() works unchanged.
# Ties on variance break on probe NAME, so selection is deterministic.
fit_preprocess_ranked <- function(x, cancer, max_features = 5000L,
                                  max_missing = 0.05,
                                  feature_rank = c("pooled", "within_tissue"),
                                  block_size = 20000L) {
  feature_rank <- match.arg(feature_rank)
  if (feature_rank == "pooled") return(fit_preprocess(x, max_features, max_missing))

  # Keep every probe that survives missingness + non-constancy, then re-rank.
  pp_all <- fit_preprocess(x, max_features = ncol(x), max_missing = max_missing)
  v <- within_tissue_var(x[, pp_all$features, drop = FALSE], cancer,
                         med = pp_all$median, block_size = block_size)
  ii <- which(is.finite(v) & v > 0)
  if (length(ii) < 2L) stop("Fewer than two probes with non-zero within-tissue variance")
  ii <- ii[order(-v[ii], names(v)[ii])]
  ii <- head(ii, max_features)
  feats <- names(v)[ii]
  list(features = feats, median = pp_all$median[feats],
       center = pp_all$center[feats], scale = pp_all$scale[feats],
       max_missing = max_missing)
}


# -----------------------------------------------------------------------------
# macro_relative_mae(): the tuning loss. Relative space, macro over tissues.
# -----------------------------------------------------------------------------
# truth is a RELATIVE target (centered HRD or normal score); pred is the model's
# relative score. Per-tissue MAE, then an unweighted mean over tissues so a large
# tissue cannot dominate. Returns NA if nothing usable.
macro_relative_mae <- function(truth, pred, cancer) {
  ok <- is.finite(truth) & is.finite(pred)
  if (!any(ok)) return(NA_real_)
  mean(tapply(abs(truth[ok] - pred[ok]), as.character(cancer)[ok], mean))
}

# Secondary diagnostic only: macro (1 - Spearman) in relative space. Tissues with
# fewer than 3 usable rows or no variation are skipped.
macro_one_minus_spearman <- function(truth, pred, cancer) {
  ok <- is.finite(truth) & is.finite(pred)
  if (!any(ok)) return(NA_real_)
  vals <- vapply(split(which(ok), as.character(cancer)[ok]), function(i) {
    if (length(i) < 3L || sd(truth[i]) == 0 || sd(pred[i]) == 0) return(NA_real_)
    1 - suppressWarnings(cor(truth[i], pred[i], method = "spearman"))
  }, numeric(1))
  vals <- vals[is.finite(vals)]
  if (!length(vals)) NA_real_ else mean(vals)
}


# -----------------------------------------------------------------------------
# calibration_knots(): percentile curve from cross-fitted SOURCE-tissue scores.
# -----------------------------------------------------------------------------
# Input is the inner-CV OUT-OF-FOLD relative scores for the training (source)
# tissues at the selected hyperparameters - i.e. every score was produced by a
# model that had not seen that sample. Held-out-tissue data never enters.
#
# Stored as at most `max_knots` (score, percentile) pairs so the bundle stays
# small; percentiles are midranks of the OOF distribution. score_to_percentile()
# interpolates this curve, so a single sample can be placed without any cohort.
calibration_knots <- function(oof_scores, max_knots = 501L) {
  s <- sort(oof_scores[is.finite(oof_scores)])
  n <- length(s)
  if (n < 5L) return(list(score = numeric(0), percentile = numeric(0), n = n))
  p <- (seq_len(n) - 0.5) / n
  if (n > max_knots) {
    sel <- unique(round(seq(1, n, length.out = max_knots)))
    s <- s[sel]; p <- p[sel]
  }
  keep <- !duplicated(s)          # approx() needs strictly increasing x
  list(score = s[keep], percentile = p[keep], n = n)
}


# -----------------------------------------------------------------------------
# score_to_percentile(): raw relative score -> percentile in (0,1).
# -----------------------------------------------------------------------------
# Uses ONLY the knots frozen into the bundle at fit time (source-tissue
# cross-fitted predictions). Single-sample safe: no cohort statistic is computed
# here, and a length-1 input is legal. rule=2 clamps beyond the observed range to
# the extreme knots rather than returning NA, which is the honest behaviour for a
# score outside everything seen in development (it is reported as "at or beyond
# the most extreme source-tissue score", not extrapolated).
score_to_percentile <- function(bundle, raw_scores) {
  k <- bundle$calibration
  if (is.null(k) || length(k$score) < 2L) return(rep(NA_real_, length(raw_scores)))
  as.numeric(approx(x = k$score, y = k$percentile, xout = as.numeric(raw_scores),
                    rule = 2, ties = "ordered")$y)
}


# -----------------------------------------------------------------------------
# fit_rank_en(): the Phase C model. Variants V1/V2/V3 are argument settings.
# -----------------------------------------------------------------------------
#   V1 = feature_rank "pooled",        weight_mode "sample"
#   V2 = feature_rank "pooled",        weight_mode "tissue"
#   V3 = feature_rank "within_tissue", weight_mode "tissue"
#
# Everything passed in is TRAINING data for the current outer fold: the held-out
# tissue is not visible to this function, and neither are GBM/LGG (the caller
# asserts that with assert_cns_locked()).
#
# Structure mirrors fit_en():
#   Phase 1  inner CV over (alpha, lambda), loss = macro MAE in RELATIVE space.
#            The relative target is RECOMPUTED inside each inner fold from the
#            inner-TRAINING rows only; the inner-validation relative target is
#            computed from the validation rows' own tissue statistics and used
#            EVAL-ONLY (marked at the call site).
#   Phase 2  refit preprocessing + glmnet on all supplied rows with the winning
#            hyperparameters and the full-training relative target.
#   Phase 3  freeze the calibration curve from the Phase 1 out-of-fold scores.
fit_rank_en <- function(x, y, patient, cancer,
                        target_type = c("normal_score", "centered"),
                        weight_mode = c("sample", "tissue"),
                        feature_rank = c("pooled", "within_tissue"),
                        max_features = 5000L, seed = 260910L,
                        block_size = 20000L) {
  target_type  <- match.arg(target_type)
  weight_mode  <- match.arg(weight_mode)
  feature_rank <- match.arg(feature_rank)
  if (!requireNamespace("glmnet", quietly = TRUE)) stop("Install glmnet via scripts/setup.R")
  stopifnot(is.matrix(x), length(y) == nrow(x), all(is.finite(y)),
            length(patient) == length(y), length(cancer) == length(y))
  # The locked partition can never reach a fit, even if a caller forgets.
  stopifnot(!any(as.character(cancer) %in% CNS_LOCKED))
  cancer <- as.character(cancer)

  folds <- inner_folds(patient, cancer, seed)
  fold_ids <- sort(unique(folds))

  # --- Relative target for the FINAL REFIT: full training set ----------------
  y_rel_full <- relative_target(y, cancer, target_type)
  if (sd(y_rel_full) == 0) stop("Constant relative target on the training set")

  # Preprocessing + lambda paths on all training rows (same hoisting trick as
  # fit_en: this object is reused for Phase 2 and is never used to score an inner
  # validation fold).
  pp <- fit_preprocess_ranked(x, cancer, max_features, feature_rank = feature_rank,
                              block_size = block_size)
  z  <- apply_preprocess(x, pp, FALSE)

  alphas <- c(0.1, 0.5, 1)
  lam_by_alpha <- lapply(alphas, function(a) lambda_path(z, y_rel_full, a))
  names(lam_by_alpha) <- as.character(alphas)
  grid <- do.call(rbind, lapply(alphas, function(a)
    data.frame(alpha = a, lambda = lam_by_alpha[[as.character(a)]])))

  losses  <- matrix(NA_real_, nrow(grid), length(fold_ids))
  rho_los <- matrix(NA_real_, nrow(grid), length(fold_ids))
  oof     <- matrix(NA_real_, length(y), nrow(grid))   # cross-fitted scores

  # ---- Phase 1: inner cross-validation --------------------------------------
  for (f in fold_ids) {
    itr <- folds != f; iva <- !itr

    # Preprocessing relearned from this inner fold's TRAINING rows only.
    pp_in <- fit_preprocess_ranked(x[itr, , drop = FALSE], cancer[itr], max_features,
                                   feature_rank = feature_rank, block_size = block_size)
    ztr <- apply_preprocess(x[itr, , drop = FALSE], pp_in, FALSE)
    zva <- apply_preprocess(x[iva, , drop = FALSE], pp_in, FALSE)

    # FIT TARGET: recomputed from the inner-TRAINING rows only. The validation
    # rows' labels are not involved in any way.
    y_tr <- relative_target(y[itr], cancer[itr], target_type)

    # EVAL-ONLY TARGET: the validation rows' relative target, computed from the
    # VALIDATION rows' own tissue statistics. This is the label we score against
    # and it is never fed into glmnet, the preprocessing, the weights, or the
    # lambda path. It exists because a relative score can only be graded against
    # a relative truth, exactly as at deployment.
    y_va_eval <- relative_target(y[iva], cancer[iva], target_type)

    w_tr <- tissue_weights(cancer[itr], weight_mode)

    if (sd(y_tr) == 0) next
    for (a in alphas) {
      lam_a <- sort(lam_by_alpha[[as.character(a)]], decreasing = TRUE)
      mod <- glmnet::glmnet(ztr, y_tr, alpha = a, lambda = lam_a,
                            weights = w_tr, standardize = FALSE)
      ids <- which(grid$alpha == a)
      for (g in ids) {
        pred <- as.numeric(predict(mod, zva, s = grid$lambda[g]))
        oof[iva, g]  <- pred
        losses[g, match(f, fold_ids)]  <- macro_relative_mae(y_va_eval, pred, cancer[iva])
        rho_los[g, match(f, fold_ids)] <- macro_one_minus_spearman(y_va_eval, pred, cancer[iva])
      }
    }
  }

  # ---- Select hyperparameters. RELATIVE-SPACE LOSS ONLY. --------------------
  grid[[RELATIVE_LOSS_COLUMN]]          <- rowMeans(losses, na.rm = TRUE)
  grid$inner_macro_one_minus_spearman   <- rowMeans(rho_los, na.rm = TRUE)
  # EXECUTABLE GUARD: the objective must be the relative-space column and must
  # not be any absolute-HRDsum quantity. If someone renames or repoints the loss,
  # this fails loudly instead of silently tuning on the wrong scale.
  loss_col <- RELATIVE_LOSS_COLUMN
  stopifnot(identical(loss_col, "inner_macro_relative_MAE"),
            loss_col %in% names(grid),
            !grepl("HRDsum", loss_col, fixed = TRUE),
            !grepl("absolute", loss_col, fixed = TRUE),
            target_type %in% c("centered", "normal_score"))
  if (all(!is.finite(grid[[loss_col]]))) stop("Inner CV produced no usable relative-space loss")

  best <- which.min(grid[[loss_col]])
  best_alpha <- grid$alpha[best]; best_lambda <- grid$lambda[best]
  best_path <- lam_by_alpha[[as.character(best_alpha)]]
  boundary <- check_lambda_boundary(best_lambda, best_path, best_alpha,
                                    label = sprintf("fit_rank_en (n=%d)", length(y)))

  # ---- Phase 2: refit on all training rows ----------------------------------
  w_full <- tissue_weights(cancer, weight_mode)
  stopifnot(abs(sum(w_full) - length(y)) < 1e-8)
  mod <- glmnet::glmnet(z, y_rel_full, alpha = best_alpha,
                        lambda = sort(best_path, decreasing = TRUE),
                        weights = w_full, standardize = FALSE)

  dist <- sqrt(rowMeans(z^2)); ood_cut <- as.numeric(quantile(dist, 0.99, names = FALSE))

  # ---- Phase 3: freeze the calibration curve --------------------------------
  # SOURCE TISSUES, CROSS-FITTED ONLY: column `best` of `oof` holds, for each
  # training sample, the score from a model fitted without that sample's inner
  # fold. No held-out-tissue information and no in-sample scores.
  cal <- calibration_knots(oof[, best])

  bundle <- list(
    preprocess   = pp,
    model        = mod,
    alpha        = best_alpha,
    lambda       = best_lambda,
    target_type  = target_type,
    weight_mode  = weight_mode,
    feature_rank = feature_rank,
    tuning       = grid,
    tuning_objective = loss_col,
    lambda_paths = lam_by_alpha,
    lambda_at_boundary = if (is.null(boundary)) NA else (boundary$at_top || boundary$at_bottom),
    lambda_boundary_side = if (is.null(boundary)) NA_character_
                           else if (boundary$at_top) "top" else if (boundary$at_bottom) "bottom" else "interior",
    inner_folds  = data.frame(patient_id = patient, cancer_type = cancer, fold = folds,
                              stringsAsFactors = FALSE),
    oof_score    = oof[, best],
    calibration  = cal,
    ood_cut      = ood_cut,
    seed         = seed,
    training_n   = length(y),
    training_tissues = sort(unique(cancer)),
    # The output is a RELATIVE score. No tissue means, SDs or quantiles of y are
    # stored anywhere in this bundle, so absolute HRDsum cannot be reconstructed.
    score_type   = "relative_score",
    target       = "tissue_relative_HRD_rank"
  )
  # Executable statement of the no-back-transformation property.
  stopifnot(identical(bundle$score_type, "relative_score"),
            !any(c("training_mean", "tissue_means", "target_transform",
                   "y_center", "y_scale") %in% names(bundle)))
  bundle
}


# -----------------------------------------------------------------------------
# predict_rank(): score new samples with a frozen rank bundle. Learns nothing.
# -----------------------------------------------------------------------------
# Returns one row per input sample. `relative_score` is the model output;
# `percentile` places it on the frozen source-tissue calibration curve. Works on
# a 1-row matrix and requires no cohort statistic of any kind - that property is
# asserted in tests/test_c1_rank_model.R (test 6).
predict_rank <- function(bundle, x) {
  if (!is.matrix(x)) stop("predict_rank() expects a matrix with probe colnames")
  z <- apply_preprocess(x, bundle$preprocess, FALSE)
  raw <- as.numeric(predict(bundle$model, z, s = bundle$lambda))
  missing <- attr(z, "missing_fraction")
  score <- sqrt(rowMeans(z^2))
  fail <- missing > bundle$preprocess$max_missing
  ood  <- score > bundle$ood_cut
  data.frame(sample_id = rownames(x),
             relative_score = raw,
             percentile = score_to_percentile(bundle, raw),
             missing_fraction = missing,
             ood_score = score, ood = ood, qc_fail = fail,
             reportable = !(fail | ood),
             stringsAsFactors = FALSE)
}


# -----------------------------------------------------------------------------
# rank_null_panel(): is the within-tissue ordering better than chance?
# -----------------------------------------------------------------------------
# Holds the predicted relative scores fixed and permutes the TRUE relative target
# within tissue. Reports the observed Spearman against the permutation null.
# Like permutation_test_within_tissue() in R/model.R this cannot detect a pipeline
# that manufactures signal from noise (that needs a refit per permutation); it
# asks whether THESE scores carry within-tissue information.
rank_null_panel <- function(truth, score, cancer, n_perm = 500L, seed = 260910L) {
  ok <- is.finite(truth) & is.finite(score)
  truth <- truth[ok]; score <- score[ok]; cancer <- as.character(cancer)[ok]
  n <- length(truth)
  empty <- data.frame(n = n, observed_spearman = NA_real_, null_spearman_median = NA_real_,
                      null_spearman_q025 = NA_real_, null_spearman_q975 = NA_real_,
                      perm_p_value = NA_real_, n_perm = 0L)
  if (n < 5L || sd(truth) == 0) return(empty)
  # A CONSTANT score claims no ordering at all. Spearman is undefined there, and
  # the honest summary is "zero observed association, p = 1" rather than NA -
  # otherwise a fully shrunk model (the correct answer on pure noise) looks like
  # a broken one. This branch is exercised by test 10.
  if (sd(score) == 0) {
    return(data.frame(n = n, observed_spearman = 0, null_spearman_median = 0,
                      null_spearman_q025 = 0, null_spearman_q975 = 0,
                      perm_p_value = 1, n_perm = 0L))
  }
  obs <- suppressWarnings(cor(truth, score, method = "spearman"))
  set.seed(seed)
  nulls <- vapply(seq_len(n_perm), function(k)
    suppressWarnings(cor(permute_within_tissue(truth, cancer), score, method = "spearman")),
    numeric(1))
  nulls <- nulls[is.finite(nulls)]
  data.frame(n = n, observed_spearman = obs,
             null_spearman_median = median(nulls),
             null_spearman_q025 = as.numeric(quantile(nulls, 0.025, names = FALSE)),
             null_spearman_q975 = as.numeric(quantile(nulls, 0.975, names = FALSE)),
             perm_p_value = (1 + sum(abs(nulls) >= abs(obs))) / (length(nulls) + 1),
             n_perm = length(nulls))
}


# -----------------------------------------------------------------------------
# rank_metrics(): evaluation panel for one held-out tissue.
# -----------------------------------------------------------------------------
# truth is the held-out rows' own relative target (EVAL-ONLY, computed after the
# fact from the held-out tissue's labels - unavailable at deployment).
rank_metrics <- function(truth, score, cancer) {
  ok <- is.finite(truth) & is.finite(score)
  truth <- truth[ok]; score <- score[ok]; cancer <- as.character(cancer)[ok]
  n <- length(truth)
  variable <- n > 2 && sd(truth) > 0 && sd(score) > 0
  data.frame(n = n,
             within_Spearman = if (variable) suppressWarnings(cor(truth, score, method = "spearman")) else NA_real_,
             within_Pearson  = if (variable) suppressWarnings(cor(truth, score)) else NA_real_,
             macro_relative_MAE = macro_relative_mae(truth, score, cancer),
             score_sd = if (n > 1) sd(score) else NA_real_)
}


# =============================================================================
# PHASE D - WITHIN-TISSUE PAIRWISE RANKER (CONTINGENCY PROTOTYPE)
# =============================================================================
# STATUS: PROTOTYPE. Exercised on synthetic data only (tests 11). Not wired into
# any production script and not validated at scale. It exists as a fallback in
# case the pointwise relative target in Phase C proves too weak.
#
# IDEA: the clinical question is ordinal and within-tissue ("of two tumours of
# the same type, which has more scar?"). Train directly on that: features are
# x_i - x_j for two samples of the SAME tissue, outcome is I(y_i > y_j),
# penalised logistic regression. Because every pair is within-tissue, any
# tissue-level offset cancels in the difference exactly - the tissue confound is
# removed by construction rather than by a learned centring.
#
# A single new sample is scored as beta'x. The intercept cancels in differences,
# so the model is fitted with intercept=FALSE and its absence is asserted; the
# resulting score is only meaningful up to a monotone transform, which is why it
# is reported as a percentile against source-tissue cross-fitted scores.
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# sample_within_tissue_pairs(): bounded, balanced, strictly within-tissue pairs.
# -----------------------------------------------------------------------------
# For each tissue: enumerate (if small) or randomly sample (if large) at most
# max_pairs_per_tissue pairs of DISTINCT samples with DISTINCT y. Tissues with
# fewer than 2 distinct y values are skipped - they can supply no ordered pair.
#
# BALANCE: each sampled pair is first oriented so that y[i] > y[j] (outcome 1),
# then a balanced random half are flipped. That is a random orientation of every
# pair, but with the number of flips fixed at floor(n/2) so the outcome mean is
# 0.5 by construction rather than 0.5 in expectation - important for small
# tissues where binomial noise would otherwise unbalance the classes.
sample_within_tissue_pairs <- function(y, cancer, max_pairs_per_tissue = 2000L,
                                       seed = 260910L) {
  stopifnot(length(y) == length(cancer), all(is.finite(y)))
  cancer <- as.character(cancer)
  set.seed(seed)
  out <- list()

  for (g in sort(unique(cancer))) {
    idx <- which(cancer == g)
    if (length(idx) < 2L || length(unique(y[idx])) < 2L) next   # no ordered pair
    n <- length(idx)
    total <- n * (n - 1L) / 2L

    if (total <= max_pairs_per_tissue) {
      cmb <- utils::combn(n, 2L)
      a <- idx[cmb[1, ]]; b <- idx[cmb[2, ]]
    } else {
      # Oversample random index pairs, drop self-pairs and duplicates, then trim.
      m <- min(as.integer(max_pairs_per_tissue * 3L), 200000L)
      ia <- sample.int(n, m, replace = TRUE); ib <- sample.int(n, m, replace = TRUE)
      keep <- ia != ib
      ia <- ia[keep]; ib <- ib[keep]
      lo <- pmin(ia, ib); hi <- pmax(ia, ib)
      key <- paste(lo, hi, sep = "_")
      dup <- duplicated(key)
      lo <- lo[!dup]; hi <- hi[!dup]
      take <- seq_len(min(length(lo), max_pairs_per_tissue))
      a <- idx[lo[take]]; b <- idx[hi[take]]
    }

    # Drop tied pairs: I(y_i > y_j) is not defined for them.
    keep <- y[a] != y[b]
    a <- a[keep]; b <- b[keep]
    if (!length(a)) next
    if (length(a) > max_pairs_per_tissue) {
      s <- sample.int(length(a), max_pairs_per_tissue)
      a <- a[s]; b <- b[s]
    }

    # Orient to outcome 1, then flip a balanced random half.
    hi_first <- y[a] > y[b]
    i <- ifelse(hi_first, a, b); j <- ifelse(hi_first, b, a)
    np <- length(i)
    flip <- sample(rep(c(TRUE, FALSE), length.out = np))
    ii <- ifelse(flip, j, i); jj <- ifelse(flip, i, j)

    out[[g]] <- data.frame(i = ii, j = jj,
                           outcome = as.integer(y[ii] > y[jj]),
                           cancer_type = g, stringsAsFactors = FALSE)
  }
  if (!length(out)) stop("No usable within-tissue pairs")
  res <- do.call(rbind, out); rownames(res) <- NULL
  # Executable invariants: within-tissue only, bounded, no self-pairs.
  stopifnot(all(cancer[res$i] == cancer[res$j]),
            all(res$i != res$j),
            all(table(res$cancer_type) <= max_pairs_per_tissue))
  res
}


# -----------------------------------------------------------------------------
# fit_pairwise_ranker(): penalised logistic regression on within-tissue
#                        difference vectors. PROTOTYPE.
# -----------------------------------------------------------------------------
# Observation weights give each tissue equal TOTAL weight regardless of how many
# pairs it contributed, so a large tissue cannot dominate the ranker.
# Hyperparameters are chosen by leave-one-tissue-out inner CV on PAIR
# classification error, macro-averaged over held-out tissues - a relative-space
# objective by construction (every pair is within-tissue).
# The calibration curve is built from cross-fitted SAMPLE-level scores on source
# tissues, reusing calibration_knots()/score_to_percentile().
fit_pairwise_ranker <- function(x, y, cancer, patient = NULL,
                                max_features = 5000L, max_pairs_per_tissue = 2000L,
                                seed = 260910L) {
  if (!requireNamespace("glmnet", quietly = TRUE)) stop("Install glmnet via scripts/setup.R")
  stopifnot(is.matrix(x), length(y) == nrow(x), all(is.finite(y)))
  stopifnot(!any(as.character(cancer) %in% CNS_LOCKED))
  cancer <- as.character(cancer)

  pp <- fit_preprocess(x, max_features)
  z  <- apply_preprocess(x, pp, FALSE)
  pairs <- sample_within_tissue_pairs(y, cancer, max_pairs_per_tissue, seed)

  pair_weights <- function(ct) {
    cnt <- table(ct); w <- 1 / as.numeric(cnt[ct]); w * (length(ct) / sum(w))
  }

  alphas <- c(0.1, 0.5, 1)
  d_all <- z[pairs$i, , drop = FALSE] - z[pairs$j, , drop = FALSE]
  lam_by_alpha <- lapply(alphas, function(a) lambda_path(d_all, pairs$outcome, a))
  names(lam_by_alpha) <- as.character(alphas)
  grid <- do.call(rbind, lapply(alphas, function(a)
    data.frame(alpha = a, lambda = lam_by_alpha[[as.character(a)]])))

  tissues <- sort(unique(cancer))
  losses <- matrix(NA_real_, nrow(grid), length(tissues))
  oof <- matrix(NA_real_, nrow(x), nrow(grid))

  for (k in seq_along(tissues)) {
    tt <- tissues[k]
    tr_rows <- cancer != tt
    if (!any(tr_rows)) next
    p_tr <- pairs[pairs$cancer_type != tt, , drop = FALSE]
    p_va <- pairs[pairs$cancer_type == tt, , drop = FALSE]
    if (!nrow(p_tr) || !nrow(p_va)) next

    # Preprocessing relearned on inner-training rows only.
    pp_in <- fit_preprocess(x[tr_rows, , drop = FALSE], max_features)
    ztr <- apply_preprocess(x[tr_rows, , drop = FALSE], pp_in, FALSE)
    zall <- apply_preprocess(x, pp_in, FALSE)
    # Map original row indices onto positions within the inner-training block.
    row_map <- which(tr_rows)
    pos <- match(p_tr$i, row_map); neg <- match(p_tr$j, row_map)
    stopifnot(!anyNA(pos), !anyNA(neg))
    d_tr <- ztr[pos, , drop = FALSE] - ztr[neg, , drop = FALSE]
    w_tr <- pair_weights(p_tr$cancer_type)
    d_va <- zall[p_va$i, , drop = FALSE] - zall[p_va$j, , drop = FALSE]

    for (a in alphas) {
      lam_a <- sort(lam_by_alpha[[as.character(a)]], decreasing = TRUE)
      mod <- try(glmnet::glmnet(d_tr, p_tr$outcome, family = "binomial", alpha = a,
                                lambda = lam_a, weights = w_tr,
                                standardize = FALSE, intercept = FALSE), silent = TRUE)
      if (inherits(mod, "try-error")) next
      ids <- which(grid$alpha == a)
      for (g in ids) {
        pv <- as.numeric(predict(mod, d_va, s = grid$lambda[g], type = "link"))
        losses[g, k] <- mean((pv > 0) != (p_va$outcome == 1L))   # pair error rate
        # Sample-level cross-fitted score for the held-out tissue's own rows.
        oof[cancer == tt, g] <- as.numeric(predict(mod, zall[cancer == tt, , drop = FALSE],
                                                   s = grid$lambda[g], type = "link"))
      }
    }
  }

  grid$inner_macro_pair_error <- rowMeans(losses, na.rm = TRUE)
  if (all(!is.finite(grid$inner_macro_pair_error))) stop("Pairwise inner CV produced no usable loss")
  best <- which.min(grid$inner_macro_pair_error)
  best_alpha <- grid$alpha[best]; best_lambda <- grid$lambda[best]

  w_all <- pair_weights(pairs$cancer_type)
  mod <- glmnet::glmnet(d_all, pairs$outcome, family = "binomial", alpha = best_alpha,
                        lambda = sort(lam_by_alpha[[as.character(best_alpha)]], decreasing = TRUE),
                        weights = w_all, standardize = FALSE, intercept = FALSE)
  # The intercept must be exactly zero: it cancels in the differences and would
  # otherwise add a constant to every single-sample score.
  a0 <- as.numeric(coef(mod, s = best_lambda)[1])
  stopifnot(is.finite(a0), abs(a0) < 1e-10)

  bundle <- list(preprocess = pp, model = mod, alpha = best_alpha, lambda = best_lambda,
                 tuning = grid, tuning_objective = "inner_macro_pair_error",
                 n_pairs = nrow(pairs),
                 pairs_per_tissue = as.data.frame(table(pairs$cancer_type),
                                                  stringsAsFactors = FALSE),
                 calibration = calibration_knots(oof[, best]),
                 oof_score = oof[, best], seed = seed, training_n = nrow(x),
                 max_pairs_per_tissue = max_pairs_per_tissue,
                 score_type = "relative_score", model_kind = "pairwise_prototype")
  stopifnot(!any(c("training_mean", "tissue_means") %in% names(bundle)))
  bundle
}


# -----------------------------------------------------------------------------
# predict_pairwise(): raw beta'x for one or many samples. PROTOTYPE.
# -----------------------------------------------------------------------------
# No intercept term (see fit_pairwise_ranker), no cohort statistic, single-sample
# safe. Percentile comes from the frozen source-tissue cross-fitted curve.
predict_pairwise <- function(bundle, x) {
  if (!is.matrix(x)) stop("predict_pairwise() expects a matrix with probe colnames")
  z <- apply_preprocess(x, bundle$preprocess, FALSE)
  raw <- as.numeric(predict(bundle$model, z, s = bundle$lambda, type = "link"))
  data.frame(sample_id = rownames(x), relative_score = raw,
             percentile = score_to_percentile(bundle, raw),
             missing_fraction = attr(z, "missing_fraction"),
             stringsAsFactors = FALSE)
}

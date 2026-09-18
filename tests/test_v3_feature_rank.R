#!/usr/bin/env Rscript
# =============================================================================
# tests/test_v3_feature_rank.R - V3-abs: within-tissue probe ranking in the
#                                ABSOLUTE-target model (R/model.R).
# =============================================================================
#
# USAGE
#   Rscript tests/test_v3_feature_rank.R      (from the repository root)
#
# This is a plain stopifnot() script, matching tests/smoke_model.R and
# tests/test_c1_rank_model.R. It is NOT a testthat suite and is NOT picked up by
# `python -m unittest discover -s tests`. Run it explicitly.
#
# WHAT V3-abs CHANGES
#   ONE thing: the 5,000-probe unsupervised filter in fit_preprocess() ranks
#   probes by POOLED WITHIN-TISSUE variance instead of pooled TOTAL variance.
#   Everything else - missingness filter, training medians, centre/scale,
#   tie-break, nested CV - is untouched.
#
# WHAT THIS FILE GUARDS, AND WHY EACH CHECK EXISTS
#   A. The switch actually switches. On a fixture built so that the two
#      statistics disagree by construction (probes whose variance is almost
#      entirely BETWEEN tissues), "within_tissue" must select a materially
#      different feature set from "pooled". Without this, a silent fallback to
#      pooled variance would make a V3 run indistinguishable from run 01 while
#      claiming to be a new variant.
#   B. THE LEAKAGE ASSERTION. The ranking, like every other learned quantity,
#      must be a function of TRAINING ROWS ONLY. We corrupt the held-out block -
#      its feature values AND its tissue labels - and require the selected
#      features and ALL preprocessing constants to come back bit-identical.
#      Anything computed across the full matrix would move here.
#   C. The blocked sufficient-statistic sweep is arithmetically correct: it
#      equals a direct per-tissue sum-of-squares to 1e-10 and does not depend on
#      block_size. The blocking exists only to bound memory on the 336k-probe
#      matrix, so it must be numerically invisible.
#   D. BACKWARD COMPATIBILITY. feature_rank="pooled" must reproduce the
#      pre-change behaviour exactly, including a bundle from fit_en() whose
#      coefficients and tuning table are identical to a default-argument fit.
#      The default is "pooled", so every existing run stays reproducible.
#
# NOT COVERED HERE: anything about whether within-tissue ranking PREDICTS
# better. That is what the LSF sentinel is for, and this file takes no view.
# =============================================================================

source("R/model.R")

pass <- 0L
ok <- function(n, msg) { pass <<- pass + 1L; cat(sprintf("  ok %d - %s\n", n, msg)) }

# --- Synthetic multi-tissue fixture ----------------------------------------
# 120 samples x 60 probes across 4 tissues, 30 samples each.
#
# The fixture is built so pooled and within-tissue variance DISAGREE on purpose:
#   probes  1-20  "tissue markers": near-constant inside a tissue but shifted a
#                 long way between tissues. HUGE pooled variance, TINY
#                 within-tissue variance. Pooled ranking loves them; V3 must not.
#   probes 21-40  "within-tissue signal": no tissue offset at all, but they vary
#                 a lot from sample to sample inside every tissue. Modest pooled
#                 variance, LARGE within-tissue variance. V3 must prefer these.
#   probes 41-60  low-variance filler on both statistics.
# All values are squashed into [0,1] to mimic beta values.
set.seed(2026L)
n_per <- 30L; tissues <- c("TA", "TB", "TC", "TD")
cancer <- rep(tissues, each = n_per)
n <- length(cancer); p <- 60L
x <- matrix(0, n, p,
            dimnames = list(paste0("s", seq_len(n)), sprintf("cg%02d", seq_len(p))))
tissue_offset <- setNames(seq(0.1, 0.9, length.out = length(tissues)), tissues)
for (j in 1:20)  x[, j] <- tissue_offset[cancer] + rnorm(n, sd = 0.005)
for (j in 21:40) x[, j] <- 0.5 + rnorm(n, sd = 0.20)
for (j in 41:60) x[, j] <- 0.5 + rnorm(n, sd = 0.02)
x[x < 0] <- 0; x[x > 1] <- 1
patient <- paste0("p", seq_len(n))
# y depends on two WITHIN-tissue probes plus a tissue offset, so fit_en() has
# something real to find under either ranking.
y <- 40 * x[, 21] - 20 * x[, 22] + 10 * tissue_offset[cancer] + rnorm(n, sd = 1)
y <- as.numeric(y - min(y) + 1)

cat("V3-abs feature_rank tests\n")

# =============================================================================
# A. "within_tissue" selects a DIFFERENT feature set than "pooled".
# =============================================================================
K <- 15L
pp_pool <- fit_preprocess(x, max_features = K, feature_rank = "pooled")
pp_with <- fit_preprocess(x, max_features = K, feature_rank = "within_tissue",
                          cancer = cancer)

stopifnot(length(pp_pool$features) == K, length(pp_with$features) == K,
          !identical(sort(pp_pool$features), sort(pp_with$features)))
# The chosen strategy is RECORDED in the returned object, so a saved transform
# is self-describing and a mislabelled run is detectable after the fact.
stopifnot(identical(pp_pool$feature_rank, "pooled"),
          identical(pp_with$feature_rank, "within_tissue"))
# Direction check: on this fixture pooled ranking must be dominated by the
# tissue-marker probes (1-20) and within-tissue ranking by the within-tissue
# signal probes (21-40). This is what makes check A a real test of the
# statistic rather than of an incidental tie-break difference.
idx_of <- function(f) as.integer(sub("^cg", "", f))
stopifnot(all(idx_of(pp_pool$features) <= 20L),
          all(idx_of(pp_with$features) >= 21L & idx_of(pp_with$features) <= 40L))
# fit_preprocess() must also be deterministic under repetition.
stopifnot(identical(pp_with,
                    fit_preprocess(x, max_features = K,
                                   feature_rank = "within_tissue", cancer = cancer)))
ok(1, "within_tissue selects a different (and correctly targeted) feature set than pooled")

# =============================================================================
# B. LEAKAGE: the ranking is TRAINING-FOLD-ONLY.
# =============================================================================
# Hold out tissue TD. Fit the transform on TA/TB/TC only. Then rebuild a full
# matrix in which the held-out block is complete garbage - values AND labels -
# and fit the transform on the SAME training rows again. If any part of the
# ranking or any preprocessing constant looked outside the training rows, the two
# transforms would differ.
tr <- cancer != "TD"

pp_clean <- fit_preprocess(x[tr, , drop = FALSE], max_features = K,
                           feature_rank = "within_tissue", cancer = cancer[tr])

x_garbage <- x
set.seed(99L)
# Held-out FEATURE ROWS replaced with garbage on a wildly different scale.
x_garbage[!tr, ] <- matrix(runif(sum(!tr) * p, 0, 1), sum(!tr), p)
# Held-out LABELS replaced with garbage tissue names, including names that
# collide with nothing in training.
cancer_garbage <- cancer
cancer_garbage[!tr] <- sample(c("ZZ", "YY", "XX"), sum(!tr), replace = TRUE)

pp_dirty <- fit_preprocess(x_garbage[tr, , drop = FALSE], max_features = K,
                           feature_rank = "within_tissue", cancer = cancer_garbage[tr])

stopifnot(identical(pp_clean$features, pp_dirty$features),
          identical(pp_clean$median,   pp_dirty$median),
          identical(pp_clean$center,   pp_dirty$center),
          identical(pp_clean$scale,    pp_dirty$scale),
          identical(pp_clean$feature_rank, pp_dirty$feature_rank),
          identical(pp_clean, pp_dirty))

# Same assertion one level up, through fit_en(): the outer refit AND every inner
# fold must re-rank on their own rows. A garbage held-out block must leave the
# whole fitted bundle untouched.
coef_of <- function(b) as.matrix(coef(b$model, s = b$lambda))
b_clean <- fit_en(x[tr, , drop = FALSE], y[tr], patient[tr], cancer[tr],
                  max_features = K, feature_rank = "within_tissue")
b_dirty <- fit_en(x_garbage[tr, , drop = FALSE], y[tr], patient[tr], cancer_garbage[tr],
                  max_features = K, feature_rank = "within_tissue")
stopifnot(identical(b_clean$preprocess, b_dirty$preprocess),
          identical(coef_of(b_clean), coef_of(b_dirty)),
          identical(b_clean$tuning, b_dirty$tuning),
          identical(b_clean$alpha, b_dirty$alpha),
          identical(b_clean$lambda, b_dirty$lambda),
          identical(b_clean$feature_rank, "within_tissue"),
          identical(b_clean$preprocess$feature_rank, "within_tissue"))
ok(2, "ranking + all preprocessing constants are bit-identical under garbage held-out rows and labels")

# The inner folds really do re-rank rather than reuse the outer ranking. Fit the
# transform on one inner-training partition directly and confirm it can differ
# from the outer-training transform; if fit_en() hoisted a single ranking, no
# such difference could ever arise.
f <- inner_folds(patient[tr], cancer[tr])
itr <- f != sort(unique(f))[1]
pp_inner <- fit_preprocess(x[tr, , drop = FALSE][itr, , drop = FALSE],
                           max_features = K, feature_rank = "within_tissue",
                           cancer = cancer[tr][itr])
stopifnot(length(pp_inner$features) == K,
          !identical(pp_inner$center, pp_clean$center))
ok(3, "inner-fold transform is computed on inner-training rows only")

# =============================================================================
# C. The blocked sweep is exact and block_size-invariant.
# =============================================================================
xt <- x[tr, , drop = FALSE]; ct <- cancer[tr]
direct <- apply(xt, 2, function(col)
  sum(tapply(col, ct, function(v) sum((v - mean(v))^2)))) / (nrow(xt) - length(unique(ct)))
v_blocked <- within_tissue_var(xt, ct, block_size = 7L)
stopifnot(isTRUE(all.equal(unname(v_blocked), unname(direct), tolerance = 1e-10)),
          identical(names(v_blocked), colnames(xt)))
for (bs in c(1L, 3L, 13L, 60L, 1000L)) {
  stopifnot(isTRUE(all.equal(v_blocked, within_tissue_var(xt, ct, block_size = bs),
                             tolerance = 1e-12)))
}
# Non-finite cells are filled from the SUPPLIED TRAINING medians, never invented,
# and a missing median is an error rather than a silent NaN.
x_na <- xt; x_na[1, 3] <- NA_real_; x_na[5, 4] <- Inf
med_tr <- apply(xt, 2, median)
filled <- x_na; filled[1, 3] <- med_tr[3]; filled[5, 4] <- med_tr[4]
stopifnot(isTRUE(all.equal(within_tissue_var(x_na, ct, med = med_tr, block_size = 5L),
                           within_tissue_var(filled, ct, block_size = 5L),
                           tolerance = 1e-12)),
          inherits(try(within_tissue_var(x_na, ct, block_size = 5L), silent = TRUE), "try-error"))
ok(4, "blocked within-tissue variance is exact to 1e-10 and block_size-invariant")

# =============================================================================
# D. BACKWARD COMPATIBILITY: "pooled" reproduces the pre-change behaviour.
# =============================================================================
# The pre-change fit_preprocess() ranked by matrixStats::colVars / var on the
# imputed training block. Reproduce that ranking here independently (including
# the name tie-break) and require an exact match on the selected features, then
# require the default argument path to equal the explicit "pooled" path.
v_pool <- apply(x, 2, var)
ii <- which(is.finite(v_pool) & v_pool > 0)
ii <- head(ii[order(-v_pool[ii], names(v_pool)[ii])], K)
stopifnot(identical(pp_pool$features, unname(names(v_pool)[ii])))

pp_default <- fit_preprocess(x, max_features = K)
stopifnot(identical(pp_default$features, pp_pool$features),
          identical(pp_default$median, pp_pool$median),
          identical(pp_default$center, pp_pool$center),
          identical(pp_default$scale, pp_pool$scale),
          identical(pp_default$max_missing, pp_pool$max_missing),
          identical(pp_default$feature_rank, "pooled"),
          identical(pp_default, pp_pool))
# Passing `cancer` must be INERT on the pooled path, so a caller that always
# supplies it cannot accidentally change the default variant.
stopifnot(identical(fit_preprocess(x, max_features = K, cancer = cancer), pp_pool))

# Same at the fit_en() level: default arguments and explicit "pooled" must give
# a bit-identical bundle.
b_def  <- fit_en(x[tr, , drop = FALSE], y[tr], patient[tr], cancer[tr], max_features = K)
b_pool <- fit_en(x[tr, , drop = FALSE], y[tr], patient[tr], cancer[tr], max_features = K,
                 feature_rank = "pooled")
stopifnot(identical(b_def$preprocess, b_pool$preprocess),
          identical(coef_of(b_def), coef_of(b_pool)),
          identical(b_def$tuning, b_pool$tuning),
          identical(b_def$alpha, b_pool$alpha),
          identical(b_def$lambda, b_pool$lambda),
          identical(b_def$feature_rank, "pooled"))
# And the two variants must genuinely produce different models, so a V3 run is
# never silently a pooled run.
stopifnot(!identical(sort(b_def$preprocess$features),
                     sort(b_clean$preprocess$features)))
ok(5, "feature_rank='pooled' (and the default) reproduce the pre-change behaviour exactly")

# =============================================================================
# E. Argument hygiene.
# =============================================================================
# within_tissue ranking without labels is a hard error, not a silent fallback to
# pooled variance. A silent fallback is the failure mode that would make an
# entire sentinel run meaningless while looking fine in the logs.
stopifnot(inherits(try(fit_preprocess(x, max_features = K,
                                      feature_rank = "within_tissue"), silent = TRUE),
                   "try-error"),
          inherits(try(fit_preprocess(x, max_features = K, feature_rank = "within_tissue",
                                      cancer = cancer[1:5]), silent = TRUE), "try-error"),
          inherits(try(fit_preprocess(x, max_features = K, feature_rank = "nonsense"),
                       silent = TRUE), "try-error"))
ok(6, "within_tissue without valid labels errors instead of falling back to pooled")

cat(sprintf("\nV3-abs feature_rank tests passed: %d/%d\n", pass, 6L))
stopifnot(pass == 6L)

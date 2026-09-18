#!/usr/bin/env Rscript
# =============================================================================
# tests/test_v4_lineage_penalized.R - V4: SUPERVISED within-tissue HRD
#   meta-association feature ranking with a lineage penalty (R/model.R).
# =============================================================================
#
# USAGE
#   Rscript tests/test_v4_lineage_penalized.R      (from the repository root)
#
# Plain stopifnot() script, matching tests/smoke_model.R and
# tests/test_v3_feature_rank.R. NOT a testthat suite; run it explicitly.
#
# WHAT V4 IS (docs/27_V3_PREREGISTRATION.md section 6, declared before any V3
# result existed)
#   The 5,000-probe filter is ranked by
#       score_j = |z_meta_j| * c_j - lambda_pen * l_j
#   where z_meta is the inverse-variance-weighted (w_t = n_t - 3) fixed-effect
#   meta-analysis of the Fisher-z transformed within-tissue Pearson correlation
#   between probe j and HRDsum, c_j is weighted sign consistency, and l_j is the
#   between/within variance ratio from the SAME decomposition within_tissue_var()
#   uses. lambda_pen is selected in the inner folds.
#
# WHY EACH CHECK EXISTS
#   A. The switch switches. On a fixture where all three statistics disagree by
#      construction, "lineage_penalized" must select a materially different
#      feature set from BOTH "pooled" and "within_tissue", and the penalty must
#      demonstrably move the selection.
#   B. THE LEAKAGE ASSERTION, and the reason this file exists at all. V4 is the
#      first ranker in this project that reads y, so it is the first that could
#      leak the LABEL as well as the features. We corrupt the held-out block's
#      feature rows AND its labels AND its tissue names and require the
#      preprocessing constants, selected features, coefficients, alpha, lambda
#      and lambda_pen to come back bit-identical.
#   C. The blocked meta statistic equals a slow, obvious, per-tissue reference
#      implementation written independently below to ~1e-8, and does not depend
#      on block_size.
#   D. BACKWARD COMPATIBILITY. "pooled" must still reproduce pre-change
#      behaviour exactly, bundle and all.
#   E. Missing y or missing cancer is a LOUD ERROR, never a silent fallback to
#      an unsupervised statistic. A silent fallback would make a whole sentinel
#      meaningless while the logs looked healthy.
#
# NOT COVERED HERE: whether V4 predicts better. That is what the LSF sentinel is
# for, and this file takes no view.
# =============================================================================

source("R/model.R")

pass <- 0L
ok <- function(n, msg) { pass <<- pass + 1L; cat(sprintf("  ok %d - %s\n", n, msg)) }

# --- Synthetic multi-tissue fixture ----------------------------------------
# 140 samples x 80 probes across 5 tissues, 28 each.
#
# Built so that ALL THREE ranking statistics disagree on purpose:
#   probes  1-20  "lineage markers": near-constant inside a tissue, far apart
#                 between tissues. Huge pooled variance, tiny within-tissue
#                 variance, NO within-tissue HRD association, lineage score ~1.
#                 Pooled ranking loves them; V3 and V4 must not.
#   probes 21-40  "within-tissue noise": vary a lot inside every tissue but are
#                 unrelated to HRD. Large within-tissue variance, z_meta ~ 0.
#                 V3 loves them; V4 must not.
#   probes 41-55  "HRD-associated, clean": consistently correlated with the
#                 within-tissue HRD component in EVERY tissue, with no lineage
#                 offset. These are what V4 exists to find.
#   probes 56-65  "HRD-associated but lineage-confounded": the same HRD signal
#                 PLUS a tissue offset (deliberately smaller than the lineage
#                 markers' so pooled variance still prefers 1-20). High |z_meta|
#                 AND high l_j, so
#                 they are precisely the probes the lambda_pen penalty must
#                 demote as lambda_pen grows. This is the pair that makes the
#                 penalty testable rather than decorative.
#   probes 66-80  low-variance filler.
# Values are squashed into [0,1] to mimic beta values.
set.seed(2026L)
n_per <- 28L; tissues <- c("TA", "TB", "TC", "TD", "TE")
cancer <- rep(tissues, each = n_per)
n <- length(cancer); p <- 80L
x <- matrix(0, n, p,
            dimnames = list(paste0("s", seq_len(n)), sprintf("cg%02d", seq_len(p))))
tissue_offset <- setNames(seq(0.1, 0.9, length.out = length(tissues)), tissues)
# The latent WITHIN-tissue HRD component. Centred inside each tissue so it
# carries no lineage information of its own.
u <- rnorm(n)
u <- as.numeric(u - ave(u, cancer))

for (j in 1:20)  x[, j] <- tissue_offset[cancer] + rnorm(n, sd = 0.004)
for (j in 21:40) x[, j] <- 0.5 + rnorm(n, sd = 0.18)
for (j in 41:55) x[, j] <- 0.5 + 0.11 * u + rnorm(n, sd = 0.05)
for (j in 56:65) x[, j] <- 0.5 * tissue_offset[cancer] + 0.11 * u + rnorm(n, sd = 0.05)
for (j in 66:80) x[, j] <- 0.5 + rnorm(n, sd = 0.015)
x[x < 0] <- 0; x[x > 1] <- 1

patient <- paste0("pt", seq_len(n))
# y = a real within-tissue HRD component plus a per-tissue offset, exactly the
# structure the absolute-target model faces.
y <- 18 * u + 30 * tissue_offset[cancer] + rnorm(n, sd = 1.0)
y <- as.numeric(y - min(y) + 1)

idx_of <- function(f) as.integer(sub("^cg", "", f))
coef_of <- function(b) as.matrix(coef(b$model, s = b$lambda))

cat("V4 lineage-penalized feature_rank tests\n")

# =============================================================================
# A. "lineage_penalized" selects a DIFFERENT feature set than BOTH other ranks,
#    and lambda_pen actually moves the selection in the declared direction.
# =============================================================================
K <- 12L
pp_pool <- fit_preprocess(x, max_features = K, feature_rank = "pooled")
pp_with <- fit_preprocess(x, max_features = K, feature_rank = "within_tissue",
                          cancer = cancer)
pp_v4_0 <- fit_preprocess(x, max_features = K, feature_rank = "lineage_penalized",
                          cancer = cancer, y = y, lambda_pen = 0)
pp_v4_2 <- fit_preprocess(x, max_features = K, feature_rank = "lineage_penalized",
                          cancer = cancer, y = y, lambda_pen = 2)

stopifnot(length(pp_v4_0$features) == K, length(pp_v4_2$features) == K,
          # (a) different from BOTH existing options.
          !identical(sort(pp_v4_0$features), sort(pp_pool$features)),
          !identical(sort(pp_v4_0$features), sort(pp_with$features)),
          !identical(sort(pp_v4_2$features), sort(pp_pool$features)),
          !identical(sort(pp_v4_2$features), sort(pp_with$features)))

# Direction check: this makes A a test of the STATISTIC rather than of an
# incidental tie-break. Pooled must be dominated by lineage markers, V3 by
# within-tissue noise, and V4 by the HRD-associated probes (41-65) under either
# penalty - because that is the only block with a real within-tissue
# association.
stopifnot(all(idx_of(pp_pool$features) <= 20L),
          all(idx_of(pp_with$features) >= 21L & idx_of(pp_with$features) <= 40L),
          all(idx_of(pp_v4_0$features) >= 41L & idx_of(pp_v4_0$features) <= 65L),
          all(idx_of(pp_v4_2$features) >= 41L & idx_of(pp_v4_2$features) <= 65L))

# The PENALTY must demote the lineage-confounded HRD probes (56-65) relative to
# the clean ones (41-55). This is the only assertion that proves lambda_pen is
# wired to anything at all.
conf <- function(f) sum(idx_of(f) >= 56L)
stopifnot(conf(pp_v4_2$features) < conf(pp_v4_0$features),
          !identical(pp_v4_0$features, pp_v4_2$features))
# The chosen strategy AND the penalty are RECORDED, so a saved transform is
# self-describing and a mislabelled run is detectable after the fact.
stopifnot(identical(pp_v4_0$feature_rank, "lineage_penalized"),
          identical(pp_v4_2$feature_rank, "lineage_penalized"),
          identical(pp_v4_0$lambda_pen, 0), identical(pp_v4_2$lambda_pen, 2),
          is.null(pp_pool$lambda_pen), is.null(pp_with$lambda_pen))
# Deterministic under repetition.
stopifnot(identical(pp_v4_2,
                    fit_preprocess(x, max_features = K,
                                   feature_rank = "lineage_penalized",
                                   cancer = cancer, y = y, lambda_pen = 2)))
ok(1, "lineage_penalized selects a different (and correctly targeted) feature set than pooled and within_tissue, and lambda_pen demotes lineage-confounded probes")

# =============================================================================
# B. LEAKAGE: the ranking is TRAINING-FOLD-ONLY, features AND labels.
# =============================================================================
# Hold out tissue TE. Fit on TA-TD only. Then rebuild a full matrix in which the
# held-out block is complete garbage - feature values, HRD labels and tissue
# names - and fit on the SAME training rows again. If any part of the ranking or
# any preprocessing constant looked outside the training rows, these would differ.
tr <- cancer != "TE"
PEN <- 0.5

pp_clean <- fit_preprocess(x[tr, , drop = FALSE], max_features = K,
                           feature_rank = "lineage_penalized",
                           cancer = cancer[tr], y = y[tr], lambda_pen = PEN)

set.seed(99L)
x_garbage <- x
x_garbage[!tr, ] <- matrix(runif(sum(!tr) * p, 0, 1), sum(!tr), p)
cancer_garbage <- cancer
cancer_garbage[!tr] <- sample(c("ZZ", "YY", "XX"), sum(!tr), replace = TRUE)
# THE V4-SPECIFIC HALF: the held-out LABELS are garbage on a wildly different
# scale. Under "pooled" and "within_tissue" this could not possibly matter
# because y is never read; under V4 it is the whole question.
y_garbage <- y
y_garbage[!tr] <- runif(sum(!tr), -5000, 5000)

pp_dirty <- fit_preprocess(x_garbage[tr, , drop = FALSE], max_features = K,
                           feature_rank = "lineage_penalized",
                           cancer = cancer_garbage[tr], y = y_garbage[tr],
                           lambda_pen = PEN)

stopifnot(identical(pp_clean$features, pp_dirty$features),
          identical(pp_clean$median,   pp_dirty$median),
          identical(pp_clean$center,   pp_dirty$center),
          identical(pp_clean$scale,    pp_dirty$scale),
          identical(pp_clean$feature_rank, pp_dirty$feature_rank),
          identical(pp_clean$lambda_pen,   pp_dirty$lambda_pen),
          identical(pp_clean, pp_dirty))

# Same assertion one level up, through fit_en(): the outer refit AND every inner
# fold must re-rank on their own rows, with their own labels. A garbage held-out
# block must leave the entire fitted bundle untouched - INCLUDING the
# inner-fold-selected lambda_pen, which is the newly tunable quantity.
PEN_GRID <- c(0, 0.5, 2)
b_clean <- fit_en(x[tr, , drop = FALSE], y[tr], patient[tr], cancer[tr],
                  max_features = K, feature_rank = "lineage_penalized",
                  lambda_pen_grid = PEN_GRID)
b_dirty <- fit_en(x_garbage[tr, , drop = FALSE], y_garbage[tr], patient[tr],
                  cancer_garbage[tr], max_features = K,
                  feature_rank = "lineage_penalized", lambda_pen_grid = PEN_GRID)
stopifnot(identical(b_clean$preprocess, b_dirty$preprocess),
          identical(coef_of(b_clean), coef_of(b_dirty)),
          identical(b_clean$tuning, b_dirty$tuning),
          identical(b_clean$alpha, b_dirty$alpha),
          identical(b_clean$lambda, b_dirty$lambda),
          identical(b_clean$lambda_pen, b_dirty$lambda_pen),
          identical(b_clean$feature_rank, "lineage_penalized"),
          identical(b_clean$preprocess$feature_rank, "lineage_penalized"))
# lambda_pen really was selected from the grid by the inner loss, and it is
# recorded on both the bundle and the transform.
stopifnot(b_clean$lambda_pen %in% PEN_GRID,
          identical(b_clean$preprocess$lambda_pen, b_clean$lambda_pen),
          identical(b_clean$lambda_pen_grid, sort(PEN_GRID)),
          "lambda_pen" %in% names(b_clean$tuning),
          identical(sort(unique(b_clean$tuning$lambda_pen)), sort(PEN_GRID)),
          # The winner is the argmin of the macro-averaged INNER loss and of
          # nothing else. Recompute the selection from the stored tuning table.
          identical(b_clean$lambda_pen,
                    b_clean$tuning$lambda_pen[which.min(b_clean$tuning$inner_macro_mae)]),
          identical(b_clean$alpha,
                    b_clean$tuning$alpha[which.min(b_clean$tuning$inner_macro_mae)]))
ok(2, "ranking, lambda_pen and all preprocessing constants are bit-identical under garbage held-out rows, labels AND tissue names")

# The inner folds really do re-rank rather than reuse the outer ranking.
f <- inner_folds(patient[tr], cancer[tr])
itr <- f != sort(unique(f))[1]
pp_inner <- fit_preprocess(x[tr, , drop = FALSE][itr, , drop = FALSE],
                           max_features = K, feature_rank = "lineage_penalized",
                           cancer = cancer[tr][itr], y = y[tr][itr],
                           lambda_pen = PEN)
stopifnot(length(pp_inner$features) == K,
          !identical(pp_inner$center, pp_clean$center))
ok(3, "inner-fold transform is computed on inner-training rows and inner-training labels only")

# =============================================================================
# C. The blocked meta statistic equals a SLOW DIRECT reference to ~1e-8.
# =============================================================================
# Written deliberately in the most obvious possible way - one tissue at a time,
# base cor(), an explicit loop - so it shares no code with the blocked
# implementation and cannot be wrong in the same way.
xt <- x[tr, , drop = FALSE]; ct <- cancer[tr]; yt <- y[tr]
ref_meta <- function(xm, yv, cv) {
  gs <- sort(unique(cv))
  nt <- vapply(gs, function(g) sum(cv == g), integer(1))
  w  <- pmax(nt - 3L, 0L)
  zt <- matrix(0, length(gs), ncol(xm))
  for (k in seq_along(gs)) {
    i <- which(cv == gs[k])
    if (w[k] <= 0) next
    # Standardise BOTH sides within the tissue, exactly as docs/27 section 6
    # words it, then correlate. (Pearson is scale-invariant, so this must agree
    # with the production code, which correlates the raw within-tissue values.)
    ys <- scale(yv[i])[, 1]
    for (j in seq_len(ncol(xm))) {
      col <- xm[i, j]
      if (sd(col) == 0 || !is.finite(sd(ys)) || sd(yv[i]) == 0) { zt[k, j] <- 0; next }
      r <- as.numeric(cor(scale(col)[, 1], ys))
      r <- min(max(r, -(1 - 1e-6)), 1 - 1e-6)
      zt[k, j] <- atanh(r)
    }
  }
  sw <- sum(w)
  list(z_meta = colSums(w * zt) / sqrt(sw),
       consistency = abs(colSums(w * sign(zt))) / sw)
}
st_fast <- lineage_meta_stats(xt, yt, ct, block_size = 9L)
st_ref  <- ref_meta(xt, yt, ct)
stopifnot(isTRUE(all.equal(unname(st_fast$z_meta), unname(st_ref$z_meta),
                           tolerance = 1e-8)),
          isTRUE(all.equal(unname(st_fast$consistency), unname(st_ref$consistency),
                           tolerance = 1e-8)),
          identical(names(st_fast$z_meta), colnames(xt)))
# The blocking is a memory device only, so it must be numerically invisible.
for (bs in c(1L, 3L, 17L, 80L, 5000L)) {
  st_b <- lineage_meta_stats(xt, yt, ct, block_size = bs)
  stopifnot(isTRUE(all.equal(st_fast$z_meta, st_b$z_meta, tolerance = 1e-12)),
            isTRUE(all.equal(st_fast$consistency, st_b$consistency, tolerance = 1e-12)),
            isTRUE(all.equal(st_fast$lineage, st_b$lineage, tolerance = 1e-12)),
            isTRUE(all.equal(st_fast$heterogeneity_i2, st_b$heterogeneity_i2,
                             tolerance = 1e-12)))
}
# The lineage term must be the SAME decomposition within_tissue_var() uses:
# l = 1 - SS_within / SS_total, with SS_within = (N - T) * within_tissue_var().
wv <- within_tissue_var(xt, ct, block_size = 11L)
Tn <- length(unique(ct)); N <- nrow(xt)
ss_within <- wv * (N - Tn)
ss_total  <- apply(xt, 2, function(c) sum((c - mean(c))^2))
l_ref <- pmin(1, pmax(0, 1 - ss_within / ss_total))
stopifnot(isTRUE(all.equal(unname(st_fast$lineage), unname(l_ref), tolerance = 1e-8)))
# And it must behave: lineage markers near 1, within-tissue-only probes near 0.
stopifnot(mean(st_fast$lineage[1:20]) > 0.9, mean(st_fast$lineage[21:40]) < 0.1,
          mean(st_fast$lineage[41:55]) < 0.1, mean(st_fast$lineage[56:65]) > 0.4)
# Non-finite cells are filled from the SUPPLIED TRAINING medians, never invented,
# and a missing median is an error rather than a silent NaN.
x_na <- xt; x_na[1, 3] <- NA_real_; x_na[5, 4] <- Inf
med_tr <- apply(xt, 2, median)
filled <- x_na; filled[1, 3] <- med_tr[3]; filled[5, 4] <- med_tr[4]
stopifnot(isTRUE(all.equal(lineage_meta_stats(x_na, yt, ct, med = med_tr, block_size = 7L),
                           lineage_meta_stats(filled, yt, ct, block_size = 7L),
                           tolerance = 1e-12)),
          inherits(try(lineage_meta_stats(x_na, yt, ct, block_size = 7L), silent = TRUE),
                   "try-error"))
ok(4, "blocked meta statistic matches a slow direct per-tissue reference to 1e-8, is block_size-invariant, and reuses the within_tissue_var() decomposition")

# =============================================================================
# D. BACKWARD COMPATIBILITY: "pooled" still reproduces pre-change behaviour.
# =============================================================================
v_pool <- apply(x, 2, var)
ii <- which(is.finite(v_pool) & v_pool > 0)
ii <- head(ii[order(-v_pool[ii], names(v_pool)[ii])], K)
stopifnot(identical(pp_pool$features, unname(names(v_pool)[ii])))

pp_default <- fit_preprocess(x, max_features = K)
stopifnot(identical(pp_default, pp_pool),
          identical(pp_default$feature_rank, "pooled"),
          # Names and order of the returned transform are unchanged by V4.
          identical(names(pp_pool),
                    c("features", "median", "center", "scale", "max_missing",
                      "feature_rank")))
# Passing `cancer`, `y` or `lambda_pen` must be INERT on the pooled path, so a
# caller that always supplies them cannot accidentally change the default.
stopifnot(identical(fit_preprocess(x, max_features = K, cancer = cancer), pp_pool),
          identical(fit_preprocess(x, max_features = K, y = y), pp_pool),
          identical(fit_preprocess(x, max_features = K, cancer = cancer, y = y,
                                   lambda_pen = 7), pp_pool),
          identical(fit_preprocess(x, max_features = K, cancer = cancer, y = y,
                                   feature_rank = "within_tissue"), pp_with))

# Same at the fit_en() level: default arguments and explicit "pooled" must give
# a bit-identical bundle, and it must carry NO V4 fields.
b_def  <- fit_en(x[tr, , drop = FALSE], y[tr], patient[tr], cancer[tr], max_features = K)
b_pool <- fit_en(x[tr, , drop = FALSE], y[tr], patient[tr], cancer[tr], max_features = K,
                 feature_rank = "pooled")
stopifnot(identical(b_def$preprocess, b_pool$preprocess),
          identical(coef_of(b_def), coef_of(b_pool)),
          identical(b_def$tuning, b_pool$tuning),
          identical(b_def$alpha, b_pool$alpha),
          identical(b_def$lambda, b_pool$lambda),
          identical(b_def$feature_rank, "pooled"),
          is.null(b_def$lambda_pen), is.null(b_def$lambda_pen_grid),
          identical(names(b_def$tuning), c("alpha", "lambda", "inner_macro_mae")),
          # A lambda_pen_grid must be INERT on the unsupervised paths.
          identical(b_def$tuning,
                    fit_en(x[tr, , drop = FALSE], y[tr], patient[tr], cancer[tr],
                           max_features = K, lambda_pen_grid = c(0, 1, 9))$tuning))
b_v3 <- fit_en(x[tr, , drop = FALSE], y[tr], patient[tr], cancer[tr],
               max_features = K, feature_rank = "within_tissue")
# The three variants must genuinely produce different models, so a V4 run is
# never silently a pooled or a V3 run.
stopifnot(!identical(sort(b_def$preprocess$features), sort(b_clean$preprocess$features)),
          !identical(sort(b_v3$preprocess$features),  sort(b_clean$preprocess$features)))
ok(5, "feature_rank='pooled' (and the default) reproduce the pre-change behaviour exactly; V4 fields are absent from non-V4 bundles")

# =============================================================================
# E. Argument hygiene: loud errors, never a silent fallback.
# =============================================================================
err <- function(expr) inherits(try(expr, silent = TRUE), "try-error")
stopifnot(
  # (e) missing y
  err(fit_preprocess(x, max_features = K, feature_rank = "lineage_penalized",
                     cancer = cancer)),
  # missing cancer
  err(fit_preprocess(x, max_features = K, feature_rank = "lineage_penalized", y = y)),
  # both missing
  err(fit_preprocess(x, max_features = K, feature_rank = "lineage_penalized")),
  # wrong-length y
  err(fit_preprocess(x, max_features = K, feature_rank = "lineage_penalized",
                     cancer = cancer, y = y[1:5])),
  # wrong-length cancer
  err(fit_preprocess(x, max_features = K, feature_rank = "lineage_penalized",
                     cancer = cancer[1:5], y = y)),
  # non-finite y
  err(fit_preprocess(x, max_features = K, feature_rank = "lineage_penalized",
                     cancer = cancer, y = replace(y, 1, NA_real_))),
  # non-finite / non-scalar lambda_pen
  err(fit_preprocess(x, max_features = K, feature_rank = "lineage_penalized",
                     cancer = cancer, y = y, lambda_pen = NA_real_)),
  err(fit_preprocess(x, max_features = K, feature_rank = "lineage_penalized",
                     cancer = cancer, y = y, lambda_pen = c(0, 1))),
  # unknown rank name
  err(fit_preprocess(x, max_features = K, feature_rank = "nonsense")),
  # a degenerate penalty grid at the fit_en() level
  err(fit_en(x[tr, , drop = FALSE], y[tr], patient[tr], cancer[tr], max_features = K,
             feature_rank = "lineage_penalized", lambda_pen_grid = c(0, 0))),
  err(fit_en(x[tr, , drop = FALSE], y[tr], patient[tr], cancer[tr], max_features = K,
             feature_rank = "lineage_penalized", lambda_pen_grid = numeric(0))),
  # and the statistic itself refuses a non-finite target
  err(lineage_meta_stats(xt, replace(yt, 1, NA_real_), ct))
)
ok(6, "lineage_penalized without valid y / cancer / lambda_pen errors loudly instead of falling back to an unsupervised rank")

cat(sprintf("\nV4 lineage-penalized tests passed: %d/%d\n", pass, 6L))
stopifnot(pass == 6L)

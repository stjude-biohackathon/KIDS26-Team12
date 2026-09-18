# =============================================================================
# tests/test_c1_rank_model.R - correctness + leakage checks for the
#                              tissue-relative (rank) model in R/rank_model.R.
# =============================================================================
#
# USAGE
#   Rscript tests/test_c1_rank_model.R        (from the repository root)
#
# Style matches tests/smoke_model.R: a plain stopifnot() script, NOT testthat.
# Everything runs on SMALL synthetic matrices (300 samples x 400 probes,
# 6 tissues) so the whole suite finishes in well under two minutes. The real
# 336k-probe matrix is never read.
#
# WHAT IS COVERED (numbered to match the cat() lines below)
#    1  GBM/LGG cannot enter fitting or summaries; the assertions ERROR when
#       locked samples are forced in.
#    2  Held-out tissue labels cannot affect target transforms or calibration.
#    3  Permuting held-out labels leaves predictions bit-identical.
#    4  Midrank behaviour of the normal-score target, including constant tissues.
#    5  weight_mode="tissue" gives every tissue equal TOTAL weight, sum(w) == n.
#    6  Single-sample scoring == batch scoring.
#    7  Column order does not change predictions.
#    8  Pure tissue intercept: no manufactured within-tissue performance.
#    9  Shared within-tissue signal transfers to an UNSEEN tissue.
#   10  Negative control: pure noise features give chance-level performance.
#   11  Pairwise ranker sanity + bounded/balanced/within-tissue pair sampling.
#   12  The tuning objective is a relative-space loss, not absolute HRDsum MAE.
# =============================================================================

source("R/rank_model.R")

ok <- function(n, msg) cat(sprintf("  [%2s] PASS  %s\n", n, msg))
fail_expected <- function(expr) inherits(try(expr, silent = TRUE), "try-error")

# -----------------------------------------------------------------------------
# Shared synthetic generators. Small, seeded, and deliberately boring.
# -----------------------------------------------------------------------------
# 6 tissues x 50 samples = 300 rows, 400 probes, one patient per row.
# The first `n_signal` probes are drawn BIMODAL (a 0.1/0.9 mixture, like a real
# CpG that is methylated in some samples and not in others) and therefore have
# roughly twice the variance of the uniform background probes. That matters:
# the pipeline's feature filter is UNSUPERVISED (top-variance), so a fixture in
# which the informative probes are indistinguishable from background by variance
# tests nothing but the filter's luck. Real informative CpGs are variable ones.
make_x <- function(seed, n = 300L, p = 400L, n_signal = 10L) {
  set.seed(seed)
  m <- matrix(runif(n * p), n, p,
              dimnames = list(paste0("s", seq_len(n)), paste0("cg", seq_len(p))))
  for (j in seq_len(n_signal)) {
    m[, j] <- ifelse(runif(n) < 0.5, runif(n, 0, 0.25), runif(n, 0.75, 1))
  }
  m
}
TISSUES <- paste0("T", 1:6)
cancer6 <- rep(TISSUES, each = 50L)
patients300 <- paste0("p", seq_len(300L))

# Scenario A: tissue intercept only, pure within-tissue noise (no real signal).
scenario_intercept <- function(seed) {
  x <- make_x(seed)
  set.seed(seed + 1L)
  offs <- setNames(seq(10, 60, length.out = 6), TISSUES)
  y <- as.numeric(offs[cancer6]) + rnorm(300, sd = 5)
  list(x = x, y = y, cancer = cancer6, patient = patients300)
}

# Scenario B: tissue intercept PLUS a within-tissue linear signal on the first
# 10 probes, with beta shared across every tissue. This is the situation in which
# a relative-target model should transfer to an unseen tissue.
scenario_signal <- function(seed) {
  x <- make_x(seed)
  set.seed(seed + 1L)
  offs <- setNames(seq(10, 60, length.out = 6), TISSUES)
  beta <- c(rep(c(8, -8), 5), rep(0, ncol(x) - 10L))
  y <- as.numeric(offs[cancer6]) + as.numeric(x %*% beta) + rnorm(300, sd = 2)
  list(x = x, y = y, cancer = cancer6, patient = patients300)
}

# Scenario C: y completely independent of x (negative control).
scenario_noise <- function(seed) {
  x <- make_x(seed)
  set.seed(seed + 2L)
  offs <- setNames(seq(10, 60, length.out = 6), TISSUES)
  y <- as.numeric(offs[cancer6]) + rnorm(300, sd = 5)
  # Deliberately break any x-y association that could survive by chance.
  y <- y[sample.int(300L)] + as.numeric(offs[cancer6])
  list(x = x, y = y, cancer = cancer6, patient = patients300)
}

# Run one leave-one-tissue-out fold and report held-out Spearman of the relative
# score against the held-out tissue's OWN relative target (EVAL-ONLY).
run_fold <- function(d, held_out = "T6", target_type = "normal_score",
                     weight_mode = "tissue", feature_rank = "pooled",
                     max_features = 120L) {
  tr <- d$cancer != held_out; te <- !tr
  b <- fit_rank_en(d$x[tr, , drop = FALSE], d$y[tr], d$patient[tr], d$cancer[tr],
                   target_type = target_type, weight_mode = weight_mode,
                   feature_rank = feature_rank, max_features = max_features,
                   block_size = 100L)
  p <- predict_rank(b, d$x[te, , drop = FALSE])
  truth <- as.numeric(relative_target(d$y[te], d$cancer[te], target_type))
  # A fully shrunk model emits a CONSTANT score: it declines to order the
  # held-out tissue at all. Spearman is undefined there, and the honest summary
  # of "no ordering claimed" is rho = 0, not NA. We record the fact separately so
  # a test can distinguish "no claim" from "a wrong claim".
  const <- sd(p$relative_score) == 0
  list(bundle = b, pred = p, truth = truth, constant_score = const,
       rho = if (const) 0 else suppressWarnings(cor(truth, p$relative_score, method = "spearman")))
}

cat("tests/test_c1_rank_model.R\n")

# =============================================================================
# 1. GBM/LGG are locked out of fitting and summaries.
# =============================================================================
# Build metadata that DOES contain GBM and LGG, then run the fold script's own
# masking logic and assert those rows land in neither train nor test. Then force
# them in and confirm the assertions ERROR rather than passing quietly.
meta_ct <- c(rep(TISSUES, each = 45L), rep("GBM", 15L), rep("LGG", 15L))
cns <- meta_ct %in% c("GBM", "LGG")
types <- sort(unique(meta_ct[!cns]))
stopifnot(!any(c("GBM", "LGG") %in% types))
type1 <- types[1]
tr <- !cns & meta_ct != type1
te <- !cns & meta_ct == type1
stopifnot(!any(meta_ct[tr] %in% c("GBM", "LGG")),
          !any(meta_ct[te] %in% c("GBM", "LGG")),
          !type1 %in% c("GBM", "LGG"),
          sum(tr | te) == sum(!cns),            # every non-CNS row used
          all(which(cns) %in% which(!(tr | te))))  # every CNS row dropped
stopifnot(isTRUE(assert_cns_locked(meta_ct[tr], meta_ct[te], type1)))
# 1b. Forcing GBM/LGG into either side must ERROR.
stopifnot(fail_expected(assert_cns_locked(c(meta_ct[tr], "GBM"), meta_ct[te], type1)),
          fail_expected(assert_cns_locked(meta_ct[tr], c(meta_ct[te], "LGG"), type1)),
          fail_expected(assert_cns_locked(meta_ct[tr], meta_ct[te], "GBM")))
# 1c. Even if the script-level guard were bypassed, the model itself refuses.
d1 <- scenario_signal(101L)
c_bad <- d1$cancer; c_bad[1:20] <- "GBM"
stopifnot(fail_expected(fit_rank_en(d1$x, d1$y, d1$patient, c_bad, max_features = 60L)),
          fail_expected(fit_pairwise_ranker(d1$x, d1$y, c_bad, max_features = 60L)))
ok(1, "GBM/LGG excluded from train and test; assertions fire when forced in")

# =============================================================================
# 2. Held-out tissue labels cannot influence the fit.
# =============================================================================
# Fit the SAME fold twice. In the second fit the held-out tissue's y is replaced
# by garbage. Because the held-out tissue never enters fit_rank_en(), the
# preprocess statistics, the coefficients and the calibration knots must be
# byte-identical.
d2 <- scenario_signal(202L)
tr2 <- d2$cancer != "T6"; te2 <- !tr2
fit2 <- function(yy) fit_rank_en(d2$x[tr2, , drop = FALSE], yy[tr2],
                                 d2$patient[tr2], d2$cancer[tr2],
                                 target_type = "normal_score", weight_mode = "tissue",
                                 max_features = 120L)
b_clean <- fit2(d2$y)
y_garbage <- d2$y
set.seed(2021L); y_garbage[te2] <- runif(sum(te2), -1e6, 1e6)
b_garbage <- fit2(y_garbage)
coef_of <- function(b) as.numeric(coef(b$model, s = b$lambda))
stopifnot(identical(b_clean$preprocess, b_garbage$preprocess),
          identical(coef_of(b_clean), coef_of(b_garbage)),
          identical(b_clean$calibration, b_garbage$calibration),
          identical(b_clean$alpha, b_garbage$alpha),
          identical(b_clean$lambda, b_garbage$lambda),
          identical(b_clean$tuning, b_garbage$tuning))
ok(2, "held-out labels do not change preprocess, coefficients or calibration knots")

# =============================================================================
# 3. Permutation invariance of predictions.
# =============================================================================
# Permuting the held-out tissue's labels cannot change its predictions, because
# predictions depend only on x and the frozen bundle. tolerance=0 => bit-identical.
p_clean <- predict_rank(b_clean, d2$x[te2, , drop = FALSE])
p_perm  <- predict_rank(b_garbage, d2$x[te2, , drop = FALSE])
stopifnot(isTRUE(all.equal(p_clean$relative_score, p_perm$relative_score, tolerance = 0)),
          isTRUE(all.equal(p_clean$percentile, p_perm$percentile, tolerance = 0)))
ok(3, "permuting held-out labels leaves predictions bit-identical")

# =============================================================================
# 4. Midrank behaviour of the normal-score target.
# =============================================================================
y4 <- c(1, 2, 2, 3,          # tissue A: one tie pair
        5, 5, 5, 5,          # tissue B: constant
        10, 20, 30, 40)      # tissue C: no ties
c4 <- rep(c("A", "B", "C"), each = 4L)
z4 <- relative_target(y4, c4, "normal_score")
q4 <- pnorm(z4)
# A: ranks 1, 2.5, 2.5, 4 over n=4 -> q = (r-0.5)/4
stopifnot(isTRUE(all.equal(q4[1:4], (c(1, 2.5, 2.5, 4) - 0.5) / 4, tolerance = 1e-12)),
          isTRUE(all.equal(q4[2], q4[3], tolerance = 1e-12)),   # ties share q
          isTRUE(all.equal(z4[2], z4[3], tolerance = 1e-12)),   # ties share z
          isTRUE(all.equal(q4[5:8], rep(0.5, 4), tolerance = 1e-12)),  # constant tissue
          isTRUE(all.equal(z4[5:8], rep(0, 4), tolerance = 1e-12)),
          isTRUE(all.equal(q4[9:12], (1:4 - 0.5) / 4, tolerance = 1e-12)),
          all(is.finite(z4)))
# Centered target: a constant tissue centres to exactly 0; each tissue sums to 0.
ct4 <- relative_target(y4, c4, "centered")
stopifnot(all(is.finite(ct4)), all(ct4[5:8] == 0),
          isTRUE(all.equal(as.numeric(tapply(ct4, c4, sum)), rep(0, 3), tolerance = 1e-12)))
# The recorded tissue statistics are the ones actually used.
ts <- attr(z4, "tissue_stats")
stopifnot(identical(ts$cancer_type, c("A", "B", "C")), ts$constant[2], !ts$constant[1])
ok(4, "midranks, ties, constant tissue (q=0.5, z=0) and finite z")

# =============================================================================
# 5. Tissue weights give equal TOTAL weight per tissue.
# =============================================================================
c5 <- c(rep("A", 10L), rep("B", 40L), rep("C", 150L))
w5 <- tissue_weights(c5, "tissue")
tot <- tapply(w5, c5, sum)
stopifnot(max(abs(tot - mean(tot))) < 1e-10,
          abs(sum(w5) - length(c5)) < 1e-10,
          all(tissue_weights(c5, "sample") == 1))
# The bundle's own refit weights obey it too.
w_b <- tissue_weights(b_clean$inner_folds$cancer_type, b_clean$weight_mode)
tot_b <- tapply(w_b, b_clean$inner_folds$cancer_type, sum)
stopifnot(max(abs(tot_b - mean(tot_b))) < 1e-10,
          abs(sum(w_b) - b_clean$training_n) < 1e-10)
ok(5, "weight_mode='tissue': equal total weight per tissue and sum(w) == n")

# =============================================================================
# 6. Single-sample scoring == batch scoring.
# =============================================================================
x_te <- d2$x[te2, , drop = FALSE]
p_one <- predict_rank(b_clean, x_te[1, , drop = FALSE])
p_all <- predict_rank(b_clean, x_te)
stopifnot(nrow(p_one) == 1L,
          isTRUE(all.equal(p_one$relative_score, p_all$relative_score[1], tolerance = 1e-12)),
          isTRUE(all.equal(p_one$percentile, p_all$percentile[1], tolerance = 1e-12)),
          is.finite(p_one$percentile), p_one$percentile > 0, p_one$percentile < 1)
# score_to_percentile is itself single-sample safe.
stopifnot(length(score_to_percentile(b_clean, p_one$relative_score)) == 1L)
ok(6, "one-row prediction equals its row in a batch; percentile needs no cohort")

# =============================================================================
# 7. Column order does not change predictions.
# =============================================================================
set.seed(707L)
perm <- sample(ncol(x_te))
p_shuf <- predict_rank(b_clean, x_te[, perm, drop = FALSE])
stopifnot(isTRUE(all.equal(p_all$relative_score, p_shuf$relative_score, tolerance = 0)),
          isTRUE(all.equal(p_all$percentile, p_shuf$percentile, tolerance = 0)))
# Also with an extra unseen probe present and one row only.
x_extra <- cbind(x_te[1, , drop = FALSE], cg_unseen = 0.5)[, sample(ncol(x_te) + 1L), drop = FALSE]
stopifnot(isTRUE(all.equal(predict_rank(b_clean, x_extra)$relative_score,
                           p_all$relative_score[1], tolerance = 1e-12)))
ok(7, "column permutation (and unseen probes) leave predictions identical")

# =============================================================================
# 8. Pure tissue intercept -> no manufactured within-tissue performance.
# =============================================================================
# y = tissue intercept + pure noise, so there is NOTHING within tissue to learn.
# Seed 808. A relative-target model must not invent a within-tissue ordering.
d8 <- scenario_intercept(808L)
r8_v1 <- suppressWarnings(run_fold(d8, weight_mode = "sample"))
r8_v2 <- suppressWarnings(run_fold(d8, weight_mode = "tissue"))
cat(sprintf("       scenario intercept-only (seed 808): V1 rho=%.3f (constant score: %s)  V2 rho=%.3f (constant score: %s)\n",
            r8_v1$rho, r8_v1$constant_score, r8_v2$rho, r8_v2$constant_score))
stopifnot(is.finite(r8_v1$rho), is.finite(r8_v2$rho),
          abs(r8_v1$rho) < 0.25, abs(r8_v2$rho) < 0.25)
ok(8, "pure tissue intercept: |held-out rho| < 0.25 (seed 808)")

# =============================================================================
# 9. Shared within-tissue signal transfers to an UNSEEN tissue.
# =============================================================================
# Seed 909. beta is shared across tissues, so the ordering learned on T1..T5
# should hold on the never-seen T6.
d9 <- scenario_signal(909L)
r9_v1 <- run_fold(d9, weight_mode = "sample", feature_rank = "pooled")
r9_v2 <- run_fold(d9, weight_mode = "tissue", feature_rank = "pooled")
r9_v3 <- run_fold(d9, weight_mode = "tissue", feature_rank = "within_tissue")
r9_ct <- run_fold(d9, weight_mode = "tissue", target_type = "centered")
cat(sprintf("       scenario signal (seed 909): V1 rho=%.3f  V2 rho=%.3f  V3 rho=%.3f  centered rho=%.3f\n",
            r9_v1$rho, r9_v2$rho, r9_v3$rho, r9_ct$rho))
stopifnot(r9_v1$rho > 0.4, r9_v2$rho > 0.4, r9_v3$rho > 0.4, r9_ct$rho > 0.4)
# V3 really did use a within-tissue variance ranking (different feature set).
stopifnot(identical(r9_v3$bundle$feature_rank, "within_tissue"),
          !identical(sort(r9_v3$bundle$preprocess$features),
                     sort(r9_v2$bundle$preprocess$features)))
# Blocked within-tissue variance must equal a direct (small-scale) computation.
xx <- d9$x[d9$cancer != "T6", 1:40, drop = FALSE]; cc <- d9$cancer[d9$cancer != "T6"]
direct <- apply(xx, 2, function(col) sum(tapply(col, cc, function(v) sum((v - mean(v))^2)))) /
  (nrow(xx) - length(unique(cc)))
stopifnot(isTRUE(all.equal(unname(within_tissue_var(xx, cc, block_size = 7L)),
                           unname(direct), tolerance = 1e-10)),
          isTRUE(all.equal(within_tissue_var(xx, cc, block_size = 7L),
                           within_tissue_var(xx, cc, block_size = 40L), tolerance = 1e-12)))
ok(9, "within-tissue signal transfers to an unseen tissue (rho > 0.4) for V1/V2/V3")

# =============================================================================
# 10. Negative control: features unrelated to y.
# =============================================================================
d10 <- scenario_noise(1010L)
r10 <- run_fold(d10, weight_mode = "tissue")
np10 <- rank_null_panel(r10$truth, r10$pred$relative_score,
                        rep("T6", length(r10$truth)), n_perm = 500L)
cat(sprintf("       negative control (seed 1010): rho=%.3f  perm p=%.3f  null 95%% [%.3f, %.3f]\n",
            r10$rho, np10$perm_p_value, np10$null_spearman_q025, np10$null_spearman_q975))
stopifnot(abs(r10$rho) < 0.25,
          np10$perm_p_value > 0.05,
          r10$rho >= np10$null_spearman_q025, r10$rho <= np10$null_spearman_q975)
ok(10, "negative control indistinguishable from the within-tissue permutation null")

# =============================================================================
# 11. Pairwise ranker prototype (Phase D).
# =============================================================================
# 11a. Pair sampling: bounded, balanced, strictly within tissue.
d11 <- scenario_signal(1111L)
prs <- sample_within_tissue_pairs(d11$y, d11$cancer, max_pairs_per_tissue = 60L, seed = 11L)
cnt <- table(prs$cancer_type)
stopifnot(all(cnt <= 60L), length(cnt) == 6L,
          all(d11$cancer[prs$i] == d11$cancer[prs$j]),
          all(prs$i != prs$j),
          all(prs$outcome == as.integer(d11$y[prs$i] > d11$y[prs$j])))
bal <- tapply(prs$outcome, prs$cancer_type, mean)
stopifnot(mean(prs$outcome) >= 0.4, mean(prs$outcome) <= 0.6,
          all(bal >= 0.4 & bal <= 0.6))
# A tissue with fewer than two distinct y values supplies no pairs.
y_const <- d11$y; y_const[d11$cancer == "T3"] <- 7
stopifnot(!"T3" %in% sample_within_tissue_pairs(y_const, d11$cancer, 60L, seed = 11L)$cancer_type)
# 11b. Recovers the within-tissue ordering on the signal scenario, ~zero on noise.
pw_fold <- function(d, held_out = "T6") {
  tr <- d$cancer != held_out; te <- !tr
  b <- fit_pairwise_ranker(d$x[tr, , drop = FALSE], d$y[tr], d$cancer[tr],
                           max_features = 120L, max_pairs_per_tissue = 300L, seed = 11L)
  p <- predict_pairwise(b, d$x[te, , drop = FALSE])
  truth <- as.numeric(relative_target(d$y[te], d$cancer[te], "normal_score"))
  list(b = b, p = p, rho = suppressWarnings(cor(truth, p$relative_score, method = "spearman")))
}
pw_sig <- pw_fold(d11)
pw_noise <- pw_fold(scenario_noise(1112L))
cat(sprintf("       pairwise ranker: signal rho=%.3f   noise rho=%.3f\n",
            pw_sig$rho, pw_noise$rho))
stopifnot(pw_sig$rho > 0.4, abs(pw_noise$rho) < 0.25,
          identical(pw_sig$b$score_type, "relative_score"),
          nrow(predict_pairwise(pw_sig$b, d11$x[1, , drop = FALSE])) == 1L,
          isTRUE(all.equal(predict_pairwise(pw_sig$b, d11$x[which(d11$cancer == "T6")[1], , drop = FALSE])$relative_score,
                           pw_sig$p$relative_score[1], tolerance = 1e-12)))
ok(11, "pairwise ranker: bounded/balanced/within-tissue pairs; signal recovered, noise not")

# =============================================================================
# 12. The tuning objective is a relative-space loss.
# =============================================================================
stopifnot(identical(RELATIVE_LOSS_COLUMN, "inner_macro_relative_MAE"),
          identical(b_clean$tuning_objective, "inner_macro_relative_MAE"),
          RELATIVE_LOSS_COLUMN %in% names(b_clean$tuning),
          "inner_macro_one_minus_spearman" %in% names(b_clean$tuning),
          !any(grepl("HRDsum", names(b_clean$tuning), ignore.case = TRUE)),
          !any(grepl("absolute", names(b_clean$tuning), ignore.case = TRUE)),
          # The selected row really is the argmin of the relative-space column.
          which.min(b_clean$tuning[[RELATIVE_LOSS_COLUMN]]) ==
            which(b_clean$tuning$alpha == b_clean$alpha & b_clean$tuning$lambda == b_clean$lambda)[1])
# The bundle carries nothing that would permit back-transformation to HRDsum.
stopifnot(identical(b_clean$score_type, "relative_score"),
          !any(c("training_mean", "tissue_means", "y_center", "y_scale",
                 "target_transform") %in% names(b_clean)),
          !any(grepl("HRDsum", unlist(lapply(b_clean, function(e)
            if (is.character(e)) e else character(0))), fixed = TRUE)))
ok(12, "tuning objective is relative-space; bundle cannot back-transform to HRDsum")

cat("\nALL 12 TEST CASES PASSED - R/rank_model.R\n")
cat("Relative-target model: locked CNS partition enforced, held-out labels inert,\n")
cat("  midranks correct, tissue weights balanced, single-sample scoring exact,\n")
cat("  signal transfers to an unseen tissue and noise does not.\n")

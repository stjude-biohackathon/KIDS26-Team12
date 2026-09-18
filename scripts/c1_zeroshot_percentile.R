#!/usr/bin/env Rscript
# =============================================================================
# c1_zeroshot_percentile.R
#
# PURPOSE
#   Phase B of investigation C1. Zero-shot percentile triage: WITHOUT any
#   target-tissue label or target-tissue cohort statistic, can the raw
#   prediction be mapped through a frozen source-tissue-only mapping to a
#   within-tissue percentile q_hat that beats the constant-0.5 null?
#
#   Leave-one-tissue-out over all 30 tissues. For each target tissue t the
#   mapping f: raw prediction -> percentile is fit on the other 29 tissues only,
#   with equal TOTAL weight per source tissue, in two variants:
#     (i)  LINEAR   : weighted lm(q ~ prediction), clipped to [0,1]
#     (ii) MONOTONE : weighted isotonic regression (hand-rolled weighted PAVA)
#   A third frozen source-only mapping (weighted logistic regression) gives
#   calibrated probabilities for the top-quartile event, scored by Brier.
#
# INPUT (the ONLY data file read by this script)
#   results/loco_run01/loco_predictions.tsv
#
# OUTPUTS (results/c1_rank_probe/)
#   zeroshot_percentile_by_tissue.tsv
#   zeroshot_percentile_macro.tsv
#   zeroshot_reliability_bins.tsv
#   zeroshot_tissue_r2.tsv
#
# LEAKAGE BOUNDARIES
#   * The mapping for target tissue t uses source-tissue rows ONLY; tissue t
#     contributes nothing (no labels, no cohort mean, no rank, no n).
#   * A sample's q_hat depends only on its own raw prediction and the frozen
#     mapping. This is asserted: one-row-at-a-time scoring must equal block
#     scoring (stopifnot + all.equal).
#   * Tissue t's own labels are used for EVALUATION ONLY (true q, AUROC, Brier).
#   * GBM and LGG are excluded by construction (asserted below).
#   * No pediatric / PBTA / PBTP data is touched.
#
# !! STACKED-DIAGNOSTIC WARNING !!
#   The `predicted_reference_HRDsum` column is itself a LOCO prediction, i.e. it
#   comes from 30 DIFFERENT outer-fold models. Stacking a percentile mapping on
#   top of heterogeneous outer-fold predictions means this is a FEASIBILITY
#   SCREEN ONLY, not a validated frozen-model result. A deployable version would
#   require a single frozen genomic model plus a mapping fit inside the same
#   nested CV structure.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

set.seed(4871)

IN_FILE <- "results/loco_run01/loco_predictions.tsv"
OUT_DIR <- "results/c1_rank_probe"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("\n*** STACKED DIAGNOSTIC WARNING ***\n",
    "Input predictions come from 30 different LOCO outer-fold models. The\n",
    "percentile mapping is stacked on top of those heterogeneous predictions,\n",
    "so these results are a FEASIBILITY SCREEN ONLY -- not a validated\n",
    "frozen-model result.\n\n", sep = "")

# ---------------------------------------------------------------- load + guard
dat <- fread(IN_FILE, sep = "\t")
stopifnot(!any(dat$cancer_type %in% c("GBM", "LGG")))
stopifnot(!any(is.na(dat$predicted_reference_HRDsum)), !any(is.na(dat$actual)))

tissues <- sort(unique(dat$cancer_type))

# within-tissue midrank percentile of the TRUE label
dat[, q_true := (frank(actual, ties.method = "average") - 0.5) / .N, by = cancer_type]
# equal total weight per tissue (renormalised per fit below)
dat[, wt_raw := 1 / .N, by = cancer_type]

# ------------------------------------------------------- weighted PAVA (isotonic)
# Returns block boundaries + fitted values for a non-decreasing weighted fit.
weighted_pava <- function(x, y, w) {
  stopifnot(length(x) == length(y), length(y) == length(w), all(w > 0))
  o <- order(x)
  x <- x[o]; y <- y[o]; w <- w[o]
  # pool exact x ties first so the fit is a well-defined function of x
  if (anyDuplicated(x)) {
    agg <- data.table(x = x, y = y, w = w)[, .(y = sum(y * w) / sum(w), w = sum(w)),
                                           by = x][order(x)]
    x <- agg$x; y <- agg$y; w <- agg$w
  }
  n <- length(x)
  bv <- numeric(n); bw <- numeric(n); bx <- numeric(n); bn <- integer(n)
  j <- 0L
  for (i in seq_len(n)) {
    j <- j + 1L
    bv[j] <- y[i]; bw[j] <- w[i]; bx[j] <- x[i]; bn[j] <- 1L
    while (j > 1L && bv[j - 1L] > bv[j]) {          # violation -> pool
      nw <- bw[j - 1L] + bw[j]
      bv[j - 1L] <- (bv[j - 1L] * bw[j - 1L] + bv[j] * bw[j]) / nw
      bw[j - 1L] <- nw
      bx[j - 1L] <- bx[j]                           # right edge of pooled block
      bn[j - 1L] <- bn[j - 1L] + bn[j]
      j <- j - 1L
    }
  }
  fitted <- rep(bv[seq_len(j)], bn[seq_len(j)])     # per (unique) x, in x order
  list(x = x, fitted = fitted,
       block_x = bx[seq_len(j)], block_value = bv[seq_len(j)])
}

# --- unit check: with equal weights the PAVA must reproduce stats::isoreg ------
local({
  set.seed(1)
  xx <- sort(runif(60)); yy <- xx + rnorm(60, sd = 0.4)
  p  <- weighted_pava(xx, yy, rep(1, 60))
  ir <- stats::isoreg(xx, yy)
  stopifnot(isTRUE(all.equal(as.numeric(p$fitted), as.numeric(ir$yf),
                             tolerance = 1e-10)))
  # monotonicity guarantee
  stopifnot(all(diff(p$fitted) >= -1e-12))
})
cat("PAVA unit check vs stats::isoreg (equal weights): PASS\n")

pava_predict <- function(fit, xout) {
  as.numeric(approx(x = fit$x, y = fit$fitted, xout = xout,
                    method = "linear", rule = 2, ties = "ordered")$y)
}

# --------------------------------------------------------------------- AUROC
auroc_mw <- function(score, event) {
  pos <- score[event]; neg <- score[!event]
  if (length(pos) < 5L || length(neg) < 5L) return(NA_real_)
  r <- rank(c(pos, neg))
  (sum(r[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) /
    (length(pos) * length(neg))
}

clip01 <- function(x) pmin(1, pmax(0, x))

# ------------------------------------------------------------- main LOTO loop
per_tissue <- list(); pooled <- list()

for (t in tissues) {
  src <- dat[cancer_type != t]
  tgt <- dat[cancer_type == t]
  w   <- src$wt_raw / sum(src$wt_raw)

  # (i) linear
  fit_lin <- lm(q_true ~ predicted_reference_HRDsum, data = src, weights = w)
  # (ii) monotone
  fit_iso <- weighted_pava(src$predicted_reference_HRDsum, src$q_true, w)
  # (iii) frozen probability model for the top-quartile event
  fit_glm <- suppressWarnings(
    glm(I(q_true > 0.75) ~ predicted_reference_HRDsum, data = src,
        family = binomial, weights = w))

  score_block <- function(x) {
    list(lin = clip01(as.numeric(predict(fit_lin,
            newdata = data.frame(predicted_reference_HRDsum = x)))),
         iso = clip01(pava_predict(fit_iso, x)),
         prb = as.numeric(predict(fit_glm,
            newdata = data.frame(predicted_reference_HRDsum = x),
            type = "response")))
  }

  blk <- score_block(tgt$predicted_reference_HRDsum)
  # CRITICAL leakage assertion: per-sample scoring == block scoring
  one <- lapply(tgt$predicted_reference_HRDsum, score_block)
  stopifnot(isTRUE(all.equal(blk$lin, sapply(one, `[[`, "lin"))),
            isTRUE(all.equal(blk$iso, sapply(one, `[[`, "iso"))),
            isTRUE(all.equal(blk$prb, sapply(one, `[[`, "prb"))))

  qt   <- tgt$q_true                      # EVALUATION ONLY
  ev   <- qt > 0.75
  mae  <- function(a, b) mean(abs(a - b))

  cal <- function(qh) {
    if (sd(qh) < 1e-12) return(c(NA_real_, NA_real_))
    cf <- coef(lm(qt ~ qh)); c(cf[1], cf[2])
  }
  cl <- cal(blk$lin); ci <- cal(blk$iso)

  per_tissue[[t]] <- data.table(
    cancer_type = t, n = nrow(tgt),
    mae_q_linear   = mae(blk$lin, qt),
    mae_q_monotone = mae(blk$iso, qt),
    mae_q_null_0.5 = mae(rep(0.5, nrow(tgt)), qt),
    spearman_linear   = suppressWarnings(cor(blk$lin, qt, method = "spearman")),
    spearman_monotone = suppressWarnings(cor(blk$iso, qt, method = "spearman")),
    spearman_raw_vs_actual = suppressWarnings(
      cor(tgt$predicted_reference_HRDsum, tgt$actual, method = "spearman")),
    n_events_top_quartile = sum(ev), n_nonevents = sum(!ev),
    auroc_linear   = auroc_mw(blk$lin, ev),
    auroc_monotone = auroc_mw(blk$iso, ev),
    brier_frozen_logistic = mean((blk$prb - as.numeric(ev))^2),
    brier_null_0.25       = mean((0.25 - as.numeric(ev))^2),
    calib_intercept_linear = cl[1], calib_slope_linear = cl[2],
    calib_intercept_monotone = ci[1], calib_slope_monotone = ci[2],
    fail_linear_vs_null   = mae(blk$lin, qt) > mae(rep(0.5, nrow(tgt)), qt),
    fail_monotone_vs_null = mae(blk$iso, qt) > mae(rep(0.5, nrow(tgt)), qt))

  pooled[[t]] <- data.table(cancer_type = t,
                            raw_pred = tgt$predicted_reference_HRDsum,
                            q_hat_linear = blk$lin, q_hat_monotone = blk$iso,
                            p_hat = blk$prb, q_true = qt, event = ev,
                            wt = 1 / nrow(tgt))
}

by_tissue <- rbindlist(per_tissue)[order(cancer_type)]
pool      <- rbindlist(pooled)

write.table(by_tissue, file.path(OUT_DIR, "zeroshot_percentile_by_tissue.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------------------------------------------ macro
pc <- function(f) coef(lm(q_true ~ qh, data = data.frame(q_true = pool$q_true, qh = f)))
cl_p <- pc(pool$q_hat_linear); ci_p <- pc(pool$q_hat_monotone)

macro <- data.table(
  n_tissues = nrow(by_tissue),
  macro_mae_q_linear   = mean(by_tissue$mae_q_linear),
  macro_mae_q_monotone = mean(by_tissue$mae_q_monotone),
  macro_mae_q_null_0.5 = mean(by_tissue$mae_q_null_0.5),
  n_tissues_linear_beats_null   = sum(!by_tissue$fail_linear_vs_null),
  n_tissues_monotone_beats_null = sum(!by_tissue$fail_monotone_vs_null),
  macro_spearman_linear   = mean(by_tissue$spearman_linear),
  macro_spearman_monotone = mean(by_tissue$spearman_monotone),
  macro_spearman_raw      = mean(by_tissue$spearman_raw_vs_actual),
  max_abs_spearman_gap_linear = max(abs(by_tissue$spearman_linear -
                                        by_tissue$spearman_raw_vs_actual)),
  max_abs_spearman_gap_monotone = max(abs(by_tissue$spearman_monotone -
                                          by_tissue$spearman_raw_vs_actual)),
  macro_auroc_linear   = mean(by_tissue$auroc_linear, na.rm = TRUE),
  macro_auroc_monotone = mean(by_tissue$auroc_monotone, na.rm = TRUE),
  n_tissues_auroc_evaluable = sum(!is.na(by_tissue$auroc_monotone)),
  macro_brier_frozen_logistic = mean(by_tissue$brier_frozen_logistic),
  macro_brier_null_0.25       = mean(by_tissue$brier_null_0.25),
  n_tissues_brier_beats_null  = sum(by_tissue$brier_frozen_logistic <
                                    by_tissue$brier_null_0.25),
  macro_calib_intercept_linear   = mean(by_tissue$calib_intercept_linear, na.rm = TRUE),
  macro_calib_slope_linear       = mean(by_tissue$calib_slope_linear, na.rm = TRUE),
  macro_calib_intercept_monotone = mean(by_tissue$calib_intercept_monotone, na.rm = TRUE),
  macro_calib_slope_monotone     = mean(by_tissue$calib_slope_monotone, na.rm = TRUE),
  pooled_calib_intercept_linear   = cl_p[1], pooled_calib_slope_linear   = cl_p[2],
  pooled_calib_intercept_monotone = ci_p[1], pooled_calib_slope_monotone = ci_p[2])
write.table(macro, file.path(OUT_DIR, "zeroshot_percentile_macro.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# -------------------------------------------------------- reliability bins
brks <- seq(0, 1, by = 0.1)
relbin <- function(qh, lab) {
  b <- cut(qh, breaks = brks, include.lowest = TRUE)
  data.table(method = lab, bin = levels(b),
             bin_lower = brks[-11], bin_upper = brks[-1])[
    data.table(bin = as.character(b), q_hat = qh, q_true = pool$q_true)[
      , .(n = .N, mean_q_hat = mean(q_hat), mean_q_true = mean(q_true)),
      by = bin], on = "bin"]
}
bins <- rbindlist(list(relbin(pool$q_hat_linear, "linear"),
                       relbin(pool$q_hat_monotone, "monotone")),
                  use.names = TRUE)
bins <- bins[order(method, bin_lower)]
write.table(bins, file.path(OUT_DIR, "zeroshot_reliability_bins.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------------- tissue-identity R^2 on the pooled set
r2 <- function(y) summary(lm(y ~ factor(pool$cancer_type)))$r.squared
tr2 <- data.table(
  quantity = c("q_hat_linear", "q_hat_monotone", "raw_predicted_reference_HRDsum",
               "q_true"),
  r2_tissue_identity = c(r2(pool$q_hat_linear), r2(pool$q_hat_monotone),
                         r2(pool$raw_pred), r2(pool$q_true)),
  note = c("frozen source-only linear mapping",
           "frozen source-only monotone mapping",
           "apples-to-apples comparator for the reported 56.2% figure",
           "by construction ~0: within-tissue percentile removes tissue level"))
write.table(tr2, file.path(OUT_DIR, "zeroshot_tissue_r2.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------------------------------------------ stdout
cat("\n=== PHASE B macro summary (equal weight per tissue, 30 tissues) ===\n")
print(t(round(as.matrix(macro[, !"n_tissues"]), 4)))

cat("\n=== Worst tissues by monotone MAE vs 0.5 null (top 8) ===\n")
print(head(by_tissue[order(-(mae_q_monotone - mae_q_null_0.5)),
                     .(cancer_type, n, mae_mono = round(mae_q_monotone, 3),
                       mae_null = round(mae_q_null_0.5, 3),
                       auroc = round(auroc_monotone, 3),
                       rho = round(spearman_monotone, 3))], 8), row.names = FALSE)

cat("\n=== Tissue-identity R^2 on pooled held-out rows ===\n")
print(tr2[, .(quantity, r2 = round(r2_tissue_identity, 4))], row.names = FALSE)

cat("\n=== Reliability bins (monotone) ===\n")
print(bins[method == "monotone", .(bin, n, mean_q_hat = round(mean_q_hat, 3),
                                   mean_q_true = round(mean_q_true, 3))],
      row.names = FALSE)

cat("\nCAVEAT: both mappings are monotone (non-decreasing) in the raw prediction,\n",
    "so they cannot change within-tissue ranking. For the STRICTLY increasing\n",
    "linear map, Spearman(q_hat, q_true) EQUALS Spearman(raw prediction, actual)\n",
    "within each tissue exactly (max abs gap = ",
    signif(macro$max_abs_spearman_gap_linear, 3), "). The isotonic map is only\n",
    "WEAKLY increasing -- its flat/clipped blocks create ties, which can only\n",
    "LOSE ranking information (max abs gap = ",
    signif(macro$max_abs_spearman_gap_monotone, 3), ", always downward). Either\n",
    "way the mapping cannot improve within-tissue ranking, only its absolute\n",
    "placement on the percentile scale.\n", sep = "")
cat("CAVEAT: STACKED DIAGNOSTIC -- 30 different outer-fold models feed these\n",
    "predictions; feasibility screen only.\n", sep = "")
cat("\nWrote outputs to ", OUT_DIR, "\n", sep = "")

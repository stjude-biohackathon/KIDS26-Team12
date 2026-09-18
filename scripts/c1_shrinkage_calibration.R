#!/usr/bin/env Rscript
# =============================================================================
# c1_shrinkage_calibration.R
#
# PURPOSE
#   Phase A of investigation C1. Tests hypothesis H1: empirical-Bayes shrinkage
#   of the few-shot tissue offset toward a cross-tissue prior makes small-k
#   calibration SAFE (i.e. never worse than the uncorrected model), whereas the
#   existing unshrunk few-shot calibrator is harmful at k=3.
#
#   For each held-out target tissue t and each k in {1,3,5,10,20}, we draw B=1000
#   random k-sample "shot" sets, estimate the tissue offset, shrink it toward a
#   source-tissue-only prior, and evaluate MAE on the remaining ("eval") rows
#   against six comparators, including a FAIR HIERARCHICAL NULL (m6) that is
#   given exactly the same prior information as the shrunk model.
#
# INPUT (the ONLY data file read by this script)
#   results/loco_run01/loco_predictions.tsv
#     leave-one-cancer-type-out predictions: each sample scored by a model
#     trained on the other 29 tissues.
#
# OUTPUTS (results/c1_shrinkage/)
#   shrinkage_by_tissue_k.tsv        per (cancer_type, k) summary
#   shrinkage_macro_by_k.tsv         macro-over-tissue summary per k
#   shrinkage_draw_distribution.tsv  per-draw delta quantiles (+ macro pool)
#   shrinkage_by_offset_stratum.tsv  macro summary stratified by |true offset|
#   shrinkage_priors.tsv             leave-one-tissue-out prior parameters
#
# LEAKAGE BOUNDARIES
#   * Priors (mu0, tau2, sigma2, mu0_y, tau2_y, sigma2_y) for target tissue t are
#     estimated from SOURCE tissues only (all tissues except t). Tissue t
#     contributes nothing to its own prior.
#   * The k "shot" rows are the ONLY target-tissue labels used to form an offset.
#   * Evaluation uses only the held-out eval rows (disjoint from the shot set).
#   * b_oracle (full-tissue offset) and the true-offset column in
#     shrinkage_priors.tsv are BENCHMARK / REFERENCE ONLY and never feed a fit.
#   * GBM and LGG are excluded by construction (asserted below).
#   * No pediatric / PBTA / PBTP data is touched.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

set.seed(4871)

IN_FILE  <- "results/loco_run01/loco_predictions.tsv"
OUT_DIR  <- "results/c1_shrinkage"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

K_GRID <- c(1, 3, 5, 10, 20)
B      <- 1000L
MIN_EVAL <- 10L

# ---------------------------------------------------------------- load + guard
dat <- fread(IN_FILE, sep = "\t")

stopifnot(!any(dat$cancer_type %in% c("GBM", "LGG")))
stopifnot(nrow(dat) == 7065L)
stopifnot(all(c("predicted_reference_HRDsum", "actual", "cancer_type") %in% names(dat)))
stopifnot(!any(is.na(dat$predicted_reference_HRDsum)), !any(is.na(dat$actual)))

dat[, resid := predicted_reference_HRDsum - actual]

tissues <- sort(unique(dat$cancer_type))
cat(sprintf("Loaded %d rows, %d cancer types.\n", nrow(dat), length(tissues)))

# per-tissue sufficient statistics (used only to build LOO priors)
ts <- dat[, .(n   = .N,
              rbar = mean(resid),
              s2   = var(resid),
              ybar = mean(actual),
              s2y  = var(actual)), by = cancer_type]
stopifnot(!any(is.na(ts$s2)), !any(is.na(ts$s2y)))   # every tissue has n >= 2

# ------------------------------------------------------- LOO prior computation
prior_for <- function(t) {
  s <- ts[cancer_type != t]
  mu0      <- mean(s$rbar)                    # equal weight per source tissue
  sigma2   <- mean(s$s2)                      # macro pooled within-tissue var
  tau2_raw <- var(s$rbar)                     # equal-weight across-tissue var
  tau2     <- max(0, tau2_raw - mean(s$s2 / s$n))
  mu0_y    <- mean(s$ybar)
  sigma2_y <- mean(s$s2y)
  tau2y_raw<- var(s$ybar)
  tau2_y   <- max(0, tau2y_raw - mean(s$s2y / s$n))
  list(mu0 = mu0, sigma2 = sigma2, tau2_raw = tau2_raw, tau2 = tau2,
       mu0_y = mu0_y, sigma2_y = sigma2_y, tau2_y_raw = tau2y_raw, tau2_y = tau2_y)
}

priors <- rbindlist(lapply(tissues, function(t) {
  p <- prior_for(t)
  data.table(cancer_type = t,
             n_t = ts[cancer_type == t, n],
             mu0 = p$mu0, tau2_raw = p$tau2_raw, tau2 = p$tau2, sigma2 = p$sigma2,
             mu0_y = p$mu0_y, tau2_y_raw = p$tau2_y_raw, tau2_y = p$tau2_y,
             sigma2_y = p$sigma2_y,
             true_offset_reference_only = ts[cancer_type == t, rbar])
}))

write.table(priors, file.path(OUT_DIR, "shrinkage_priors.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --------------------------------------------------------------- draw engine
draw_rows <- list()   # per (tissue,k) per-draw records

for (t in tissues) {
  sub  <- dat[cancer_type == t]
  n    <- nrow(sub)
  r    <- sub$resid
  y    <- sub$actual
  p    <- prior_for(t)
  b_oracle_full <- mean(r)

  for (k in K_GRID) {
    if (n < k + MIN_EVAL) next
    w   <- p$tau2   / (p$tau2   + p$sigma2   / k)
    w_u <- p$tau2_raw / (p$tau2_raw + p$sigma2 / k)   # sensitivity: uncorrected tau2
    w_y <- p$tau2_y / (p$tau2_y + p$sigma2_y / k)

    m1 <- m2 <- m3 <- m3u <- m4 <- m5 <- m6 <- numeric(B)
    eu <- es <- numeric(B)
    b4e <- numeric(B)

    for (bi in seq_len(B)) {
      idx  <- sample.int(n, k)
      rE   <- r[-idx]; yE <- y[-idx]
      rbar_t <- mean(r[idx]); ybar_t <- mean(y[idx])

      b_un <- rbar_t
      b_sh <- w   * rbar_t + (1 - w)   * p$mu0
      b_su <- w_u * rbar_t + (1 - w_u) * p$mu0
      yhat <- w_y * ybar_t + (1 - w_y) * p$mu0_y

      m1[bi]  <- mean(abs(rE))
      m2[bi]  <- mean(abs(rE - b_un))
      m3[bi]  <- mean(abs(rE - b_sh))
      m3u[bi] <- mean(abs(rE - b_su))
      m4[bi]  <- mean(abs(rE - b_oracle_full))
      m5[bi]  <- mean(abs(yE - ybar_t))
      m6[bi]  <- mean(abs(yE - yhat))
      eu[bi]  <- abs(b_un - b_oracle_full)
      es[bi]  <- abs(b_sh - b_oracle_full)
      b4e[bi] <- mean(rE)                       # b_oracle_eval (recorded)
    }

    draw_rows[[paste(t, k)]] <- data.table(
      cancer_type = t, k = k, n = n, n_eval = n - k,
      w_mean = w, w_uncorrected = w_u, w_y = w_y,
      mu0 = p$mu0, tau2 = p$tau2, tau2_raw = p$tau2_raw, sigma2 = p$sigma2,
      b_oracle_full = b_oracle_full,
      b_oracle_eval_mean = mean(b4e),
      m1 = m1, m2 = m2, m3 = m3, m3_uncorrected_tau2 = m3u,
      m4 = m4, m5 = m5, m6 = m6,
      offset_err_unshrunk = eu, offset_err_shrunk = es)
  }
}

draws <- rbindlist(draw_rows)
stopifnot(nrow(draws) > 0)

# ----------------------------------------------------------- (a) per tissue,k
by_tk <- draws[, .(
  n = n[1], n_eval = n_eval[1],
  w_mean = w_mean[1], w_uncorrected = w_uncorrected[1], w_y = w_y[1],
  mu0 = mu0[1], tau2 = tau2[1], tau2_raw = tau2_raw[1], sigma2 = sigma2[1],
  b_oracle_full = b_oracle_full[1],
  b_oracle_eval_mean = mean(b_oracle_eval_mean),
  m1_model_uncorrected      = mean(m1),
  m2_model_fewshot_unshrunk = mean(m2),
  m3_model_fewshot_shrunk   = mean(m3),
  m3_shrunk_uncorrected_tau2= mean(m3_uncorrected_tau2),
  m4_model_oracle           = mean(m4),
  m5_null_fewshot_flat      = mean(m5),
  m6_null_fewshot_shrunk    = mean(m6),
  offset_err_unshrunk_mean = mean(offset_err_unshrunk),
  offset_err_unshrunk_sd   = sd(offset_err_unshrunk),
  offset_err_shrunk_mean   = mean(offset_err_shrunk),
  offset_err_shrunk_sd     = sd(offset_err_shrunk),
  frac_m3_lt_m1 = mean(m3 < m1),
  frac_m3_lt_m2 = mean(m3 < m2),
  frac_m3_lt_m5 = mean(m3 < m5),
  frac_m3_lt_m6 = mean(m3 < m6)
), by = .(cancer_type, k)][order(k, cancer_type)]

write.table(by_tk, file.path(OUT_DIR, "shrinkage_by_tissue_k.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------------------------------------- (b) macro by k
macro <- by_tk[, {
  mm1 <- mean(m1_model_uncorrected); mm3 <- mean(m3_model_fewshot_shrunk)
  mm2 <- mean(m2_model_fewshot_unshrunk); mm4 <- mean(m4_model_oracle)
  .(n_tissues_evaluated = .N,
    macro_m1_model_uncorrected      = mm1,
    macro_m2_model_fewshot_unshrunk = mm2,
    macro_m3_model_fewshot_shrunk   = mm3,
    macro_m3_shrunk_uncorrected_tau2= mean(m3_shrunk_uncorrected_tau2),
    macro_m4_model_oracle           = mm4,
    macro_m5_null_fewshot_flat      = mean(m5_null_fewshot_flat),
    macro_m6_null_fewshot_shrunk    = mean(m6_null_fewshot_shrunk),
    macro_offset_err_shrunk   = mean(offset_err_shrunk_mean),
    macro_offset_err_unshrunk = mean(offset_err_unshrunk_mean),
    mean_w = mean(w_mean),
    n_tissues_helped              = sum(m3_model_fewshot_shrunk < m1_model_uncorrected),
    n_tissues_helped_unshrunk     = sum(m2_model_fewshot_unshrunk < m1_model_uncorrected),
    n_tissues_beating_flat_null   = sum(m3_model_fewshot_shrunk < m5_null_fewshot_flat),
    n_tissues_beating_shrunk_null = sum(m3_model_fewshot_shrunk < m6_null_fewshot_shrunk),
    fraction_of_oracle_gain_recovered          = (mm1 - mm3) / (mm1 - mm4),
    fraction_of_oracle_gain_recovered_unshrunk = (mm1 - mm2) / (mm1 - mm4))
}, by = k][order(k)]

write.table(macro, file.path(OUT_DIR, "shrinkage_macro_by_k.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --------------------------------------------------- (c) draw distribution
QP <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)
qlab <- paste0("q", sub("0\\.", "", format(QP)))

qsumm <- function(x) as.list(setNames(as.numeric(quantile(x, QP)), qlab))

dist_tissue <- draws[, c(.(level = "tissue", n_draws = .N, delta = "m1_minus_m3"),
                         qsumm(m1 - m3)), by = .(cancer_type, k)]
dist_tissue2 <- draws[, c(.(level = "tissue", n_draws = .N, delta = "m2_minus_m3"),
                          qsumm(m2 - m3)), by = .(cancer_type, k)]
# Equal weight per tissue: every evaluated tissue contributes exactly B draws at
# a given k, so pooling all draws at that k is already an equal-weight pool.
dist_macro <- draws[, c(.(cancer_type = "ALL", level = "macro_pooled",
                          n_draws = .N, delta = "m1_minus_m3"),
                        qsumm(m1 - m3)), by = k]
dist_macro2 <- draws[, c(.(cancer_type = "ALL", level = "macro_pooled",
                           n_draws = .N, delta = "m2_minus_m3"),
                         qsumm(m2 - m3)), by = k]

dist_all <- rbindlist(list(dist_tissue, dist_tissue2, dist_macro, dist_macro2),
                      use.names = TRUE, fill = TRUE)
setcolorder(dist_all, c("level", "cancer_type", "k", "delta", "n_draws", qlab))
write.table(dist_all, file.path(OUT_DIR, "shrinkage_draw_distribution.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------------------------------- (d) offset strata
by_tk[, offset_stratum := fifelse(abs(b_oracle_full) < 2, "small (<2)",
                          fifelse(abs(b_oracle_full) < 5, "moderate (2-5)", ">=5"))]
strat <- by_tk[, .(n_tissues = .N,
                   macro_m1 = mean(m1_model_uncorrected),
                   macro_m2 = mean(m2_model_fewshot_unshrunk),
                   macro_m3 = mean(m3_model_fewshot_shrunk),
                   macro_m4 = mean(m4_model_oracle),
                   macro_m5 = mean(m5_null_fewshot_flat),
                   macro_m6 = mean(m6_null_fewshot_shrunk),
                   n_helped = sum(m3_model_fewshot_shrunk < m1_model_uncorrected),
                   n_beating_m6 = sum(m3_model_fewshot_shrunk < m6_null_fewshot_shrunk)),
               by = .(offset_stratum, k)][order(k, offset_stratum)]
write.table(strat, file.path(OUT_DIR, "shrinkage_by_offset_stratum.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------------------------------------------ stdout
cat("\n=== PHASE A: macro-over-tissue MAE by k (equal weight per tissue) ===\n")
pr <- macro[, .(k, n_tiss = n_tissues_evaluated,
                m1 = round(macro_m1_model_uncorrected, 3),
                m2 = round(macro_m2_model_fewshot_unshrunk, 3),
                m3 = round(macro_m3_model_fewshot_shrunk, 3),
                m4 = round(macro_m4_model_oracle, 3),
                m5 = round(macro_m5_null_fewshot_flat, 3),
                m6 = round(macro_m6_null_fewshot_shrunk, 3),
                w = round(mean_w, 3),
                helped = n_tissues_helped,
                beat_m5 = n_tissues_beating_flat_null,
                beat_m6 = n_tissues_beating_shrunk_null,
                orcl_shr = round(fraction_of_oracle_gain_recovered, 3),
                orcl_uns = round(fraction_of_oracle_gain_recovered_unshrunk, 3))]
print(pr, row.names = FALSE)

cat("\n=== Priors (range across the 30 leave-one-tissue-out fits) ===\n")
print(priors[, .(mu0 = sprintf("%.3f-%.3f", min(mu0), max(mu0)),
                 tau2_raw = sprintf("%.2f-%.2f", min(tau2_raw), max(tau2_raw)),
                 tau2 = sprintf("%.2f-%.2f", min(tau2), max(tau2)),
                 sigma2 = sprintf("%.1f-%.1f", min(sigma2), max(sigma2)),
                 tau2_y = sprintf("%.1f-%.1f", min(tau2_y), max(tau2_y)),
                 sigma2_y = sprintf("%.1f-%.1f", min(sigma2_y), max(sigma2_y)))],
      row.names = FALSE)

cat("\n=== Macro MAE by |true tissue offset| stratum ===\n")
print(strat[, .(k, offset_stratum, n_tissues,
                m1 = round(macro_m1, 2), m2 = round(macro_m2, 2),
                m3 = round(macro_m3, 2), m4 = round(macro_m4, 2),
                m6 = round(macro_m6, 2), n_helped, n_beating_m6)],
      row.names = FALSE)

cat("\n=== Macro pooled per-draw delta quantiles (positive = shrinkage helps) ===\n")
print(dist_all[level == "macro_pooled",
               .(k, delta, q05 = round(q05, 3), q25 = round(q25, 3),
                 q50 = round(q50, 3), q75 = round(q75, 3), q95 = round(q95, 3))],
      row.names = FALSE)

cat("\nNOTE: m4 (oracle) uses the FULL-tissue offset, the standard 'best achievable\n",
    "constant correction' benchmark; b_oracle_eval (eval-rows-only offset) is also\n",
    "recorded in shrinkage_by_tissue_k.tsv as b_oracle_eval_mean.\n", sep = "")
cat("NOTE: shrinkage does NOT enable zero-label inference; every method m2/m3\n",
    "requires k>=1 labelled target-tissue samples.\n", sep = "")
cat("\nWrote outputs to ", OUT_DIR, "\n", sep = "")

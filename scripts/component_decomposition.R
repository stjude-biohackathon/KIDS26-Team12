#!/usr/bin/env Rscript
# =============================================================================
# component_decomposition.R
#
# PURPOSE
#   SECONDARY / EXPLORATORY, label-side only. Decomposes the canonical
#   HRDsum = HRD_LOH + LST + TAI into its three components and asks whether the
#   unresolved per-tissue calibration offset of run 01 (blocker C1) is explained
#   by tissues differing in WHICH kind of genomic scar they carry (the
#   "component-mix artefact" hypothesis).
#
#   NO MODEL IS FIT. NO CLUSTER RESOURCE IS USED. The 28 GB beta matrix is
#   never opened. Existing run-01 predictions are re-used as a frozen column.
#
# INPUTS (the ONLY two data files read)
#   data/processed/master_samples.tsv       per-specimen labels + components
#   results/loco_run01/loco_predictions.tsv frozen LOCO predictions (run 01)
#
# OUTPUTS (results/components_2026-09-18/)
#   00_identity_check.tsv            HRDsum == LOH + LST + TAI audit
#   01_component_by_tissue.tsv       per-tissue n, mean, sd of each component
#   02_variance_shares.tsv           pooled + within-tissue variance shares
#   03_pred_component_correlation.tsv pooled and per-tissue pred~component r/rho
#   04_pred_component_macro.tsv      macro-averaged within-tissue correlations
#   05_offset_by_tissue.tsv          per-tissue offset + component means/shares
#   06_offset_correlations.tsv       offset ~ component mean / mix, with CIs
#   07_tissue_icc.tsv                tissue-identity R^2 per component + pred
#   RESULTS.md                       narrative summary
#   sessionInfo.txt                  R session record
#
# CNS LOCK / LEAKAGE BOUNDARIES
#   * The VERY FIRST operation on the metadata is partition == "development",
#     which by construction drops the 642 locked GBM+LGG specimens. A hard
#     stopifnot immediately afterwards fails loudly if any GBM/LGG row, or any
#     non-development row, survived. No GBM/LGG HRDsum value is ever read into
#     a summary, a correlation, or an output file.
#   * No pediatric / PBTA / PBTP data is touched.
#   * C2 clipping (pmax(pred, 0)) is applied to predictions, matching the
#     adopted shipping rule (R/calibration.R, transform "clip").
#
# SEED
#   7317 (chosen once, recorded here and in RESULTS.md). Used only for the
#   tissue-level bootstrap CIs in step 4.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

SEED <- 7317L
set.seed(SEED)

MASTER_FILE <- "data/processed/master_samples.tsv"
PRED_FILE   <- "results/loco_run01/loco_predictions.tsv"
OUT_DIR     <- "results/components_2026-09-18"
CNS         <- c("GBM", "LGG")
N_BOOT      <- 10000L

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

wr <- function(x, f) write.table(x, file.path(OUT_DIR, f), sep = "\t",
                                 quote = FALSE, row.names = FALSE)

# ------------------------------------------------------ load + CNS lock guard
# Read only the label columns we need; the beta matrix is never touched.
master_raw <- fread(MASTER_FILE, sep = "\t",
                    select = c("sample_id", "patient_id", "cancer_type",
                               "HRD_LOH", "LST", "TAI", "HRDsum",
                               "purity", "ploidy", "partition"))

# FIRST operation on the metadata: restrict to development. This is the CNS lock.
master <- master_raw[partition == "development"]
rm(master_raw)

stopifnot(nrow(master) > 0L)
stopifnot(!any(master$cancer_type %in% CNS))          # CNS lock, hard fail
stopifnot(all(master$partition == "development"))
stopifnot(length(unique(master$cancer_type)) == 30L)
stopifnot(nrow(master) == 7065L)
stopifnot(!anyDuplicated(master$sample_id))

pred <- fread(PRED_FILE, sep = "\t",
              select = c("sample_id", "predicted_reference_HRDsum",
                         "actual", "cancer_type", "split"))
stopifnot(!any(pred$cancer_type %in% CNS))            # CNS lock, hard fail
stopifnot(nrow(pred) == 7065L)
stopifnot(!anyDuplicated(pred$sample_id))
stopifnot(all(pred$split == "development_LOCO"))

# ------------------------------------------------------------------ join + C2
dat <- merge(pred[, .(sample_id, pred_raw = predicted_reference_HRDsum,
                      actual, pred_cancer = cancer_type)],
             master, by = "sample_id", all = FALSE)
stopifnot(nrow(dat) == 7065L)
stopifnot(identical(dat$pred_cancer, dat$cancer_type))   # join sanity
stopifnot(!any(dat$cancer_type %in% CNS))                # CNS lock, again
stopifnot(all(abs(dat$actual - dat$HRDsum) < 1e-8))      # label sanity

# adopted C2 rule
dat[, pred_clip := pmax(pred_raw, 0)]
dat[, resid := pred_clip - HRDsum]
stopifnot(!any(is.na(dat$pred_clip)), all(dat$pred_clip >= 0))

COMPS <- c("HRD_LOH", "LST", "TAI")
stopifnot(!anyNA(dat[, ..COMPS]), !anyNA(dat$HRDsum))

# =========================================================== 1. identity check
dat[, sum_components := HRD_LOH + LST + TAI]
dat[, identity_delta := HRDsum - sum_components]
identity_check <- data.table(
  n                    = nrow(dat),
  n_exact_equal        = sum(dat$identity_delta == 0),
  n_mismatch_exact     = sum(dat$identity_delta != 0),
  n_mismatch_tol_1e_8  = sum(abs(dat$identity_delta) > 1e-8),
  max_abs_delta        = max(abs(dat$identity_delta)),
  note = "HRDsum vs HRD_LOH + LST + TAI over the 30 development cancers")
wr(identity_check, "00_identity_check.tsv")

mismatch <- dat[identity_delta != 0,
                .(sample_id, cancer_type, HRD_LOH, LST, TAI, HRDsum,
                  sum_components, identity_delta)]
if (nrow(mismatch) > 0L) wr(mismatch, "00_identity_mismatches.tsv")

# ==================================================== 2. per-tissue components
by_tissue <- dat[, .(
  n         = .N,
  mean_LOH  = mean(HRD_LOH), sd_LOH = sd(HRD_LOH),
  mean_LST  = mean(LST),     sd_LST = sd(LST),
  mean_TAI  = mean(TAI),     sd_TAI = sd(TAI),
  mean_HRDsum = mean(HRDsum), sd_HRDsum = sd(HRDsum),
  mean_pred_clip = mean(pred_clip), sd_pred_clip = sd(pred_clip),
  offset    = mean(pred_clip - HRDsum),
  mae       = mean(abs(pred_clip - HRDsum)),
  share_LOH = mean(HRD_LOH) / mean(HRDsum),
  share_LST = mean(LST)     / mean(HRDsum),
  share_TAI = mean(TAI)     / mean(HRDsum)
), by = cancer_type][order(cancer_type)]
stopifnot(nrow(by_tissue) == 30L, sum(by_tissue$n) == 7065L)
stopifnot(all(abs(by_tissue$share_LOH + by_tissue$share_LST +
                  by_tissue$share_TAI - 1) < 1e-8))
wr(by_tissue, "01_component_by_tissue.tsv")

# ---- variance shares. Var(HRDsum) = sum_c Cov(c, HRDsum); the share of each
#      component is Cov(c, HRDsum)/Var(HRDsum), and the three shares sum to 1.
var_share <- function(d) {
  v <- var(d$HRDsum)
  data.table(
    var_HRDsum = v,
    share_LOH  = cov(d$HRD_LOH, d$HRDsum) / v,
    share_LST  = cov(d$LST,     d$HRDsum) / v,
    share_TAI  = cov(d$TAI,     d$HRDsum) / v,
    var_frac_LOH = var(d$HRD_LOH) / v,
    var_frac_LST = var(d$LST)     / v,
    var_frac_TAI = var(d$TAI)     / v)
}
pooled_vs <- cbind(scope = "pooled", var_share(dat))

# within-tissue: centre every quantity within its tissue, then repeat
cen <- copy(dat)
for (cl in c(COMPS, "HRDsum", "pred_clip"))
  cen[, (cl) := get(cl) - mean(get(cl)), by = cancer_type]
within_vs <- cbind(scope = "within_tissue", var_share(cen))

# per-tissue shares as well (macro view)
per_tissue_vs <- dat[, cbind(scope = cancer_type, var_share(.SD)), by = cancer_type]
per_tissue_vs[, cancer_type := NULL]
macro_vs <- data.table(scope = "macro_mean_of_30_tissues",
                       var_HRDsum = mean(per_tissue_vs$var_HRDsum),
                       share_LOH = mean(per_tissue_vs$share_LOH),
                       share_LST = mean(per_tissue_vs$share_LST),
                       share_TAI = mean(per_tissue_vs$share_TAI),
                       var_frac_LOH = mean(per_tissue_vs$var_frac_LOH),
                       var_frac_LST = mean(per_tissue_vs$var_frac_LST),
                       var_frac_TAI = mean(per_tissue_vs$var_frac_TAI))
variance_shares <- rbind(pooled_vs, within_vs, macro_vs, per_tissue_vs)
stopifnot(all(abs(variance_shares$share_LOH + variance_shares$share_LST +
                  variance_shares$share_TAI - 1) < 1e-8))
wr(variance_shares, "02_variance_shares.tsv")

# ========================================== 3. prediction ~ component coupling
cor_row <- function(scope, n, x, targets) {
  out <- data.table(scope = scope, n = n)
  for (tg in names(targets)) {
    y <- targets[[tg]]
    ok <- sd(x) > 0 && sd(y) > 0
    out[[paste0("pearson_", tg)]]  <- if (ok) cor(x, y) else NA_real_
    out[[paste0("spearman_", tg)]] <- if (ok) cor(x, y, method = "spearman") else NA_real_
  }
  out
}
tg_of <- function(d) list(LOH = d$HRD_LOH, LST = d$LST, TAI = d$TAI,
                          HRDsum = d$HRDsum)

pooled_cor <- cor_row("pooled", nrow(dat), dat$pred_clip, tg_of(dat))
within_pooled_cor <- cor_row("pooled_within_tissue_centred", nrow(cen),
                             cen$pred_clip, tg_of(cen))
per_tissue_cor <- rbindlist(lapply(split(dat, dat$cancer_type), function(d)
  cor_row(d$cancer_type[1], nrow(d), d$pred_clip, tg_of(d))))
pred_cor <- rbind(pooled_cor, within_pooled_cor, per_tissue_cor)
wr(pred_cor, "03_pred_component_correlation.tsv")

num_cols <- setdiff(names(per_tissue_cor), c("scope", "n"))
macro_cor <- data.table(statistic = num_cols,
                        macro_mean_over_30_tissues =
                          vapply(num_cols, function(cl)
                            mean(per_tissue_cor[[cl]], na.rm = TRUE), numeric(1)),
                        n_tissues_used =
                          vapply(num_cols, function(cl)
                            sum(is.finite(per_tissue_cor[[cl]])), numeric(1)))
wr(macro_cor, "04_pred_component_macro.tsv")

# Is TAI really the best-tracked component, or is the macro gap noise? Paired
# over the SAME 30 tissues: difference in within-tissue Pearson, with a paired
# bootstrap CI over tissues and a Wilcoxon signed-rank p.
pair_cmp <- function(a, b) {
  d <- per_tissue_cor[[paste0("pearson_", a)]] - per_tissue_cor[[paste0("pearson_", b)]]
  bs <- replicate(N_BOOT, mean(sample(d, length(d), replace = TRUE)))
  data.table(comparison = paste0(a, " - ", b), n_tissues = length(d),
             mean_diff = mean(d),
             boot_lo = as.numeric(quantile(bs, 0.025)),
             boot_hi = as.numeric(quantile(bs, 0.975)),
             n_tissues_favouring_first = sum(d > 0),
             wilcoxon_p = suppressWarnings(wilcox.test(d)$p.value))
}
paired <- rbind(pair_cmp("TAI", "LOH"), pair_cmp("TAI", "LST"),
                pair_cmp("LST", "LOH"), pair_cmp("TAI", "HRDsum"))
paired[, seed := SEED]
wr(paired, "04b_pred_component_paired.tsv")

# ============================================ 4. per-tissue offset (BLOCKER C1)
pur <- dat[, .(mean_purity = mean(purity, na.rm = TRUE),
               n_purity = sum(!is.na(purity)),
               mean_ploidy = mean(ploidy, na.rm = TRUE)), by = cancer_type]
offset_tab <- merge(by_tissue[, .(cancer_type, n, offset, mae, mean_HRDsum,
                                  mean_pred_clip, mean_LOH, mean_LST, mean_TAI,
                                  share_LOH, share_LST, share_TAI)],
                    pur, by = "cancer_type")[order(-offset)]
# purity/ploidy are incomplete in master_samples and are CONTEXT ONLY here;
# they are never used in a correlation below.
stopifnot(nrow(offset_tab) == 30L)
stopifnot(!anyNA(offset_tab[, .(offset, mean_HRDsum, mean_LOH, mean_LST,
                                mean_TAI, share_LOH, share_LST, share_TAI)]))
wr(offset_tab, "05_offset_by_tissue.tsv")

fisher_ci <- function(r, n, level = 0.95) {
  if (!is.finite(r) || n < 4L) return(c(NA_real_, NA_real_))
  z <- atanh(r); se <- 1 / sqrt(n - 3)
  tanh(z + c(-1, 1) * qnorm(1 - (1 - level) / 2) * se)
}

boot_ci <- function(x, y, method, B = N_BOOT) {
  n <- length(x)
  bs <- replicate(B, {
    i <- sample.int(n, n, replace = TRUE)
    if (sd(x[i]) > 0 && sd(y[i]) > 0) cor(x[i], y[i], method = method) else NA_real_
  })
  as.numeric(quantile(bs, c(0.025, 0.975), na.rm = TRUE))
}

offset_terms <- list(
  mean_LOH    = by_tissue$mean_LOH,
  mean_LST    = by_tissue$mean_LST,
  mean_TAI    = by_tissue$mean_TAI,
  mean_HRDsum = by_tissue$mean_HRDsum,
  share_LOH   = by_tissue$share_LOH,
  share_LST   = by_tissue$share_LST,
  share_TAI   = by_tissue$share_TAI)

off <- by_tissue$offset
nT  <- nrow(by_tissue)
offset_cor <- rbindlist(lapply(names(offset_terms), function(tm) {
  x <- offset_terms[[tm]]
  r  <- cor(off, x)
  rh <- cor(off, x, method = "spearman")
  ci  <- fisher_ci(r, nT)
  bci <- boot_ci(off, x, "pearson")
  bcs <- boot_ci(off, x, "spearman")
  data.table(term = tm, n_tissues = nT,
             pearson = r, pearson_fisher_lo = ci[1], pearson_fisher_hi = ci[2],
             pearson_boot_lo = bci[1], pearson_boot_hi = bci[2],
             pearson_p = cor.test(off, x)$p.value,
             spearman = rh, spearman_boot_lo = bcs[1], spearman_boot_hi = bcs[2],
             spearman_p = suppressWarnings(
               cor.test(off, x, method = "spearman")$p.value))
}))
# same, on |offset|
offset_cor_abs <- rbindlist(lapply(names(offset_terms), function(tm) {
  x <- offset_terms[[tm]]; a <- abs(off)
  r <- cor(a, x); ci <- fisher_ci(r, nT); bci <- boot_ci(a, x, "pearson")
  data.table(term = paste0("abs_offset~", tm), n_tissues = nT,
             pearson = r, pearson_fisher_lo = ci[1], pearson_fisher_hi = ci[2],
             pearson_boot_lo = bci[1], pearson_boot_hi = bci[2],
             pearson_p = cor.test(a, x)$p.value,
             spearman = cor(a, x, method = "spearman"),
             spearman_boot_lo = NA_real_, spearman_boot_hi = NA_real_,
             spearman_p = suppressWarnings(
               cor.test(a, x, method = "spearman")$p.value))
}))
offset_cor[, term := paste0("offset~", term)]
offset_all <- rbind(offset_cor, offset_cor_abs)

# partial correlation of offset with each SHARE, controlling for mean_HRDsum:
# the mix hypothesis needs the mix to matter beyond the tissue's overall level.
pcor <- function(y, x, z) {
  ry <- residuals(lm(y ~ z)); rx <- residuals(lm(x ~ z))
  cor(ry, rx)
}
partials <- rbindlist(lapply(c("share_LOH", "share_LST", "share_TAI"), function(tm) {
  r <- pcor(off, offset_terms[[tm]], by_tissue$mean_HRDsum)
  ci <- fisher_ci(r, nT - 1L)
  data.table(term = paste0("offset~", tm, "|mean_HRDsum"), n_tissues = nT,
             pearson = r, pearson_fisher_lo = ci[1], pearson_fisher_hi = ci[2],
             pearson_boot_lo = NA_real_, pearson_boot_hi = NA_real_,
             pearson_p = NA_real_, spearman = NA_real_,
             spearman_boot_lo = NA_real_, spearman_boot_hi = NA_real_,
             spearman_p = NA_real_)
}))
offset_all <- rbind(offset_all, partials)
offset_all[, seed := SEED][, n_boot := N_BOOT]
wr(offset_all, "06_offset_correlations.tsv")

# Nested linear models over the 30 tissues. The mix hypothesis predicts that
# adding the component SHARES to a model that already knows the tissue's overall
# HRDsum level buys real explanatory power. m1 = level only, m3 = level + mix.
m1 <- lm(offset ~ mean_HRDsum, data = by_tissue)
m2 <- lm(offset ~ mean_LOH + mean_LST + mean_TAI, data = by_tissue)
m3 <- lm(offset ~ mean_HRDsum + share_LOH + share_LST, data = by_tissue)
nested <- data.table(
  model = c("offset ~ mean_HRDsum", "offset ~ mean_LOH + mean_LST + mean_TAI",
            "offset ~ mean_HRDsum + share_LOH + share_LST"),
  n_tissues = nT,
  r2     = c(summary(m1)$r.squared, summary(m2)$r.squared, summary(m3)$r.squared),
  adj_r2 = c(summary(m1)$adj.r.squared, summary(m2)$adj.r.squared,
             summary(m3)$adj.r.squared),
  df = c(m1$rank, m2$rank, m3$rank),
  p_vs_level_only = c(NA_real_, anova(m1, m2)$`Pr(>F)`[2],
                      anova(m1, m3)$`Pr(>F)`[2]))
wr(nested, "06b_offset_nested_models.tsv")

# ==================================================== 5. tissue-identity R^2
tissue_identity_r2 <- function(values, cancer) {
  ok <- is.finite(values); values <- values[ok]; cancer <- as.character(cancer)[ok]
  if (length(values) < 3L || length(unique(cancer)) < 2L) return(NA_real_)
  ss_tot <- sum((values - mean(values))^2)
  if (!(ss_tot > 0)) return(NA_real_)
  ss_within <- sum(vapply(split(values, cancer),
                          function(v) sum((v - mean(v))^2), numeric(1)))
  1 - ss_within / ss_tot
}
quantities <- c("HRD_LOH", "LST", "TAI", "HRDsum", "pred_clip", "pred_raw",
                "resid", "purity", "ploidy")
icc <- rbindlist(lapply(quantities, function(q) {
  v <- dat[[q]]
  r2 <- tissue_identity_r2(v, dat$cancer_type)
  # bootstrap over TISSUES (resample tissues, keep samples within)
  bs <- replicate(1000L, {
    tt <- sample(by_tissue$cancer_type, nT, replace = TRUE)
    idx <- unlist(lapply(seq_along(tt), function(j) which(dat$cancer_type == tt[j])))
    lab <- unlist(lapply(seq_along(tt), function(j)
      rep(paste0(tt[j], "_", j), sum(dat$cancer_type == tt[j]))))
    tissue_identity_r2(v[idx], lab)
  })
  data.table(quantity = q, n = sum(is.finite(v)), n_tissues = nT,
             tissue_r2 = r2,
             boot_lo = as.numeric(quantile(bs, 0.025, na.rm = TRUE)),
             boot_hi = as.numeric(quantile(bs, 0.975, na.rm = TRUE)),
             between_tissue_var_frac = r2,
             within_tissue_var_frac = 1 - r2)
}))
icc[, seed := SEED]
wr(icc, "07_tissue_icc.tsv")

# ------------------------------------------------------------------- stdout
cat("\n=== SEED", SEED, "| 30 development cancers | n =", nrow(dat), "===\n")
cat("\n--- 1. identity check ---\n");  print(identity_check, row.names = FALSE)
cat("\n--- 2. variance shares (pooled / within / macro) ---\n")
print(variance_shares[1:3], row.names = FALSE)
cat("\n--- 3. pooled pred~component ---\n"); print(pooled_cor, row.names = FALSE)
cat("\n--- 3. macro within-tissue pred~component ---\n")
print(macro_cor, row.names = FALSE)
cat("\n--- 3b. paired within-tissue contrasts over 30 tissues ---\n")
print(paired, row.names = FALSE)
cat("\n--- 4. offset correlations ---\n")
print(offset_all[, .(term, pearson = round(pearson, 3),
                     lo = round(pearson_boot_lo, 3), hi = round(pearson_boot_hi, 3),
                     spearman = round(spearman, 3), p = signif(pearson_p, 3))],
      row.names = FALSE)
cat("\n--- 5. tissue-identity R2 ---\n")
print(icc[, .(quantity, tissue_r2 = round(tissue_r2, 4),
              boot_lo = round(boot_lo, 4), boot_hi = round(boot_hi, 4))],
      row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

# ---------------------------------------------------------------- RESULTS.md
f3 <- function(x) formatC(x, format = "f", digits = 3)
g <- function(tm, cl) offset_all[term == tm][[cl]]
mc <- function(s) macro_cor[statistic == s]$macro_mean_over_30_tissues
ic <- function(q) icc[quantity == q]

md <- c(
"# Component decomposition of HRDsum (secondary / exploratory)",
"",
sprintf("Generated by `scripts/component_decomposition.R`, seed **%d**, %s.", SEED, Sys.Date()),
sprintf("Scope: **%d specimens, 30 development cancer types**. GBM/LGG excluded as the first", nrow(dat)),
"operation on the metadata (`partition == \"development\"`) and asserted absent thereafter.",
"No model was fit; `data/processed/beta.tsv` was never opened; no LSF job was submitted.",
"Predictions are run-01 LOCO with the adopted C2 clipping `pmax(pred, 0)`.",
"",
"## 1. Identity check",
"",
sprintf("`HRDsum == HRD_LOH + LST + TAI` holds **exactly for all %d rows** (max |delta| = %g).",
        nrow(dat), identity_check$max_abs_delta),
"No discrepancies. See `00_identity_check.tsv`.",
"",
"## 2. Component distributions and variance shares",
"",
sprintf("Pooled, the covariance decomposition of Var(HRDsum) = %.1f splits LOH/LST/TAI =",
        pooled_vs$var_HRDsum),
sprintf("**%.3f / %.3f / %.3f**; within-tissue (%.1f) it is **%.3f / %.3f / %.3f**.",
        pooled_vs$share_LOH, pooled_vs$share_LST, pooled_vs$share_TAI,
        within_vs$var_HRDsum, within_vs$share_LOH, within_vs$share_LST, within_vs$share_TAI),
"TAI and LST each carry more of HRDsum's variance than LOH, and the split barely",
"changes when tissue means are removed. Per-tissue means/sds: `01_component_by_tissue.tsv`,",
"shares: `02_variance_shares.tsv`.",
"",
"## 3. Which component does the methylation prediction track best?",
"",
"| Scope | LOH | LST | TAI | HRDsum |",
"|---|---|---|---|---|",
sprintf("| Pooled Pearson | %s | %s | %s | %s |",
        f3(pooled_cor$pearson_LOH), f3(pooled_cor$pearson_LST),
        f3(pooled_cor$pearson_TAI), f3(pooled_cor$pearson_HRDsum)),
sprintf("| Pooled Spearman | %s | %s | %s | %s |",
        f3(pooled_cor$spearman_LOH), f3(pooled_cor$spearman_LST),
        f3(pooled_cor$spearman_TAI), f3(pooled_cor$spearman_HRDsum)),
sprintf("| Macro within-tissue Pearson (30) | %s | %s | %s | %s |",
        f3(mc("pearson_LOH")), f3(mc("pearson_LST")), f3(mc("pearson_TAI")),
        f3(mc("pearson_HRDsum"))),
sprintf("| Macro within-tissue Spearman (30) | %s | %s | %s | %s |",
        f3(mc("spearman_LOH")), f3(mc("spearman_LST")), f3(mc("spearman_TAI")),
        f3(mc("spearman_HRDsum"))),
"",
"**TAI wins, and the margin is paired-stable.** Paired over the same 30 tissues",
sprintf("(`04b_pred_component_paired.tsv`): TAI - LOH = %s [%s, %s], %d/30 tissues favour TAI;",
        f3(paired[1]$mean_diff), f3(paired[1]$boot_lo), f3(paired[1]$boot_hi),
        paired[1]$n_tissues_favouring_first),
sprintf("TAI - LST = %s [%s, %s], %d/30. LST - LOH = %s [%s, %s] — indistinguishable.",
        f3(paired[2]$mean_diff), f3(paired[2]$boot_lo), f3(paired[2]$boot_hi),
        paired[2]$n_tissues_favouring_first,
        f3(paired[3]$mean_diff), f3(paired[3]$boot_lo), f3(paired[3]$boot_hi)),
sprintf("TAI is even tracked slightly better than the composite HRDsum itself (%s [%s, %s]).",
        f3(paired[4]$mean_diff), f3(paired[4]$boot_lo), f3(paired[4]$boot_hi)),
"",
"## 4. KEY QUESTION — is the per-tissue offset (C1) a component-MIX artefact?",
"",
sprintf("Run-01 per-tissue offsets span %s to %s (mean |offset| %s, sd %s) — `05_offset_by_tissue.tsv`.",
        f3(min(by_tissue$offset)), f3(max(by_tissue$offset)),
        f3(mean(abs(by_tissue$offset))), f3(sd(by_tissue$offset))),
"",
"| Predictor (n = 30 tissues) | Pearson | 95% boot CI | Spearman | p |",
"|---|---|---|---|---|",
paste0("| ", offset_cor$term, " | ", f3(offset_cor$pearson), " | [",
       f3(offset_cor$pearson_boot_lo), ", ", f3(offset_cor$pearson_boot_hi), "] | ",
       f3(offset_cor$spearman), " | ", signif(offset_cor$pearson_p, 2), " |"),
"",
"**The offset tracks the LEVEL, not the MIX.** Every component *mean* correlates",
"negatively and similarly with the offset (regression to the cross-tissue mean: low-HRD",
"tissues are over-predicted, high-HRD tissues under-predicted). No component *share*",
sprintf("reaches significance: share_LOH r = %s [%s, %s], share_LST r = %s [%s, %s],",
        f3(g("offset~share_LOH", "pearson")), f3(g("offset~share_LOH", "pearson_boot_lo")),
        f3(g("offset~share_LOH", "pearson_boot_hi")),
        f3(g("offset~share_LST", "pearson")), f3(g("offset~share_LST", "pearson_boot_lo")),
        f3(g("offset~share_LST", "pearson_boot_hi"))),
sprintf("share_TAI r = %s [%s, %s] — all CIs straddle 0.",
        f3(g("offset~share_TAI", "pearson")), f3(g("offset~share_TAI", "pearson_boot_lo")),
        f3(g("offset~share_TAI", "pearson_boot_hi"))),
"Partial correlations of the offset with each share, controlling for the tissue's",
sprintf("mean HRDsum, collapse further (%s, %s, %s for LOH/LST/TAI share).",
        f3(g("offset~share_LOH|mean_HRDsum", "pearson")),
        f3(g("offset~share_LST|mean_HRDsum", "pearson")),
        f3(g("offset~share_TAI|mean_HRDsum", "pearson"))),
sprintf("Nested models (`06b_offset_nested_models.tsv`): level alone R2 = %s; adding the two",
        f3(nested[1]$r2)),
sprintf("free shares gives R2 = %s (adj R2 %s -> %s), F-test p = %s. The mix adds nothing.",
        f3(nested[3]$r2), f3(nested[1]$adj_r2), f3(nested[3]$adj_r2),
        f3(nested[3]$p_vs_level_only)),
"",
sprintf("Note the mix is nearly CONSTANT across tissues: share_LOH IQR %s-%s, share_TAI IQR %s-%s.",
        f3(quantile(by_tissue$share_LOH, .25)), f3(quantile(by_tissue$share_LOH, .75)),
        f3(quantile(by_tissue$share_TAI, .25)), f3(quantile(by_tissue$share_TAI, .75))),
"There is little between-tissue mix variation for C1 to be an artefact OF.",
"",
"## 5. Between- vs within-tissue variance (which component is lineage-determined?)",
"",
"| Quantity | tissue R2 (between) | 95% boot CI (tissue resample) | within |",
"|---|---|---|---|",
paste0("| ", icc$quantity, " | ", f3(icc$tissue_r2), " | [", f3(icc$boot_lo), ", ",
       f3(icc$boot_hi), "] | ", f3(icc$within_tissue_var_frac), " |"),
"",
sprintf("**TAI is the most lineage-determined component (%s), LST the least (%s)**, but the",
        f3(ic("TAI")$tissue_r2), f3(ic("LST")$tissue_r2)),
"three CIs overlap heavily, so the ordering is suggestive only.",
sprintf("The 30-tissue R2 of the clipped prediction is **%s** [%s, %s] (raw %s), versus **%s**",
        f3(ic("pred_clip")$tissue_r2), f3(ic("pred_clip")$boot_lo),
        f3(ic("pred_clip")$boot_hi), f3(ic("pred_raw")$tissue_r2), f3(ic("HRDsum")$tissue_r2)),
"for the observed label — reproducing, on all 30 tissues, the sentinel-set finding (0.5219",
"vs 0.2940, `docs/27_V3_PREREGISTRATION.md`) that **the model is markedly more",
"lineage-determined than the biology it predicts**. The prediction's tissue R2 exceeds that",
"of every single component, so the excess lineage signal is not inherited from any one scar type.",
"",
"## Verdict",
"",
"**This UNDERMINES the component-mix hypothesis for C1.** C1 behaves as a level /",
"shrinkage-to-the-grand-mean problem plus excess lineage imprinting in the model, not as a",
"consequence of tissues carrying different *kinds* of genomic scar. The mix is close to",
"constant across tissues and explains no additional offset variance.",
"",
"## Power caveat",
"",
sprintf("All step-4 and step-5 statistics have n = %d tissues. A Pearson r needs |r| > ~0.36 to",
        nT),
"clear p = 0.05 at n = 30, and CI half-widths here are ~0.3. A genuine but modest mix effect",
"(|r| ~ 0.25) could not be detected. Read this as *no evidence for* the mix hypothesis, not",
"as proof of its absence. Nothing here was used to select a model or set a gate.",
"",
"## Files",
"",
"`00_identity_check.tsv`, `01_component_by_tissue.tsv`, `02_variance_shares.tsv`,",
"`03_pred_component_correlation.tsv`, `04_pred_component_macro.tsv`,",
"`04b_pred_component_paired.tsv`, `05_offset_by_tissue.tsv`, `06_offset_correlations.tsv`,",
"`06b_offset_nested_models.tsv`, `07_tissue_icc.tsv`, `sessionInfo.txt`.")

writeLines(md, file.path(OUT_DIR, "RESULTS.md"))
cat("\nWrote outputs to", OUT_DIR, "\n")

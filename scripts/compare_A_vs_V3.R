# =============================================================================
# scripts/compare_A_vs_V3.R - Full 30-fold head-to-head, Candidate A vs V3-abs
# =============================================================================
#
# USAGE
#   Rscript scripts/compare_A_vs_V3.R
#
# INPUTS
#   results/loco_run01/loco_predictions.tsv       Candidate A (pooled variance)
#   results/v3_loco_full/folds/predictions_*.tsv  Candidate B (within-tissue)
# OUTPUTS
#   results/tables/tableA2_A_vs_V3_full30.tsv
#   results/tables/tableA2_A_vs_V3_full30.md
#   results/tables/A_vs_V3_per_tissue_full30.tsv
#
# ESTIMATOR NOTE
#   Per-tissue statistics are computed under the adopted C2 clipping and then
#   macro-averaged, identically for both candidates. This macro estimator is
#   NOT the pooled-within estimator behind the 0.612 headline in docs/24. The
#   two are different estimators of different things and must not be quoted
#   side by side.
#
# ON THE HEADLINE QUANTITY
#   delta_R2_excess = (R2_pred_A - R2_truth) - (R2_pred_B - R2_truth)
#   Both candidates predict the SAME labels on the SAME cohort, so R2_truth is
#   identical and cancels exactly. The quantity therefore reduces to
#   R2_pred_A - R2_pred_B. The excess framing is still the right way to read
#   it - it says how much lineage structure the model carries BEYOND what the
#   target legitimately carries - but it adds no information beyond the raw
#   difference in tissue R2 of predictions, and it is reported here as one
#   number with a bootstrap interval rather than two.
#
# UNCERTAINTY
#   Tissue R2 is a property of a single pooled fit, so it has no per-tissue
#   replicate to average over. It is bootstrapped by resampling SAMPLES WITHIN
#   TISSUE (stratified), which holds the 30-tissue design fixed and propagates
#   within-tissue sampling noise. Both models are recomputed on the SAME
#   resample, so the difference is paired.
#   Per-tissue correlations are compared with a paired bootstrap over the 30
#   tissues and a paired Wilcoxon signed-rank test.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

SEED <- 5813
set.seed(SEED)
B_REP <- 10000L

dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

# --- Load -------------------------------------------------------------------
A <- fread("results/loco_run01/loco_predictions.tsv")
pf <- Sys.glob("results/v3_loco_full/folds/predictions_*.tsv")
stopifnot(length(pf) == 30L)
B <- rbindlist(lapply(pf, fread), fill = TRUE)

for (d in list(A, B)) {
  if (any(d$cancer_type %in% c("GBM", "LGG")))
    stop("CNS rows present - refusing to compare on locked data.")
}
A[, pred := pmax(predicted_reference_HRDsum, 0)]
B[, pred := pmax(predicted_reference_HRDsum, 0)]

setkey(A, sample_id); setkey(B, sample_id)
stopifnot(setequal(A$sample_id, B$sample_id))

# Align on a single frame so every statistic is paired sample-for-sample.
D <- merge(A[, .(sample_id, cancer_type, actual, pred_A = pred)],
           B[, .(sample_id, pred_B = pred)], by = "sample_id")
stopifnot(nrow(D) == 7065L)

# --- Per-tissue -------------------------------------------------------------
per <- D[, .(n = .N,
             r_A   = cor(actual, pred_A), r_B   = cor(actual, pred_B),
             rho_A = cor(actual, pred_A, method = "spearman"),
             rho_B = cor(actual, pred_B, method = "spearman"),
             MAE_A = mean(abs(actual - pred_A)), MAE_B = mean(abs(actual - pred_B)),
             bias_A = mean(pred_A - actual),     bias_B = mean(pred_B - actual),
             # Oracle tissue-mean null: the held-out tissue's own mean label.
             # Not available at deployment; it is the honest bar for "did the
             # model learn anything beyond recognising the tissue".
             MAE_null = mean(abs(actual - mean(actual)))),
         by = cancer_type][order(cancer_type)]
per[, `:=`(d_r = r_B - r_A, d_rho = rho_B - rho_A,
           d_MAE = MAE_B - MAE_A, d_absbias = abs(bias_B) - abs(bias_A))]
fwrite(per, "results/tables/A_vs_V3_per_tissue_full30.tsv", sep = "\t")

r2t <- function(y, g) summary(lm(y ~ factor(g)))$r.squared

R2_truth <- r2t(D$actual, D$cancer_type)
R2_A <- r2t(D$pred_A, D$cancer_type)
R2_B <- r2t(D$pred_B, D$cancer_type)
exc_A <- R2_A - R2_truth
exc_B <- R2_B - R2_truth
d_exc <- exc_A - exc_B              # == R2_A - R2_B, see header
frac_removed <- d_exc / exc_A

# --- Bootstrap: tissue R2, stratified within tissue, paired -----------------
idx_by_t <- split(seq_len(nrow(D)), D$cancer_type)
bs <- replicate(B_REP, {
  s <- unlist(lapply(idx_by_t, function(i) sample(i, length(i), TRUE)), use.names = FALSE)
  g <- D$cancer_type[s]
  c(r2t(D$pred_A[s], g), r2t(D$pred_B[s], g), r2t(D$actual[s], g))
})
ci_R2A   <- quantile(bs[1, ], c(.025, .975), names = FALSE)
ci_R2B   <- quantile(bs[2, ], c(.025, .975), names = FALSE)
ci_dexc  <- quantile(bs[1, ] - bs[2, ], c(.025, .975), names = FALSE)
ci_excA  <- quantile(bs[1, ] - bs[3, ], c(.025, .975), names = FALSE)
ci_excB  <- quantile(bs[2, ] - bs[3, ], c(.025, .975), names = FALSE)
ci_frac  <- quantile((bs[1, ] - bs[2, ]) / (bs[1, ] - bs[3, ]), c(.025, .975), names = FALSE)
p_dexc   <- 2 * min(mean(bs[1, ] - bs[2, ] <= 0), mean(bs[1, ] - bs[2, ] >= 0))

# --- Bootstrap over tissues for the paired per-tissue deltas ----------------
bmean <- function(x) quantile(replicate(B_REP, mean(sample(x, length(x), TRUE))),
                              c(.025, .975), names = FALSE)
ci_dr   <- bmean(per$d_r)
ci_drho <- bmean(per$d_rho)
p_r   <- wilcox.test(per$r_B,   per$r_A,   paired = TRUE)$p.value
p_rho <- wilcox.test(per$rho_B, per$rho_A, paired = TRUE)$p.value

# --- Pooled MAE and skill ---------------------------------------------------
mae_A <- mean(abs(D$actual - D$pred_A))
mae_B <- mean(abs(D$actual - D$pred_B))
mae_null_pooled <- D[, mean(abs(actual - mean(actual))), by = cancer_type][
  , weighted.mean(V1, per$n[match(cancer_type, per$cancer_type)])]
# Macro skill: mean over tissues of (oracle null MAE - model MAE).
skill_A <- mean(per$MAE_null - per$MAE_A)
skill_B <- mean(per$MAE_null - per$MAE_B)

fm <- function(x, d = 4) formatC(x, format = "f", digits = d)
ci <- function(l, u, d = 4) sprintf("[%s, %s]", fm(l, d), fm(u, d))

tab <- data.table(
  metric = c("Macro Pearson (within-tissue)", "Macro Spearman (within-tissue)",
             "Pooled MAE", "Macro skill vs oracle tissue-mean null",
             "Mean absolute tissue bias", "Tissue R2 of predictions",
             "Tissue R2 of truth", "Excess tissue R2 (pred - truth)",
             "Fraction of excess removed", "Tissues with positive r",
             "Median per-tissue delta Pearson", "Median per-tissue delta Spearman",
             "Paired difference, Pearson", "Paired difference, Spearman"),
  model_A = c(fm(mean(per$r_A)), fm(mean(per$rho_A)), fm(mae_A), fm(skill_A),
              fm(mean(abs(per$bias_A))), fm(R2_A), fm(R2_truth), fm(exc_A),
              "-", sprintf("%d of 30", sum(per$r_A > 0)), "-", "-", "-", "-"),
  v3_abs = c(fm(mean(per$r_B)), fm(mean(per$rho_B)), fm(mae_B), fm(skill_B),
             fm(mean(abs(per$bias_B))), fm(R2_B), fm(R2_truth), fm(exc_B),
             sprintf("%.1f%%", 100 * frac_removed),
             sprintf("%d of 30", sum(per$r_B > 0)),
             fm(median(per$d_r)), fm(median(per$d_rho)),
             sprintf("mean %s", fm(mean(per$d_r))),
             sprintf("mean %s", fm(mean(per$d_rho)))),
  difference = c(fm(mean(per$r_B) - mean(per$r_A)), fm(mean(per$rho_B) - mean(per$rho_A)),
                 fm(mae_B - mae_A), fm(skill_B - skill_A),
                 fm(mean(abs(per$bias_B)) - mean(abs(per$bias_A))),
                 fm(R2_B - R2_A), "0 (identical by construction)",
                 fm(exc_B - exc_A), sprintf("%s removed", fm(d_exc)),
                 sprintf("%+d", sum(per$r_B > 0) - sum(per$r_A > 0)),
                 sprintf("%d of 30 improved", sum(per$d_r > 0)),
                 sprintf("%d of 30 improved", sum(per$d_rho > 0)),
                 sprintf("p = %.3f (Wilcoxon)", p_r),
                 sprintf("p = %.3f (Wilcoxon)", p_rho)),
  ci_95 = c("-", "-", "-", "-", "-",
            sprintf("A %s | B %s", ci(ci_R2A[1], ci_R2A[2]), ci(ci_R2B[1], ci_R2B[2])),
            "-",
            sprintf("A %s | B %s", ci(ci_excA[1], ci_excA[2]), ci(ci_excB[1], ci_excB[2])),
            ci(100 * ci_frac[1], 100 * ci_frac[2], 1),
            "-", ci(ci_dr[1], ci_dr[2]), ci(ci_drho[1], ci_drho[2]),
            ci(ci_dr[1], ci_dr[2]), ci(ci_drho[1], ci_drho[2]))
)

fwrite(tab, "results/tables/tableA2_A_vs_V3_full30.tsv", sep = "\t")

md <- c("# Table A2 - Candidate A vs V3-abs, full 30-fold LOCO",
        "",
        sprintf("Development LOCO, 7,065 samples, 30 non-CNS cancer types. Generated %s, seed %d, %s bootstrap replicates.",
                format(Sys.time(), "%Y-%m-%d %H:%M %Z"), SEED, format(B_REP, big.mark = ",")),
        "Per-tissue statistics under adopted C2 clipping, macro-averaged; same estimator for both candidates.",
        "",
        paste("|", paste(names(tab), collapse = " | "), "|"),
        paste("|", paste(rep("---", ncol(tab)), collapse = " | "), "|"),
        apply(tab, 1, function(r) paste("|", paste(r, collapse = " | "), "|")),
        "",
        "## Headline quantity",
        "",
        sprintf("delta_R2_excess = %.4f, 95%% CI %s, bootstrap p = %.4f",
                d_exc, ci(ci_dexc[1], ci_dexc[2]), p_dexc),
        "",
        sprintf("R2_truth (%.4f) is identical for both candidates and cancels exactly, so this equals R2_pred_A - R2_pred_B.",
                R2_truth),
        sprintf("V3-abs removes %.1f%% of the excess lineage structure that Candidate A carries beyond the target's own.",
                100 * frac_removed),
        "",
        "## Paired ranking change",
        "",
        sprintf("Pearson:  mean %+.4f, 95%% CI %s, Wilcoxon p = %.3f, %d of 30 tissues improved",
                mean(per$d_r), ci(ci_dr[1], ci_dr[2]), p_r, sum(per$d_r > 0)),
        sprintf("Spearman: mean %+.4f, 95%% CI %s, Wilcoxon p = %.3f, %d of 30 tissues improved",
                mean(per$d_rho), ci(ci_drho[1], ci_drho[2]), p_rho, sum(per$d_rho > 0)))
writeLines(md, "results/tables/tableA2_A_vs_V3_full30.md")

cat("\n=== CANDIDATE A vs V3-abs, FULL 30-FOLD ===\n\n")
print(tab[, .(metric, model_A, v3_abs, difference)])
cat("\n--- Headline ---\n")
cat(sprintf("  delta_R2_excess = %.4f   95%% CI %s   bootstrap p = %.4f\n",
            d_exc, ci(ci_dexc[1], ci_dexc[2]), p_dexc))
cat(sprintf("  R2_truth = %.4f is identical for both and cancels; this IS R2_pred_A - R2_pred_B.\n", R2_truth))
cat(sprintf("  Excess: A %+.4f -> B %+.4f, i.e. %.1f%% of A's excess removed (95%% CI %.1f%%-%.1f%%)\n",
            exc_A, exc_B, 100 * frac_removed, 100 * ci_frac[1], 100 * ci_frac[2]))
cat("\n--- Paired ranking change across the 30 tissues ---\n")
cat(sprintf("  Pearson  mean %+.4f  95%% CI %s  p = %.3f  (%d of 30 improved)\n",
            mean(per$d_r), ci(ci_dr[1], ci_dr[2]), p_r, sum(per$d_r > 0)))
cat(sprintf("  Spearman mean %+.4f  95%% CI %s  p = %.3f  (%d of 30 improved)\n",
            mean(per$d_rho), ci(ci_drho[1], ci_drho[2]), p_rho, sum(per$d_rho > 0)))
cat("\nWritten:\n  results/tables/tableA2_A_vs_V3_full30.{tsv,md}",
    "\n  results/tables/A_vs_V3_per_tissue_full30.tsv\n")

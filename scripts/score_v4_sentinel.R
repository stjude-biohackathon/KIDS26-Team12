# =============================================================================
# scripts/score_v4_sentinel.R - Mechanical scorer for the V4 sentinel gates
# =============================================================================
#
# USAGE
#   Rscript scripts/score_v4_sentinel.R
#
# WHAT THIS IS
#   V4 = supervised within-tissue HRD meta-association feature selection with a
#   lineage penalty, the ONE alternative experiment pre-registered in
#   docs/27_V3_PREREGISTRATION.md section 6. That section closes with "Same four
#   sentinel tissues, same gates G1-G7", so this script applies the SAME
#   thresholds already frozen in docs/27 section 5 (commit d72695c, 10:27 CDT,
#   written while no sentinel output existed).
#
#   A deliberate copy of scripts/score_v3_sentinel.R rather than a refactor of
#   it: that script is an audit artifact committed before the V3 results
#   existed, and editing it now would weaken its provenance. The duplication is
#   the point.
#
# ONE ADAPTATION, STATED OPENLY
#   G7 as written reads "every fold's metrics_*.tsv records
#   feature_rank=within_tissue" - that is V3's string, because docs/27 s5 was
#   written for V3 and s6 reused the gates by reference. The INTENT of G7 is
#   integrity: that each fold self-reports the protocol it actually ran. For V4
#   the correct string is "lineage_penalized". Scored that way. No threshold is
#   changed; this is the only textual adaptation and it is not discretionary.
#
# INPUTS
#   results/v4_sentinel/folds/{predictions,metrics}_{BRCA,KICH,THCA,UCEC}.tsv
#   results/loco_run01/loco_predictions.tsv   (baseline comparator)
# OUTPUTS
#   results/v4_sentinel/gate_scorecard.tsv
#   results/v4_sentinel/per_tissue_comparison.tsv
#
# THCA is a negative control (near-constant within-tissue target). Reported,
# excluded from every macro statistic used to PASS a gate, never a win.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

SENT    <- c("BRCA", "KICH", "THCA", "UCEC")
CONTROL <- "THCA"
V4_DIR  <- "results/v4_sentinel"
EXPECTED_RANK <- "lineage_penalized"

# --- THRESHOLDS: docs/27 section 5, applied to V4 by section 6. DO NOT EDIT. -
TH <- list(
  G1_macro_r_exTHCA_min = 0.7096,
  G2_BRCA_r_min         = 0.5823,
  G3_tissue_r2_pred_max = 0.4919,
  G4_pooled_MAE_max     = 12.143,
  G5_mean_abs_bias_max  = 5.674,
  G6_n_tissues_ok_min   = 2L,
  G6_per_tissue_tol     = 0.050
)
BASELINE_R <- c(BRCA = 0.63232759, KICH = 0.88805065, UCEC = 0.75844945)

summarise_set <- function(d) {
  per <- d[, .(n = .N,
               r    = cor(actual, pred_clip),
               rho  = cor(actual, pred_clip, method = "spearman"),
               MAE  = mean(abs(actual - pred_clip)),
               bias = mean(pred_clip - actual)),
           by = cancer_type][order(cancer_type)]
  keep <- per[cancer_type != CONTROL]
  list(per = per,
       macro_r_exTHCA = mean(keep$r),
       macro_rho_exTHCA = mean(keep$rho),
       pooled_MAE = mean(abs(d$actual - d$pred_clip)),
       mean_abs_bias = mean(abs(per$bias)),
       r2_tissue_pred = summary(lm(pred_clip ~ factor(cancer_type), data = d))$r.squared,
       r2_tissue_truth = summary(lm(actual ~ factor(cancer_type), data = d))$r.squared)
}

pf <- file.path(V4_DIR, "folds", paste0("predictions_", SENT, ".tsv"))
if (any(!file.exists(pf)))
  stop("Sentinel incomplete; absent:\n  ", paste(pf[!file.exists(pf)], collapse = "\n  "))
v4 <- rbindlist(lapply(pf, fread), fill = TRUE)
v4[, pred_clip := pmax(predicted_reference_HRDsum, 0)]

mf <- file.path(V4_DIR, "folds", paste0("metrics_", SENT, ".tsv"))
mm <- rbindlist(lapply(mf, fread), fill = TRUE)
if (!"feature_rank" %in% names(mm)) stop("No feature_rank column; provenance unverifiable.")
fr_ok <- all(mm$feature_rank == EXPECTED_RANK)
tt_ok <- all(mm$target_transform == "identity")

b <- fread("results/loco_run01/loco_predictions.tsv")[cancer_type %in% SENT]
b[, pred_clip := pmax(predicted_reference_HRDsum, 0)]
if (!identical(sort(v4$sample_id), sort(b$sample_id)))
  stop("V4 and baseline cohorts differ; comparison invalid.")

S4 <- summarise_set(v4); S0 <- summarise_set(b)
cmp <- merge(S0$per, S4$per, by = "cancer_type", suffixes = c("_base", "_v4"))
cmp[, `:=`(d_r = r_v4 - r_base, d_MAE = MAE_v4 - MAE_base,
           d_absbias = abs(bias_v4) - abs(bias_base))]

g6_ok <- sum(sapply(names(BASELINE_R), function(t)
  cmp[cancer_type == t]$r_v4 >= BASELINE_R[[t]] - TH$G6_per_tissue_tol))

gates <- data.table(
  gate = paste0("G", 1:7),
  description = c("Macro within-tissue Pearson, ex-THCA", "BRCA within-tissue Pearson",
                  "Tissue R2 of predictions (MECHANISM GATE)", "Pooled MAE",
                  "Mean absolute tissue bias", "Non-THCA tissues not degraded (of 3)",
                  "Provenance: feature_rank + target_transform"),
  threshold = c(sprintf(">= %.4f", TH$G1_macro_r_exTHCA_min),
                sprintf(">= %.4f", TH$G2_BRCA_r_min),
                sprintf("<= %.4f", TH$G3_tissue_r2_pred_max),
                sprintf("<= %.4f", TH$G4_pooled_MAE_max),
                sprintf("<= %.4f", TH$G5_mean_abs_bias_max),
                sprintf(">= %d of 3", TH$G6_n_tissues_ok_min), "all TRUE"),
  baseline = c(sprintf("%.4f", S0$macro_r_exTHCA),
               sprintf("%.4f", S0$per[cancer_type == "BRCA"]$r),
               sprintf("%.4f", S0$r2_tissue_pred), sprintf("%.4f", S0$pooled_MAE),
               sprintf("%.4f", S0$mean_abs_bias), "3 of 3", "-"),
  observed = c(sprintf("%.4f", S4$macro_r_exTHCA),
               sprintf("%.4f", S4$per[cancer_type == "BRCA"]$r),
               sprintf("%.4f", S4$r2_tissue_pred), sprintf("%.4f", S4$pooled_MAE),
               sprintf("%.4f", S4$mean_abs_bias), sprintf("%d of 3", g6_ok),
               sprintf("rank=%s tt=%s", fr_ok, tt_ok)),
  pass = c(S4$macro_r_exTHCA >= TH$G1_macro_r_exTHCA_min,
           S4$per[cancer_type == "BRCA"]$r >= TH$G2_BRCA_r_min,
           S4$r2_tissue_pred <= TH$G3_tissue_r2_pred_max,
           S4$pooled_MAE <= TH$G4_pooled_MAE_max,
           S4$mean_abs_bias <= TH$G5_mean_abs_bias_max,
           g6_ok >= TH$G6_n_tissues_ok_min,
           fr_ok && tt_ok)
)

fwrite(gates, file.path(V4_DIR, "gate_scorecard.tsv"), sep = "\t")
fwrite(cmp,   file.path(V4_DIR, "per_tissue_comparison.tsv"), sep = "\t")

cat("\n=== V4 LINEAGE-PENALIZED SENTINEL (gates from docs/27 s5 via s6) ===\n\n")
print(gates[, .(gate, threshold, baseline, observed, pass)])
cat("\n--- Per-tissue, run 01 baseline vs V4 (THCA = negative control) ---\n")
print(cmp[, .(cancer_type, n = n_base, r_base, r_v4, d_r, MAE_base, MAE_v4, bias_base, bias_v4)])
cat("\n--- lambda_pen selected per fold (inner-fold selection) ---\n")
print(mm[, .(cancer_type, lambda_pen, alpha, lambda_at_boundary)])
cat(sprintf("\nLineage imprinting: tissue R2 of PREDICTION  baseline %.4f -> V4 %.4f\n",
            S0$r2_tissue_pred, S4$r2_tissue_pred))
cat(sprintf("                    tissue R2 of TRUTH       %.4f\n", S0$r2_tissue_truth))
cat(sprintf("\nGates passed: %d of 7\n", sum(gates$pass)))
cat(if (all(gates$pass))
  "\nVERDICT: PASSED all gates.\n" else
  paste0("\nVERDICT: FAILED ", paste(gates[pass == FALSE]$gate, collapse = ", "),
         ".\nPer docs/27 section 6: terminate this branch. No further iterative variants.\n"))

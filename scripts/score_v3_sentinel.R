# =============================================================================
# scripts/score_v3_sentinel.R - Mechanical scorer for the V3-abs sentinel gates
# =============================================================================
#
# USAGE
#   Rscript scripts/score_v3_sentinel.R
#
# WHY THIS SCRIPT EXISTS
#   The gates G1-G7 are pre-registered in docs/27_V3_PREREGISTRATION.md section 5,
#   committed as d72695c at 10:27 CDT while results/v3_sentinel/folds/ was empty
#   and LSF 323220006 was still RUNning. This script hard-codes those thresholds
#   as constants and prints a verdict. It was written BEFORE any fold finished,
#   for one reason: so that scoring is arithmetic rather than judgement applied
#   after seeing the numbers.
#
#   Do not edit the THRESHOLDS block. If a threshold here ever disagrees with
#   docs/27 section 5, docs/27 wins and this script is the bug.
#
# INPUTS
#   results/v3_sentinel/folds/predictions_{BRCA,KICH,THCA,UCEC}.tsv
#   results/v3_sentinel/folds/metrics_{BRCA,KICH,THCA,UCEC}.tsv
#   results/loco_run01/loco_predictions.tsv          (the baseline comparator)
#
# OUTPUTS
#   results/v3_sentinel/gate_scorecard.tsv           one row per gate
#   results/v3_sentinel/per_tissue_comparison.tsv    baseline vs V3-abs
#   printed verdict: AUTHORISED or NOT AUTHORISED
#
# NOTE ON THCA
#   THCA is a negative control (near-constant within-tissue target, run-01
#   r = -0.036). It is reported but excluded from every macro statistic used to
#   PASS a gate. It cannot constitute a win. See docs/27 section 4.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

SENT <- c("BRCA", "KICH", "THCA", "UCEC")
CONTROL <- "THCA"
V3_DIR <- "results/v3_sentinel"

# --- THRESHOLDS: frozen in docs/27 section 5. DO NOT EDIT. ------------------
TH <- list(
  G1_macro_r_exTHCA_min   = 0.7096,   # baseline 0.7596 - 0.050
  G2_BRCA_r_min           = 0.5823,   # baseline 0.6323 - 0.050
  G3_tissue_r2_pred_max   = 0.4919,   # baseline 0.5219 - 0.030  (MECHANISM GATE)
  G4_pooled_MAE_max       = 12.143,   # baseline 11.039 * 1.10
  G5_mean_abs_bias_max    = 5.674,    # baseline 5.158 * 1.10
  G6_n_tissues_ok_min     = 2L,       # of BRCA/KICH/UCEC, each losing <= 0.05 r
  G6_per_tissue_tol       = 0.050
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
  list(
    per = per,
    macro_r_exTHCA   = mean(keep$r),
    macro_rho_exTHCA = mean(keep$rho),
    macro_r_all      = mean(per$r),
    pooled_MAE       = mean(abs(d$actual - d$pred_clip)),
    mean_abs_bias    = mean(abs(per$bias)),
    # Lineage imprinting: how much of the PREDICTION is explained by tissue
    # identity alone. The mechanism gate. Compare with the same quantity for
    # the truth - a model should not encode more lineage than its target does.
    r2_tissue_pred   = summary(lm(pred_clip ~ factor(cancer_type), data = d))$r.squared,
    r2_tissue_truth  = summary(lm(actual   ~ factor(cancer_type), data = d))$r.squared
  )
}

# --- Load V3-abs ------------------------------------------------------------
pf <- file.path(V3_DIR, "folds", paste0("predictions_", SENT, ".tsv"))
missing <- pf[!file.exists(pf)]
if (length(missing)) {
  stop("Sentinel incomplete; these fold outputs are absent:\n  ",
       paste(missing, collapse = "\n  "),
       "\nRefusing to score a partial sentinel.")
}
v3 <- rbindlist(lapply(pf, fread), fill = TRUE)
v3[, pred_clip := pmax(predicted_reference_HRDsum, 0)]   # C2 adopted rule

# G7a: every fold must self-report the protocol it actually ran under.
mf <- file.path(V3_DIR, "folds", paste0("metrics_", SENT, ".tsv"))
mm <- rbindlist(lapply(mf, fread), fill = TRUE)
if (!"feature_rank" %in% names(mm)) {
  stop("metrics_*.tsv carries no feature_rank column; provenance unverifiable.")
}
fr_ok <- all(mm$feature_rank == "within_tissue")
tt_ok <- all(mm$target_transform == "identity")

# --- Load baseline, same tissues, same clipping -----------------------------
b <- fread("results/loco_run01/loco_predictions.tsv")[cancer_type %in% SENT]
b[, pred_clip := pmax(predicted_reference_HRDsum, 0)]

stopifnot(nrow(v3) == nrow(b))
setkey(v3, sample_id); setkey(b, sample_id)
if (!identical(sort(v3$sample_id), sort(b$sample_id))) {
  stop("V3 and baseline sentinel cohorts differ; comparison would be invalid.")
}

S3 <- summarise_set(v3)
S0 <- summarise_set(b)

cmp <- merge(S0$per, S3$per, by = "cancer_type", suffixes = c("_base", "_v3"))
cmp[, `:=`(d_r = r_v3 - r_base, d_rho = rho_v3 - rho_base,
           d_MAE = MAE_v3 - MAE_base, d_absbias = abs(bias_v3) - abs(bias_base))]

# --- Score ------------------------------------------------------------------
g6_ok <- sum(sapply(names(BASELINE_R), function(t)
  cmp[cancer_type == t]$r_v3 >= BASELINE_R[[t]] - TH$G6_per_tissue_tol))

gates <- data.table(
  gate = paste0("G", 1:7),
  description = c(
    "Macro within-tissue Pearson, ex-THCA",
    "BRCA within-tissue Pearson",
    "Tissue R2 of predictions (MECHANISM GATE)",
    "Pooled MAE",
    "Mean absolute tissue bias",
    "Non-THCA tissues not individually degraded (of 3)",
    "Provenance: feature_rank + target_transform recorded"),
  threshold = c(
    sprintf(">= %.4f", TH$G1_macro_r_exTHCA_min),
    sprintf(">= %.4f", TH$G2_BRCA_r_min),
    sprintf("<= %.4f", TH$G3_tissue_r2_pred_max),
    sprintf("<= %.4f", TH$G4_pooled_MAE_max),
    sprintf("<= %.4f", TH$G5_mean_abs_bias_max),
    sprintf(">= %d of 3", TH$G6_n_tissues_ok_min),
    "all TRUE"),
  baseline = c(
    sprintf("%.4f", S0$macro_r_exTHCA), sprintf("%.4f", S0$per[cancer_type == "BRCA"]$r),
    sprintf("%.4f", S0$r2_tissue_pred), sprintf("%.4f", S0$pooled_MAE),
    sprintf("%.4f", S0$mean_abs_bias), "3 of 3", "-"),
  observed = c(
    sprintf("%.4f", S3$macro_r_exTHCA), sprintf("%.4f", S3$per[cancer_type == "BRCA"]$r),
    sprintf("%.4f", S3$r2_tissue_pred), sprintf("%.4f", S3$pooled_MAE),
    sprintf("%.4f", S3$mean_abs_bias), sprintf("%d of 3", g6_ok),
    sprintf("feature_rank=%s target_transform=%s", fr_ok, tt_ok)),
  pass = c(
    S3$macro_r_exTHCA >= TH$G1_macro_r_exTHCA_min,
    S3$per[cancer_type == "BRCA"]$r >= TH$G2_BRCA_r_min,
    S3$r2_tissue_pred <= TH$G3_tissue_r2_pred_max,
    S3$pooled_MAE <= TH$G4_pooled_MAE_max,
    S3$mean_abs_bias <= TH$G5_mean_abs_bias_max,
    g6_ok >= TH$G6_n_tissues_ok_min,
    fr_ok && tt_ok)
)

authorised <- all(gates$pass)

dir.create(V3_DIR, showWarnings = FALSE, recursive = TRUE)
fwrite(gates, file.path(V3_DIR, "gate_scorecard.tsv"), sep = "\t")
fwrite(cmp,   file.path(V3_DIR, "per_tissue_comparison.tsv"), sep = "\t")

cat("\n=== V3-abs SENTINEL SCORECARD (gates frozen in docs/27 s5, commit d72695c) ===\n\n")
print(gates[, .(gate, threshold, baseline, observed, pass)])
cat("\n--- Per-tissue, baseline vs V3-abs (THCA = negative control) ---\n")
print(cmp[, .(cancer_type, n = n_base, r_base, r_v3, d_r, MAE_base, MAE_v3,
              bias_base, bias_v3)])
cat(sprintf("\nLineage imprinting: tissue R2 of PREDICTION  baseline %.4f -> V3-abs %.4f\n",
            S0$r2_tissue_pred, S3$r2_tissue_pred))
cat(sprintf("                    tissue R2 of TRUTH       %.4f (unchanged by construction)\n",
            S0$r2_tissue_truth))
cat(sprintf("\nGates passed: %d of 7\n", sum(gates$pass)))
cat(if (authorised) {
  "\nVERDICT: AUTHORISED - full 30-fold V3-abs LOCO array may be submitted.\n"
} else {
  paste0("\nVERDICT: NOT AUTHORISED - failed ",
         paste(gates[pass == FALSE]$gate, collapse = ", "),
         ".\nPer docs/27 section 5: STOP V3-abs. No re-tuning against these four tissues.\n")
})
cat("\nWritten: ", file.path(V3_DIR, "gate_scorecard.tsv"), "\n", sep = "")

# =============================================================================
# scripts/score_candidate_selection.R - Mechanical scorer for S1-S7
# =============================================================================
#
# USAGE
#   Rscript scripts/score_candidate_selection.R
#
# WHAT THIS DOES
#   Scores the complete 30-fold V3-abs LOCO run (Candidate B) against the
#   selection rule S1-S7 pre-registered in
#   docs/28_CANDIDATE_SELECTION_AND_CNS_PROTOCOL.md section 4, committed as
#   28bb913 at 11:58 CDT while results/v3_loco_full/ was still empty.
#
#   The thresholds below are transcribed from that table and must not be
#   edited. If this script and docs/28 ever disagree, docs/28 wins and this
#   script is the bug.
#
#   Both candidates are scored under the adopted C2 clipping, with per-tissue
#   statistics computed then macro-averaged - the same estimator for both, so
#   the comparison is like-for-like. Note this macro estimator is NOT the
#   pooled-within estimator that produced the 0.612 headline in docs/24; the
#   two are not interchangeable.
#
# INPUTS
#   results/loco_run01/loco_predictions.tsv    Candidate A (pooled)
#   results/v3_loco_full/folds/predictions_*.tsv  Candidate B (within_tissue)
#   results/v3_loco_full/folds/metrics_*.tsv      provenance for S7
#
# OUTPUTS
#   results/v3_loco_full/selection_scorecard.tsv
#   results/v3_loco_full/per_tissue_A_vs_B.tsv
#   results/tables/tableA_model_comparison.{tsv,md}   (Table A, filled in)
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# --- THRESHOLDS: docs/28 section 4. DO NOT EDIT. ---------------------------
TH <- list(
  S1_macro_r_min        = 0.4697,   # A - 0.05
  S2_tissue_r2_max      = 0.5084,   # A - 0.05, AND strictly below A
  S3_mean_abs_bias_max  = 3.4628,   # no worse than A
  S4_pooled_MAE_max     = 9.8689,   # A * 1.10
  S5_n_pos_r_min        = 27L,      # of 30
  S6_n_tissues_min      = 20L,      # tissues where B >= A - 0.05
  S6_tol                = 0.05
)
TIE_MARGIN <- 0.01   # docs/28 "Tie-breaking and simplicity"

summarise_run <- function(d) {
  per <- d[, .(n = .N,
               r    = cor(actual, pred_clip),
               rho  = cor(actual, pred_clip, method = "spearman"),
               MAE  = mean(abs(actual - pred_clip)),
               bias = mean(pred_clip - actual)),
           by = cancer_type][order(cancer_type)]
  list(per = per,
       macro_r        = mean(per$r),
       macro_rho      = mean(per$rho),
       pooled_MAE     = mean(abs(d$actual - d$pred_clip)),
       mean_abs_bias  = mean(abs(per$bias)),
       r2_tissue_pred = summary(lm(pred_clip ~ factor(cancer_type), data = d))$r.squared,
       r2_tissue_truth= summary(lm(actual   ~ factor(cancer_type), data = d))$r.squared,
       n_pos_r        = sum(per$r > 0))
}

# --- Candidate A ------------------------------------------------------------
A <- fread("results/loco_run01/loco_predictions.tsv")
A[, pred_clip := pmax(predicted_reference_HRDsum, 0)]

# --- Candidate B ------------------------------------------------------------
pf <- Sys.glob("results/v3_loco_full/folds/predictions_*.tsv")
if (length(pf) != 30L)
  stop("S7 fails: expected 30 fold predictions, found ", length(pf))
B <- rbindlist(lapply(pf, fread), fill = TRUE)
B[, pred_clip := pmax(predicted_reference_HRDsum, 0)]

mf <- Sys.glob("results/v3_loco_full/folds/metrics_*.tsv")
mm <- rbindlist(lapply(mf, fread), fill = TRUE)
s7_rank <- all(mm$feature_rank == "within_tissue")
s7_tt   <- all(mm$target_transform == "identity")
s7_n    <- length(mf) == 30L

for (nm in c("GBM", "LGG")) {
  if (nm %in% A$cancer_type || nm %in% B$cancer_type)
    stop("CNS rows present in a development comparison - refusing to score.")
}
if (!setequal(A$sample_id, B$sample_id))
  stop("Candidate A and B cohorts differ; comparison would be invalid.")

SA <- summarise_run(A); SB <- summarise_run(B)

cmp <- merge(SA$per, SB$per, by = "cancer_type", suffixes = c("_A", "_B"))
cmp[, `:=`(d_r = r_B - r_A, d_rho = rho_B - rho_A,
           d_MAE = MAE_B - MAE_A, d_absbias = abs(bias_B) - abs(bias_A))]
setorder(cmp, -d_r)

s6_n <- sum(cmp$r_B >= cmp$r_A - TH$S6_tol)

# --- Score ------------------------------------------------------------------
res <- data.table(
  criterion = paste0("S", 1:7),
  description = c(
    "Macro within-tissue Pearson",
    "Tissue R2 of predictions (lineage imprinting)",
    "Mean absolute tissue bias (C1)",
    "Pooled MAE",
    "Tissues with positive within-tissue r",
    "Breadth: B >= A - 0.05 across tissues",
    "30 folds, feature_rank + transform recorded"),
  threshold = c(
    sprintf(">= %.4f", TH$S1_macro_r_min),
    sprintf("<= %.4f and < A", TH$S2_tissue_r2_max),
    sprintf("<= %.4f", TH$S3_mean_abs_bias_max),
    sprintf("<= %.4f", TH$S4_pooled_MAE_max),
    sprintf(">= %d of 30", TH$S5_n_pos_r_min),
    sprintf(">= %d of 30", TH$S6_n_tissues_min),
    "all TRUE"),
  candidate_A = c(
    sprintf("%.4f", SA$macro_r), sprintf("%.4f", SA$r2_tissue_pred),
    sprintf("%.4f", SA$mean_abs_bias), sprintf("%.4f", SA$pooled_MAE),
    sprintf("%d of 30", SA$n_pos_r), "-", "-"),
  candidate_B = c(
    sprintf("%.4f", SB$macro_r), sprintf("%.4f", SB$r2_tissue_pred),
    sprintf("%.4f", SB$mean_abs_bias), sprintf("%.4f", SB$pooled_MAE),
    sprintf("%d of 30", SB$n_pos_r), sprintf("%d of 30", s6_n),
    sprintf("n=%s rank=%s tt=%s", s7_n, s7_rank, s7_tt)),
  pass = c(
    SB$macro_r >= TH$S1_macro_r_min,
    SB$r2_tissue_pred <= TH$S2_tissue_r2_max && SB$r2_tissue_pred < SA$r2_tissue_pred,
    SB$mean_abs_bias <= TH$S3_mean_abs_bias_max,
    SB$pooled_MAE <= TH$S4_pooled_MAE_max,
    SB$n_pos_r >= TH$S5_n_pos_r_min,
    s6_n >= TH$S6_n_tissues_min,
    s7_n && s7_rank && s7_tt)
)

all_pass <- all(res$pass)
# Tie-break: if B passes everything but is within +/- TIE_MARGIN of A on the
# two headline criteria, A is retained on simplicity grounds.
tie <- all_pass &&
  abs(SB$macro_r - SA$macro_r) <= TIE_MARGIN &&
  abs(SB$r2_tissue_pred - SA$r2_tissue_pred) <= TIE_MARGIN
primary <- if (all_pass && !tie) "B (V3-abs, within_tissue)" else "A (run 01, pooled)"

dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
fwrite(res, "results/v3_loco_full/selection_scorecard.tsv", sep = "\t")
fwrite(cmp, "results/v3_loco_full/per_tissue_A_vs_B.tsv", sep = "\t")

# --- Table A ----------------------------------------------------------------
tabA <- data.table(
  model = c("A: run 01 baseline", "B: V3-abs"),
  feature_strategy = c("Pooled total variance (unsupervised)",
                       "Pooled within-tissue variance (unsupervised)"),
  n_features = c(5000L, 5000L),
  macro_pearson  = round(c(SA$macro_r, SB$macro_r), 4),
  macro_spearman = round(c(SA$macro_rho, SB$macro_rho), 4),
  pooled_MAE     = round(c(SA$pooled_MAE, SB$pooled_MAE), 4),
  mean_abs_tissue_bias = round(c(SA$mean_abs_bias, SB$mean_abs_bias), 4),
  tissue_R2_predicted  = round(c(SA$r2_tissue_pred, SB$r2_tissue_pred), 4),
  tissue_R2_observed   = round(c(SA$r2_tissue_truth, SB$r2_tissue_truth), 4),
  n_tissues_positive_r = c(SA$n_pos_r, SB$n_pos_r),
  notes = c("Incumbent; independently replicated in pipeline_glmnet/",
            sprintf("Selection rule S1-S7: %d of 7 passed", sum(res$pass)))
)
fwrite(tabA, "results/tables/tableA_model_comparison.tsv", sep = "\t")

md <- c("# Table A - Source-domain model comparison",
        "",
        sprintf("Development LOCO, 7,065 samples, 30 non-CNS cancer types. Scored %s.",
                format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
        "Per-tissue statistics computed under the adopted C2 clipping, then macro-averaged.",
        "This macro estimator is not interchangeable with the pooled-within estimator (0.612) in docs/24.",
        "",
        paste("|", paste(names(tabA), collapse = " | "), "|"),
        paste("|", paste(rep("---", ncol(tabA)), collapse = " | "), "|"),
        apply(tabA, 1, function(r) paste("|", paste(r, collapse = " | "), "|")),
        "",
        sprintf("**Primary candidate: %s** (S1-S7: %d of 7 passed).", primary, sum(res$pass)))
writeLines(md, "results/tables/tableA_model_comparison.md")

# --- Report -----------------------------------------------------------------
cat("\n=== S1-S7 SELECTION SCORECARD (rule frozen in docs/28 s4, commit 28bb913) ===\n\n")
print(res[, .(criterion, threshold, candidate_A, candidate_B, pass)])

cat("\n--- Lineage imprinting ---\n")
cat(sprintf("  Tissue R2 of PREDICTION : A %.4f -> B %.4f\n", SA$r2_tissue_pred, SB$r2_tissue_pred))
cat(sprintf("  Tissue R2 of TRUTH      : %.4f (a model should not exceed this by much)\n",
            SA$r2_tissue_truth))
cat(sprintf("  Excess over truth       : A +%.4f -> B %+.4f\n",
            SA$r2_tissue_pred - SA$r2_tissue_truth,
            SB$r2_tissue_pred - SB$r2_tissue_truth))

cat("\n--- Largest per-tissue gains (B - A, within-tissue Pearson) ---\n")
print(head(cmp[, .(cancer_type, n = n_A, r_A, r_B, d_r, bias_A, bias_B)], 5))
cat("\n--- Largest per-tissue losses ---\n")
print(tail(cmp[, .(cancer_type, n = n_A, r_A, r_B, d_r, bias_A, bias_B)], 5))

cat(sprintf("\nCriteria passed: %d of 7\n", sum(res$pass)))
if (tie) cat("TIE-BREAK: margins within +/-0.01; docs/28 retains A on simplicity.\n")
cat(sprintf("\nPRIMARY CANDIDATE: %s\n", primary))
cat("\nWritten:\n  results/v3_loco_full/selection_scorecard.tsv",
    "\n  results/v3_loco_full/per_tissue_A_vs_B.tsv",
    "\n  results/tables/tableA_model_comparison.{tsv,md}\n")

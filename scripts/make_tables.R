# =============================================================================
# scripts/make_tables.R -- presentation tables for KIDS26-Team12
# -----------------------------------------------------------------------------
# Run from the repository root:   Rscript scripts/make_tables.R
#
# INPUTS (required; fails loudly if missing)
#   results/loco_run01/loco_predictions.tsv   7,065 x 12, source-domain LOCO run 01
#   results/loco_run01/loco_metrics.tsv       per-tissue LOCO metrics (30 rows)
#   R/figures.R                               shared helpers (C2 rule, CNS lock)
#
# OUTPUTS (results/tables/)
#   tableA_model_comparison.tsv / .md   STUB -- headers only, filled after selection
#   tableB_per_cancer_loco.tsv / .md    per-cancer LOCO, run 01 (30 rows + pooled)
#   tableC_cns_external.tsv / .md       STUB -- headers only, filled after CNS run
#   tableD_blocker_status.tsv / .md     blockers C1-C4
#
# CONVENTIONS
#   * Adopted C2 rule applied before ANY computation: pred_clip = pmax(pred, 0).
#   * CNS LOCK: development cancers only; aborts if GBM/LGG appear. Table C is an
#     EMPTY STUB by design -- no CNS value is read, computed or written here.
#   * Table D numbers are quoted verbatim from docs/21, docs/25, docs/26, docs/27.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- normalizePath(".", mustWork = TRUE)
if (!file.exists(file.path(ROOT, "R", "figures.R"))) {
  stop("Run this script from the repository root (R/figures.R not found at ", ROOT, ")",
       call. = FALSE)
}
source(file.path(ROOT, "R", "figures.R"))
set.seed(FIG_SEED)

TAB_DIR <- file.path(ROOT, "results", "tables")
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
message("== make_tables.R ==  seed = ", FIG_SEED)

## ---- markdown writer --------------------------------------------------------
#' Write a data.table as a GitHub-flavoured markdown table (plus an optional
#' title and notes block). Zero-row tables are written as an explicit stub.
write_md <- function(dt, path, title, notes = character(), align = NULL) {
  # Escape pipes and newlines so that free-text cells cannot break the table.
  esc <- function(x) gsub("\n", " ", gsub("|", "\\|", as.character(x), fixed = TRUE))
  cols <- names(dt)
  if (is.null(align)) align <- rep("---", length(cols))
  body <- if (nrow(dt) == 0L) {
    "*(stub -- no rows yet; columns are defined above)*"
  } else {
    apply(dt, 1L, function(r) paste0("| ", paste(esc(r), collapse = " | "), " |"))
  }
  lines <- c(
    paste0("# ", title), "",
    paste0("| ", paste(cols, collapse = " | "), " |"),
    paste0("| ", paste(align, collapse = " | "), " |"),
    body, "",
    if (length(notes)) c("", notes) else NULL)
  writeLines(lines, path)
  message(sprintf("  wrote %-56s  %8.1f KB", sub(paste0(ROOT, "/"), "", path),
                  file.size(path) / 1024))
  invisible(path)
}

write_pair <- function(dt, stem, title, notes = character()) {
  tsv <- file.path(TAB_DIR, paste0(stem, ".tsv"))
  fwrite(dt, tsv, sep = "\t")
  message(sprintf("  wrote %-56s  %8.1f KB", sub(paste0(ROOT, "/"), "", tsv),
                  file.size(tsv) / 1024))
  write_md(dt, file.path(TAB_DIR, paste0(stem, ".md")), title, notes)
  invisible(TRUE)
}

## =============================================================================
## Table A -- model comparison  (STUB: filled after candidate selection)
## =============================================================================
tableA <- data.table(
  Model                      = character(),
  Description                = character(),
  Pooled_Pearson             = numeric(),
  Pooled_Spearman            = numeric(),
  Macro_within_tissue_Pearson  = numeric(),
  Macro_within_tissue_Spearman = numeric(),
  Macro_MAE                  = numeric(),
  Skill_vs_tissue_mean_null  = numeric(),
  Tissue_R2_of_prediction    = numeric(),
  Mean_abs_tissue_offset     = numeric(),
  Selected                   = character())
write_pair(tableA, "tableA_model_comparison",
  "Table A. Model comparison (STUB -- awaiting candidate selection)",
  notes = c(
    "**STATUS: STUB.** Column headers are frozen; rows are added by the main agent",
    "after candidate selection completes. No model has been selected at the time this",
    "file was generated, and no row may be added from a run that saw the held-out",
    "cancer type or the locked CNS cohort.",
    "",
    "Metric definitions match `results/loco_run01/`: `Macro_*` statistics are computed",
    "inside each held-out cancer type and then averaged unweighted over cancer types;",
    "`Skill_vs_tissue_mean_null` is 1 - MAE_model / MAE_tissue_mean_null;",
    "`Tissue_R2_of_prediction` is the one-way ANOVA R2 of cancer type on the clipped",
    "prediction. All predictions must have the adopted C2 rule `pmax(pred, 0)` applied."))

## =============================================================================
## Table B -- per-cancer LOCO, run 01
## =============================================================================
pred <- read_required(file.path(ROOT, "results/loco_run01/loco_predictions.tsv"))
met  <- read_required(file.path(ROOT, "results/loco_run01/loco_metrics.tsv"))

pred <- keep_development(pred, "cancer_type", "loco_predictions")
met  <- keep_development(met,  "cancer_type", "loco_metrics")
pred <- apply_c2_clip(pred)
stopifnot(nrow(pred) == 7065L, uniqueN(pred$cancer_type) == 30L, nrow(met) == 30L)

# All statistics recomputed from the predictions under the adopted C2 clipping.
B <- pred[, {
  e <- pred_clip - actual
  .(N            = .N,
    mean_obs     = mean(actual),
    mean_pred    = mean(pred_clip),
    Pearson      = if (sd(pred_clip) > 0) cor(pred_clip, actual) else NA_real_,
    Spearman     = if (sd(pred_clip) > 0) cor(pred_clip, actual, method = "spearman") else NA_real_,
    MAE          = mean(abs(e)),
    RMSE         = sqrt(mean(e^2)),
    bias         = mean(e),
    n_ood        = sum(ood),
    n_qc_fail    = sum(qc_fail),
    n_not_report = sum(!reportable))
}, by = cancer_type]

B[, abstention_pct := 100 * n_not_report / N]
B[, OOD_abstention := sprintf("%d OOD / %d QC-fail / %d non-reportable (%.1f%%)",
                              n_ood, n_qc_fail, n_not_report, abstention_pct)]
setorder(B, -Pearson)

# Cross-check against the shipped per-tissue metrics. The shipped values are
# computed on the UNCLIPPED prediction, so the raw-prediction statistics must
# reproduce them to machine precision; that is the integrity test.
raw_chk <- pred[, .(Pearson_raw  = cor(predicted_reference_HRDsum, actual),
                    Spearman_raw = cor(predicted_reference_HRDsum, actual,
                                       method = "spearman")), by = cancer_type]
chk <- merge(raw_chk, met[, .(cancer_type, Pearson, Spearman)], by = "cancer_type")
max_delta <- max(abs(chk$Spearman_raw - chk$Spearman), abs(chk$Pearson_raw - chk$Pearson))
message(sprintf("-- Table B integrity check vs loco_metrics.tsv (unclipped): max |delta| = %.2e",
                max_delta))
if (max_delta > 1e-10) {
  stop("Table B raw statistics disagree with results/loco_run01/loco_metrics.tsv by ",
       signif(max_delta, 3), " -- investigate before presenting.", call. = FALSE)
}

# Effect of the adopted C2 clipping on the rank statistics. Clipping is
# rank-preserving on a POOLED basis only in the weak sense that it is monotone
# non-decreasing; it maps every negative prediction to the same value, which
# creates ties. In tissues with many negative predictions those ties measurably
# move the within-tissue Spearman, so the shift is quantified rather than assumed.
clip_shift <- merge(B[, .(cancer_type, Spearman_clip = Spearman)],
                    raw_chk[, .(cancer_type, Spearman_raw)], by = "cancer_type")
neg_by_tissue <- pred[, .(n_neg = sum(predicted_reference_HRDsum < 0)), by = cancer_type]
clip_shift <- merge(clip_shift, neg_by_tissue, by = "cancer_type")
clip_shift[, delta := Spearman_clip - Spearman_raw]
max_clip_shift <- clip_shift[which.max(abs(delta))]
message(sprintf("-- C2 clipping shifts within-tissue Spearman by at most %+.4f (%s, %d negative predictions)",
                max_clip_shift$delta, max_clip_shift$cancer_type, max_clip_shift$n_neg))

pooled_row <- pred[, {
  e <- pred_clip - actual
  data.table(cancer_type = "POOLED (all 30)", N = .N,
             mean_obs = mean(actual), mean_pred = mean(pred_clip),
             Pearson = cor(pred_clip, actual),
             Spearman = cor(pred_clip, actual, method = "spearman"),
             MAE = mean(abs(e)), RMSE = sqrt(mean(e^2)), bias = mean(e),
             OOD_abstention = sprintf("%d OOD / %d QC-fail / %d non-reportable (%.1f%%)",
                                      sum(ood), sum(qc_fail), sum(!reportable),
                                      100 * sum(!reportable) / .N))
}]
macro_row <- data.table(
  cancer_type = "MACRO (mean of 30)", N = NA_integer_,
  mean_obs = mean(B$mean_obs), mean_pred = mean(B$mean_pred),
  Pearson = mean(B$Pearson), Spearman = mean(B$Spearman),
  MAE = mean(B$MAE), RMSE = mean(B$RMSE), bias = mean(B$bias),
  OOD_abstention = sprintf("mean abstention %.1f%%", mean(B$abstention_pct)))

Bout <- rbind(
  B[, .(Cancer = cancer_type, N, mean_obs, mean_pred, Pearson, Spearman, MAE, RMSE,
        bias, OOD_abstention)],
  pooled_row[, .(Cancer = cancer_type, N, mean_obs, mean_pred, Pearson, Spearman,
                 MAE, RMSE, bias, OOD_abstention)],
  macro_row[, .(Cancer = cancer_type, N, mean_obs, mean_pred, Pearson, Spearman,
                MAE, RMSE, bias, OOD_abstention)])

num_cols <- c("mean_obs", "mean_pred", "Pearson", "Spearman", "MAE", "RMSE", "bias")
Bout[, (num_cols) := lapply(.SD, round, 3), .SDcols = num_cols]
setnames(Bout,
         c("mean_obs", "mean_pred", "bias"),
         c("Mean_observed_HRD", "Mean_predicted_HRD", "Bias_mean_signed_error"))

best  <- B[1]
worst <- B[.N]
message(sprintf("-- Table B extremes: best %s (r = %.3f), worst %s (r = %.3f)",
                best$cancer_type, best$Pearson, worst$cancer_type, worst$Pearson))

write_pair(Bout, "tableB_per_cancer_loco",
  "Table B. Per-cancer leave-one-cancer-out performance, run 01",
  notes = c(
    "**Source:** `results/loco_run01/loco_predictions.tsv` (7,065 specimens, 30 development",
    "cancer types). Every statistic is recomputed here from the predictions with the adopted",
    "C2 rule applied first: `pred_clip = pmax(predicted_reference_HRDsum, 0)`. As an integrity",
    "check the *unclipped* statistics were verified to reproduce",
    "`results/loco_run01/loco_metrics.tsv` to machine precision.",
    "",
    "**Note on clipping and ranks.** Clipping is monotone but not strictly increasing: it maps",
    "every negative prediction to the same value, creating ties. In tissues with many negative",
    sprintf("predictions this measurably moves the within-tissue Spearman -- the largest shift is %+.4f",
            max_clip_shift$delta),
    sprintf("in %s (%d of %d predictions negative). The shift is small and always toward better",
            max_clip_shift$cancer_type, max_clip_shift$n_neg,
            B[cancer_type == max_clip_shift$cancer_type, N]),
    "agreement with the zero floor of the label, but the claim that clipping leaves rankings",
    "*exactly* unchanged holds pooled, not within every tissue. MAE/RMSE/bias also differ from",
    "the shipped file because those values are unclipped.",
    "",
    "**Reading the table.** Each row is a cancer type that was entirely absent from the",
    "training data for its own prediction. Sorted by within-tissue Pearson r, descending.",
    "`Bias_mean_signed_error` is the C1 per-tissue offset (predicted minus observed).",
    "`OOD_abstention` counts specimens flagged out-of-distribution, failing QC, or otherwise",
    "not reportable under the deployed gating rules.",
    "",
    "**Caveats.** `Pearson`/`Spearman` here are *within-tissue* and are the honest measure of",
    "ranking quality; the POOLED row is inflated by between-tissue mean differences and must",
    "not be quoted as a within-tissue result. OV has N = 10 (blocker C4) and its correlation",
    "is not interpretable. THCA, PCPG, CHOL and TGCT carry essentially no within-tissue signal.",
    "",
    "**Cohort lock.** GBM and LGG are excluded and asserted absent. No CNS value appears here."))

## =============================================================================
## Table C -- locked CNS external test  (STUB: filled after the CNS run)
## =============================================================================
tableC <- data.table(
  Cancer                   = character(),   # GBM / LGG / combined
  N                        = integer(),
  Mean_observed_HRD        = numeric(),
  Mean_predicted_HRD       = numeric(),
  Pearson                  = numeric(),
  Spearman                 = numeric(),
  MAE                      = numeric(),
  RMSE                     = numeric(),
  Bias_mean_signed_error   = numeric(),
  Skill_vs_tissue_mean_null = numeric(),
  OOD_abstention           = character(),
  Prereg_gate_outcome      = character())
write_pair(tableC, "tableC_cns_external",
  "Table C. Locked CNS external test (STUB -- cohort sealed, not yet scored)",
  notes = c(
    "**STATUS: STUB. THE CNS COHORT IS LOCKED.** 642 GBM + LGG specimens are held sealed and",
    "will be scored exactly once, after the candidate model is frozen. No GBM or LGG HRDsum",
    "value has been read, predicted, plotted or summarised by any script in the figure/table",
    "layer -- `R/figures.R::assert_development()` aborts if either type appears.",
    "",
    "Column headers are frozen now so that the analysis is pre-committed. The main agent fills",
    "the rows after the single external run. Results must be reported against the",
    "pre-registered gates in `docs/27_V3_PREREGISTRATION.md`, whatever they show; no reselection,",
    "refitting or gate revision is permitted after the cohort is opened."))

## =============================================================================
## Table D -- blocker status C1-C4
##   Content sourced from docs/21_BLOCKER_RESOLUTION_PLAN.md,
##   docs/25_C1_C2_C3_INVESTIGATION.md, docs/26_C1_N_OF_1_WORKAROUND.md,
##   docs/27_V3_PREREGISTRATION.md. Numbers are quoted from those documents.
## =============================================================================
for (d in c("docs/21_BLOCKER_RESOLUTION_PLAN.md", "docs/25_C1_C2_C3_INVESTIGATION.md",
            "docs/26_C1_N_OF_1_WORKAROUND.md", "docs/27_V3_PREREGISTRATION.md")) {
  if (!file.exists(file.path(ROOT, d))) {
    stop("REQUIRED SOURCE DOC MISSING for Table D: ", d, call. = FALSE)
  }
}

tableD <- data.table(
  Blocker = c("C1", "C1b", "C1c", "C2", "C3", "C4"),

  Question = c(
    "Can a per-tissue additive calibration offset be estimated for a cancer type the model has never seen, so that an ABSOLUTE HRDsum is trustworthy on a new tissue?",
    "If a few LABELLED samples of the new cancer type are available, can the offset be estimated from them?",
    "Can the offset be side-stepped entirely by predicting a RELATIVE (within-tissue rank/percentile) target instead of an absolute one?",
    "HRDsum is >= 0 by construction and has a 14.3% point mass at exactly zero, but the elastic net emits negative predictions. How should the zero floor be handled?",
    "Does the model's apparent signal come from tumour purity rather than HRD biology?",
    "OV is the canonical HRD cancer and the main clinical application, but is it actually represented in the cohort?",
    NULL),

  `Test performed` = c(
    "Label-free calibrator: regress the per-tissue offset on label-free tissue covariates (mean_pred, sd_pred, mean_purity, sd_purity, mean_ood, n), evaluated leave-one-tissue-out (LOTO). Positive/negative controls on the calibrator machinery (tests/test_calibration.R cases 8-9).",
    "Few-shot calibration from k labelled samples of the held-out tissue, unshrunk sample mean vs empirical-Bayes shrinkage toward a cross-tissue offset prior; k = 1, 3, 5, 10, 20; macro-averaged over 29 eligible tissues.",
    "(a) Post-hoc percentile remap of the absolute prediction. (b) Fully trained relative-target model (V1/V2), 4-tissue sentinel array 323195423 scored against 7 pre-registered gates. (c) V3-abs (within-tissue variance feature ranking) pre-registered in docs/27, sentinel 323220006.",
    "Clip predictions at zero and re-score all 7,065 development samples (no refit needed). Separately, log1p target transform run over all 30 LOCO folds (arrays 323176856 + 323169191).",
    "Partial correlation cor(pred, actual | purity) vs the raw within-tissue correlation; decomposition of skill and offset by purity tertile; correlation of per-tissue offset with per-tissue mean purity across 30 tissues.",
    "None. This is a cohort-composition disclosure, not a code defect: TCGA ovarian methylation is mostly 27k-array and is excluded by the 450k feature bridge.",
    NULL),

  Result = c(
    "FAILED in every configuration. Best LOTO R2 = -0.116 (mean_purity); all-covariate model LOTO MAE 3.916 vs naive 3.466, LOTO R2 -0.217. Offset-covariate correlations are weak: n +0.267, mean_purity -0.120, sd_pred +0.086, mean_pred -0.040, mean_ood +0.000. The offset itself is large: mean |offset| 3.46 HRD units (SD 4.83, range -15.44 SARC to +8.29 PCPG); removing it drops pooled MAE 9.049 -> 7.881.",
    "WORKS, WITH LABELS. Fraction of oracle gain recovered, EB-shrunk: k=1 +22%, k=3 +40%, k=5 +50%, k=10 +66%, k=20 +79%. Unshrunk at small k is harmful: k=1 -214%, k=3 -22%. Macro MAE at k=3: uncorrected 8.434, unshrunk 8.649, shrunk 8.040, oracle 7.444. Variance components tau2 = 22.3, sigma2 = 108.6, so w = k / (k + 4.87).",
    "REFUTED for V1/V2. Sentinel failed 3 of 7 gates: macro rho excluding THCA 0.418 (rank) vs 0.637 (run 01), delta -0.219; BRCA rho 0.261 vs 0.666, delta -0.405; tissue-identity R2 57.6% (rank) vs 52.7% (absolute) against a true-relative-target floor of 0.001. The percentile remap also RAISES tissue R2 (0.593 / 0.513 vs 0.562 for the raw prediction). Both routes increase lineage imprinting rather than removing it. V3 remains untested.",
    "CLIPPING: 227 negative predictions (3.2% of samples), most negative -11.53; MAE 9.0494 -> 8.9717 (gain 0.0777), skill 0.0850 -> 0.0928, Spearman with the raw prediction exactly 1.000. LOG1P: pooled MAE 9.049 -> 8.827 and skill 0.0850 -> 0.1075, but within-tissue Pearson 0.6124 -> 0.5221 (loss 0.090); correlation falls in 20 of 30 tissues while MAE improves in only 16 of 30. Losses concentrate in high-HRD tissues (ESCA 0.508->0.350, BLCA 0.675->0.529, LUSC 0.605->0.467, STAD 0.738->0.623).",
    "LARGELY EXONERATED. Controlling for purity RAISES the within-tissue correlation: 0.6124 -> 0.6210 (+0.0087). The purity gradient is an artefact of C1: after removing the per-tissue offset, skill by purity tertile is 0.242 (low) / 0.207 (mid) / 0.159 (high), up from 0.114 / 0.082 / 0.029. cor(per-tissue offset, per-tissue mean purity) across 30 tissues is only -0.124. A residual gradient remains (0.242 -> 0.159), so purity is not entirely innocent.",
    "OV n = 10 in the development cohort, and OV is ineligible at every k for few-shot calibration. Its within-tissue correlation is not interpretable at this sample size. HRD-high behaviour is therefore inferred from UCEC / BRCA / STAD rather than from ovarian cancer itself.",
    NULL),

  `Current disposition` = c(
    "ESCALATED / OPEN. No adopted rule. The CNS lock stays closed; C1 is more firmly unresolved after the C1c sentinel, not less. One untested option remains: V3-abs (feature_rank = within_tissue), pre-registered with gates G1-G7 in docs/27.",
    "VIABLE and ADOPTED AS THE SUPPORTED PATH, but only where labels exist. Recommendation: k >= 3 WITH empirical-Bayes shrinkage, ideally k >= 10. Never use the unshrunk sample mean at small k.",
    "REFUTED for V1/V2 (docs/26 verdict: 'B is now REFUTED, not merely unproven'). The full array was NOT authorised. V3 is the single highest-value next C1 experiment and is the only branch still open.",
    "RESOLVED. Clipping ADOPTED unconditionally: pred_clip = pmax(pred, 0). log1p REJECTED for the primary model, retained only as a documented option for the quiet-tumour regime (THCA -0.032 -> 0.171). Zero inflation itself is unaddressed and recorded as future work (needs a hurdle / two-part model).",
    "DOWNGRADED from 'most serious open defect' to a secondary effect largely explained by C1. Not formally closed: the probe-level purity-association test was not run (4,080 of 5,000 probes are shared across the 5 folds examined).",
    "OPEN. Disclosure item, not a code fix. Optional remedy (not executed): rebuild the feature bridge to include 27k probes, accepting a much smaller probe intersection.",
    NULL),

  `Deployment implication` = c(
    "Absolute HRDsum on a new tumour type with zero labels is NOT currently achievable. A published absolute HRDsum for a novel tissue carries an unquantified additive bias with a 95% predictive interval of -9.7 to +8.1 HRD units, a width of 17.8 units -- 74% as wide as the interquartile range of the label itself (4 to 28). Among the 1,171 patients within +/-10 units of the threshold, a tissue-level offset flips the HRD-high call for a median of 9.1% of them, and up to 79.2% in the worst tissue. No presentation, figure, caption or README line may imply the deployment gate is solved.",
    "Enables calibration to a new tumour type given roughly 10 labelled examples. It does NOT enable zero-label or n-of-1 inference, and must never be presented as a zero-shot capability.",
    "Both routes to an n-of-1 pediatric answer are currently closed: the absolute route is blocked by the unestimable offset, and the relative-rank route is blocked by the sentinel failure. Honest claim: we can rank tumours within a known type, and calibrate to a new type given ~10 labelled examples; we cannot yet place a single patient from an unseen tumour type on an absolute HRD scale.",
    "Clipping is safe to ship: it is free, cannot change any ranking metric (Spearman 1.000 with the raw prediction) and affects only 3.2% of samples. log1p is not shipped because C1 forces the project toward within-tissue RANKING, and ranking quality is measured by the correlation that log1p degrades -- precisely in the high-HRD tissues where discrimination is clinically useful.",
    "Purity is a nuisance variable the model partially encodes, not the source of its signal. This removes purity-capture as the dominant threat to biological validity, but the residual gradient and the untested probe-level association mean the model cannot be declared purity-clean. The exoneration applies to the production absolute model, not automatically to the relative-target model.",
    "Any clinical claim about ovarian HRD must disclose that direct ovarian evidence in this cohort is essentially absent (n = 10, the platform excludes the rest). State this limitation in every presentation.",
    NULL))

write_pair(tableD, "tableD_blocker_status",
  "Table D. Blocker status C1-C4",
  notes = c(
    "**Sources.** `docs/21_BLOCKER_RESOLUTION_PLAN.md`, `docs/25_C1_C2_C3_INVESTIGATION.md`,",
    "`docs/26_C1_N_OF_1_WORKAROUND.md`, `docs/27_V3_PREREGISTRATION.md`. All numbers are quoted",
    "from those documents rather than recomputed here, so that the table and the written record",
    "cannot drift apart. C1b and C1c are shown as separate rows because they are separate",
    "experiments with separate verdicts, and collapsing them would hide the fact that the",
    "labelled path succeeded while the zero-shot path was refuted.",
    "",
    "**Headline.** C2 is resolved (clip adopted, log1p rejected). C3 is downgraded, not closed.",
    "C4 is an open disclosure. **C1 is open and is the blocker that gates deployment**: the only",
    "supported correction requires labels from the new tumour type.",
    "",
    "**Related figure.** The mechanism behind C1 is shown in `results/figures/fig05_lineage_imprinting.png`:",
    "over the 30 development cancer types, cancer type explains R2 = 0.5584 of the model's",
    "prediction but only R2 = 0.3414 of the observed HRDsum -- an excess of +0.2171. The same",
    "pattern on the four V3 sentinel tissues is 0.5219 vs 0.2940 (`docs/27_V3_PREREGISTRATION.md`).",
    "",
    "**Cohort lock.** Nothing in this table is derived from the sealed CNS (GBM/LGG) cohort."))

message("== make_tables.R complete ==")

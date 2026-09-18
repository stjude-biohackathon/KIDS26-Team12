# =============================================================================
# scripts/analyze_cns.R -- one-shot locked CNS evaluation (docs/28 section 6)
# -----------------------------------------------------------------------------
# USAGE
#   Rscript scripts/analyze_cns.R \
#     <cns_predictions.tsv> <master_samples.tsv> <source_loco_predictions.tsv> \
#     <out_dir> [label]
#
#   e.g. Rscript scripts/analyze_cns.R \
#          results/locked_cns/cns_predictions.tsv \
#          data/processed/master_samples.tsv \
#          results/loco_run01/loco_predictions.tsv \
#          results/locked_cns primary_run01
#
#   Run from the repository root (source("R/model.R") is a relative path).
#
# WHAT THIS IS
#   The analysis half of the frozen CNS test. scripts/predict_frozen.R produces
#   the predictions under the locked-partition guard and the ledger; this script
#   scores them. It FITS NOTHING, RECALIBRATES NOTHING and TRANSFORMS NOTHING:
#   predictions arrive already back-transformed and C2-clipped, and that is
#   asserted rather than assumed (see "C2 GATE" below).
#
#   The metric list is frozen by docs/28 section 6, written 2026-09-18 11:58 CDT
#   before any CNS label was visible. This script implements exactly that list -
#   no more, no fewer - and was itself written before the cohort was unlocked.
#   Few-shot EB calibration (docs/28 section 6, "Few-shot calibration on CNS")
#   is a separate, separately-labelled script and is deliberately NOT here.
#
# INPUTS (all required; missing or malformed inputs abort before any output)
#   1  cns_predictions.tsv          output of scripts/predict_frozen.R on the
#                                   locked partition. Must carry sample_id,
#                                   predicted_reference_HRDsum (authoritative,
#                                   clipped), predicted_reference_HRDsum_raw,
#                                   clipping_rule, ood_score, ood, reportable.
#   2  master_samples.tsv           sample table; supplies cancer_type, partition,
#                                   HRDsum (the locked labels) and purity.
#   3  source_loco_predictions.tsv  results/loco_run01/loco_predictions.tsv --
#                                   the 30 development tissues, held out one at a
#                                   time, that CNS is placed within.
#   4  out_dir                      created if absent.
#   5  label                        optional free-text run label, recorded in the
#                                   outputs (e.g. "primary_run01" / "secondary_v3").
#
# OUTPUTS (written to out_dir)
#   cns_metrics_by_group.tsv        GBM, LGG (primary, separate) then CNS_pooled
#                                   (secondary): ranking, absolute, calibration
#                                   and ORACLE-null metrics, one row per group.
#   cns_placement_in_loco.tsv       THE KEY FRAMING. Each of Pearson, Spearman,
#                                   MAE, bias, calibration slope for GBM and LGG,
#                                   with its empirical percentile and rank inside
#                                   the 30 source LOCO tissues (rank out of 32).
#   cns_ood_summary.tsv             ood_score distribution and metrics computed
#                                   separately on reportable vs flagged subsets.
#   cns_threshold_exploratory.tsv   EXPLORATORY ONLY. Rank-based AUC at the
#                                   exploratory HRDsum >= 42 cut, n near the cut,
#                                   and predicted-vs-true crossings.
#   cns_analysis_summary.md         human-readable report with the numbers in it.
#   sessionInfo.txt                 seed, inputs, checksums, sessionInfo().
#
# CONVENTIONS AND GUARD RAILS
#   * C2 GATE. Every row must carry a clipping_rule beginning "C2_ADOPTED:" and
#     a non-negative authoritative prediction. A predictions file produced
#     without the C2 back-transform/clip is refused, so it cannot be analysed by
#     accident. This script never re-clips and never rescales.
#   * COHORT GATE. Every row must have partition == "locked_CNS" and cancer_type
#     in {GBM, LGG}. Anything else is a fatal error, not a filter.
#   * ORACLE LABELLING. The tissue-mean null is computed WITHIN GBM and WITHIN
#     LGG. It uses CNS labels and is unavailable at deployment for a genuinely
#     new lineage. Every output that reports it says so in a column.
#   * NO PASS/FAIL LINE. docs/28 section 6 draws none. The interpretive frame is
#     placement inside the source LOCO distribution.
#   * Seeds are fixed here and recorded in sessionInfo.txt.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

## ---- frozen constants (docs/28 section 6) -----------------------------------

# Bootstrap for the ranking CIs. docs/28 fixes 10,000 resamples and requires the
# seed to be recorded; 3391 is that seed and is written into every output.
BOOT_R          <- 10000L
BOOT_SEED       <- 3391L
BOOT_COVERAGE   <- 0.95

# Exploratory threshold. NOT a validated cutoff (see the section 5 block below).
EXPLORATORY_THRESHOLD <- 42
NEAR_THRESHOLD_WINDOW <- 10

CNS_TYPES     <- c("GBM", "LGG")
LOCKED_LABEL  <- "locked_CNS"
PRED_COL      <- "predicted_reference_HRDsum"
RAW_COL       <- "predicted_reference_HRDsum_raw"

# --- operationalising docs/28 section 7 --------------------------------------
# Section 7 pre-declares THREE outcomes ("rank transfers / calibration transfers
# / neither") in prose but does not attach numbers to them. The numeric readings
# below are declared HERE, in the script written before unlock, so the verdict
# printed at the end is computed from the data rather than chosen after seeing
# it. They are deliberately conservative and empirically anchored:
#
#   RANK TRANSFERS        if, in BOTH GBM and LGG, the Spearman point estimate
#                         is > 0 and the lower bound of its bootstrap 95% CI is
#                         > 0. (Rank is the primary endpoint; a CI that excludes
#                         zero in each lineage separately is the minimum.)
#   CALIBRATION TRANSFERS if, in BOTH GBM and LGG, (a) the model beats the
#                         ORACLE within-tissue mean null (skill > 0), (b) |bias|
#                         is no worse than the 90th percentile of |bias| across
#                         the 30 source LOCO tissues, and (c) the calibration
#                         slope lies inside the observed range of source-tissue
#                         calibration slopes. (b) and (c) are anchored to the
#                         source distribution rather than to invented constants,
#                         which is the same frame as section 6's placement test.
#
# Anything else prints as "neither transfers". The three statements of claim are
# quoted verbatim from docs/28 section 7 in the console summary and the .md.
BIAS_SOURCE_QUANTILE <- 0.90

## ---- CLI --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L || length(args) > 5L) {
  stop("Expected <cns_predictions.tsv> <master_samples.tsv> ",
       "<source_loco_predictions.tsv> <out_dir> [label]", call. = FALSE)
}
P_CNS    <- args[1]
P_MASTER <- args[2]
P_SRC    <- args[3]
OUT_DIR  <- args[4]
RUN_LABEL <- if (length(args) == 5L) args[5] else "unlabelled"

if (!file.exists("R/model.R")) {
  stop("Run this script from the repository root (R/model.R not found in '",
       normalizePath(".", mustWork = FALSE), "')", call. = FALSE)
}
# metrics(), tissue_mean_null(), null_panel(), within_tissue_metrics() and
# tissue_identity_r2() all live here and are REUSED rather than reimplemented,
# so the CNS numbers are produced by the same code that produced every source
# number in the project.
source("R/model.R")

read_required <- function(path, what) {
  if (!file.exists(path)) {
    stop("REQUIRED INPUT MISSING (", what, "): '", path, "'\n  (cwd = ",
         normalizePath(".", mustWork = FALSE), ")", call. = FALSE)
  }
  dt <- data.table::fread(path, data.table = TRUE, check.names = FALSE)
  if (nrow(dt) == 0L) stop("REQUIRED INPUT EMPTY (", what, "): '", path, "'", call. = FALSE)
  message(sprintf("  read %-56s %6d x %3d", path, nrow(dt), ncol(dt)))
  dt[]
}

need_cols <- function(dt, cols, what) {
  miss <- setdiff(cols, names(dt))
  if (length(miss)) {
    stop(what, " is missing required column(s): ", paste(miss, collapse = ", "),
         "\n  present: ", paste(names(dt), collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

fmt <- function(x, d = 3) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = d), "NA")
}

#' Ordinal suffix, so "GBM MAE is 12th of 32" reads as English.
ordinal <- function(k) {
  if (!is.finite(k)) return("?")
  k <- as.integer(k)
  sfx <- if (k %% 100L %in% 11:13) "th" else
    switch(as.character(k %% 10L), "1" = "st", "2" = "nd", "3" = "rd", "th")
  paste0(k, sfx)
}

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

message("== analyze_cns.R ==  label = ", RUN_LABEL,
        "  |  bootstrap seed = ", BOOT_SEED, ", R = ", BOOT_R)
message("-- reading inputs")

cns_raw <- read_required(P_CNS,    "CNS predictions")
master  <- read_required(P_MASTER, "master sample table")
src_raw <- read_required(P_SRC,    "source LOCO predictions")

## =============================================================================
## 1. JOIN AND COHORT GATE
## =============================================================================
need_cols(cns_raw, c("sample_id", PRED_COL, RAW_COL, "clipping_rule",
                     "ood_score", "ood", "reportable"), "CNS predictions file")
need_cols(master,  c("sample_id", "cancer_type", "partition", "HRDsum"),
          "master sample table")

if (anyDuplicated(cns_raw$sample_id)) {
  stop("Duplicate sample_id in the CNS predictions file; refusing to score an ",
       "ambiguous cohort.", call. = FALSE)
}

keep <- c("sample_id", "cancer_type", "partition", "HRDsum",
          if ("purity" %in% names(master)) "purity")
m <- unique(master[, ..keep])
if (anyDuplicated(m$sample_id)) {
  stop("Duplicate sample_id in the master sample table with conflicting ",
       "cancer_type/partition/HRDsum; cannot join unambiguously.", call. = FALSE)
}

idx <- match(cns_raw$sample_id, m$sample_id)
if (anyNA(idx)) {
  stop(sprintf("%d of %d predicted samples are absent from the master table (%s). ",
               sum(is.na(idx)), nrow(cns_raw), basename(P_MASTER)),
       "Refusing to analyse samples whose partition and label cannot be established.",
       call. = FALSE)
}

d <- data.table::copy(cns_raw)
d[, cancer_type   := m$cancer_type[idx]]
d[, master_partition := m$partition[idx]]
d[, observed_HRDsum  := as.numeric(m$HRDsum[idx])]
d[, purity := if ("purity" %in% names(m)) as.numeric(m$purity[idx]) else NA_real_]

# Cohort gate. These are assertions, NOT filters: a row that does not belong in
# the locked CNS evaluation means the wrong file was passed, and silently
# dropping it would silently change the cohort.
bad_part <- d[master_partition != LOCKED_LABEL]
if (nrow(bad_part)) {
  stop(sprintf("COHORT GATE: %d of %d rows have partition != '%s' (seen: %s). ",
               nrow(bad_part), nrow(d), LOCKED_LABEL,
               paste(sort(unique(bad_part$master_partition)), collapse = ", ")),
       "This script analyses the locked CNS partition only.", call. = FALSE)
}
if ("partition" %in% names(cns_raw)) {
  pp <- as.character(cns_raw$partition)
  disagree <- pp != "not_supplied" & pp != d$master_partition
  if (any(disagree)) {
    stop(sprintf(paste0("COHORT GATE: partition column of the predictions file disagrees ",
                        "with the master table for %d row(s). Two sources of truth for the ",
                        "lock must agree before either is trusted."), sum(disagree)), call. = FALSE)
  }
}
bad_type <- d[!cancer_type %chin% CNS_TYPES]
if (nrow(bad_type)) {
  stop(sprintf("COHORT GATE: %d of %d rows have cancer_type outside {%s} (seen: %s).",
               nrow(bad_type), nrow(d), paste(CNS_TYPES, collapse = ", "),
               paste(sort(unique(bad_type$cancer_type)), collapse = ", ")), call. = FALSE)
}
if (!all(CNS_TYPES %chin% unique(d$cancer_type))) {
  stop("COHORT GATE: expected both GBM and LGG in the locked cohort; found only: ",
       paste(sort(unique(d$cancer_type)), collapse = ", "), call. = FALSE)
}
if (!any(is.finite(d$observed_HRDsum))) {
  stop("No finite HRDsum values joined for the CNS cohort - check the master table.",
       call. = FALSE)
}
n_no_label <- sum(!is.finite(d$observed_HRDsum))
if (n_no_label > 0L) {
  message(sprintf(paste0("  NOTE: %d CNS row(s) have no finite HRDsum; metrics() drops ",
                         "non-finite pairs and reports the surviving n."), n_no_label))
}

## =============================================================================
## 2. C2 GATE -- predictions arrive clipped; this script applies NO transform
## =============================================================================
# docs/28 section 6: "no post-hoc clipping changes". predict_frozen.R applies the
# bundle's target_transform inverse and then pmax(., 0), writing the rule into
# every row. We verify that happened rather than repeating it: re-clipping here
# would mask a predictions file that had never been clipped at all.
cr <- unique(as.character(d$clipping_rule))
if (any(is.na(cr)) || !all(nzchar(cr))) {
  stop("C2 GATE: clipping_rule is blank/NA for at least one row.", call. = FALSE)
}
if (!all(grepl("^C2_ADOPTED:", cr))) {
  stop("C2 GATE: clipping_rule does not record the adopted C2 rule (found: ",
       paste(cr, collapse = " | "), "). These predictions were not produced by a ",
       "C2-compliant inference run and must not be analysed.", call. = FALSE)
}
if (length(cr) > 1L) {
  stop("C2 GATE: more than one clipping_rule in a single predictions file (",
       paste(cr, collapse = " | "), "); the cohort was not scored uniformly.",
       call. = FALSE)
}
neg <- sum(is.finite(d[[PRED_COL]]) & d[[PRED_COL]] < 0)
if (neg > 0L) {
  stop(sprintf(paste0("C2 GATE: %d authoritative prediction(s) are negative despite ",
                      "clipping_rule = '%s'. HRDsum is bounded below at zero."), neg, cr),
       call. = FALSE)
}
CLIP_RULE <- cr
n_clipped <- sum(is.finite(d[[RAW_COL]]) & d[[RAW_COL]] < 0)
message(sprintf("  C2 gate passed: rule = %s;  %d of %d raw predictions were negative pre-clip (%.2f%%)",
                CLIP_RULE, n_clipped, nrow(d), 100 * n_clipped / nrow(d)))

d[, predicted := as.numeric(get(PRED_COL))]   # authoritative, already clipped
d[, residual  := predicted - observed_HRDsum] # signed: positive = over-prediction

## =============================================================================
## 3. METRIC PANEL -- GBM and LGG first (primary), pooled CNS second
## =============================================================================

#' Percentile bootstrap CIs for Pearson and Spearman from ONE set of resamples.
#' Pairs are resampled with replacement; the seed is set once per call and is
#' recorded in every output row.
boot_rank_ci <- function(y, p, R = BOOT_R, seed = BOOT_SEED, coverage = BOOT_COVERAGE) {
  ok <- is.finite(y) & is.finite(p); y <- y[ok]; p <- p[ok]; n <- length(y)
  out <- list(pearson_lo = NA_real_, pearson_hi = NA_real_,
              spearman_lo = NA_real_, spearman_hi = NA_real_, boot_n_ok = 0L)
  if (n < 5L) return(out)
  set.seed(seed)
  pr <- numeric(R); sr <- numeric(R)
  for (b in seq_len(R)) {
    i <- sample.int(n, n, replace = TRUE)
    yb <- y[i]; pb <- p[i]
    if (sd(yb) > 0 && sd(pb) > 0) {
      pr[b] <- cor(yb, pb)
      sr[b] <- cor(yb, pb, method = "spearman")
    } else {
      pr[b] <- NA_real_; sr[b] <- NA_real_
    }
  }
  a <- (1 - coverage) / 2
  pq <- quantile(pr, c(a, 1 - a), na.rm = TRUE, names = FALSE)
  sq <- quantile(sr, c(a, 1 - a), na.rm = TRUE, names = FALSE)
  list(pearson_lo = pq[1], pearson_hi = pq[2],
       spearman_lo = sq[1], spearman_hi = sq[2],
       boot_n_ok = sum(is.finite(pr)))
}

#' One row of the frozen metric panel for one group.
#' `cancer` is passed through so the tissue-mean null is computed WITHIN GBM and
#' WITHIN LGG even for the pooled row (docs/28 section 6).
group_panel <- function(y, p, cancer, group_label, scope) {
  ok <- is.finite(y) & is.finite(p)
  y <- y[ok]; p <- p[ok]; cancer <- as.character(cancer)[ok]; n <- length(y)
  mm <- metrics(y, p)                     # R/model.R -- shared with every source result
  np <- null_panel(y, p, cancer)          # R/model.R -- ORACLE within-tissue mean null
  ci <- boot_rank_ci(y, p)
  data.frame(
    group = group_label, scope = scope, label = RUN_LABEL, n = n,
    # --- ranking (primary endpoint) -----------------------------------------
    Pearson = mm$Pearson, Pearson_lo95 = ci$pearson_lo, Pearson_hi95 = ci$pearson_hi,
    Spearman = mm$Spearman, Spearman_lo95 = ci$spearman_lo, Spearman_hi95 = ci$spearman_hi,
    boot_R = BOOT_R, boot_seed = BOOT_SEED, boot_method = "percentile_pairs",
    # --- absolute ------------------------------------------------------------
    MAE = mm$MAE, RMSE = mm$RMSE,
    median_absolute_error = if (n) median(abs(y - p)) else NA_real_,
    # --- calibration ---------------------------------------------------------
    bias_mean_signed_error = mm$bias,
    calibration_intercept = mm$calibration_intercept,
    calibration_slope = mm$calibration_slope,
    observed_mean = if (n) mean(y) else NA_real_,
    observed_sd   = if (n > 1L) sd(y) else NA_real_,
    observed_min  = if (n) min(y) else NA_real_,
    observed_max  = if (n) max(y) else NA_real_,
    predicted_mean = if (n) mean(p) else NA_real_,
    predicted_sd   = if (n > 1L) sd(p) else NA_real_,
    predicted_min  = if (n) min(p) else NA_real_,
    predicted_max  = if (n) max(p) else NA_real_,
    range_compression_ratio = if (n > 1L && sd(y) > 0) sd(p) / sd(y) else NA_real_,
    # --- null (ORACLE) -------------------------------------------------------
    null_type = "within_CNS_group_tissue_mean",
    null_is_oracle = TRUE,
    null_oracle_note = "uses CNS labels; NOT available at deployment for a new lineage",
    MAE_oracle_tissue_mean_null = np$MAE_tissue_mean_null,
    skill_vs_oracle_tissue_mean = np$skill_vs_tissue_mean,
    R2_vs_group_mean = mm$R2,
    clipping_rule = CLIP_RULE,
    stringsAsFactors = FALSE)
}

message("-- metric panel (GBM, LGG primary; pooled CNS secondary)")
panels <- list()
for (g in CNS_TYPES) {
  s <- d[cancer_type == g]
  panels[[g]] <- group_panel(s$observed_HRDsum, s$predicted, s$cancer_type, g,
                             "primary_per_lineage")
}
panels[["CNS_pooled"]] <- group_panel(d$observed_HRDsum, d$predicted, d$cancer_type,
                                      "CNS_pooled", "secondary_pooled")
metrics_by_group <- data.table::rbindlist(panels)

# Pooled-only diagnostics that are meaningless within a single lineage: how much
# of the pooled CNS prediction (and of the pooled CNS label) is explained by
# GBM-vs-LGG identity alone, and the within-lineage-centred agreement.
wt <- within_tissue_metrics(d$observed_HRDsum, d$predicted, d$cancer_type)
R2_TISSUE_PRED <- tissue_identity_r2(d$predicted, d$cancer_type)
R2_TISSUE_OBS  <- tissue_identity_r2(d$observed_HRDsum, d$cancer_type)
metrics_by_group[, within_lineage_Pearson  := c(NA_real_, NA_real_, wt$within_Pearson)]
metrics_by_group[, within_lineage_Spearman := c(NA_real_, NA_real_, wt$within_Spearman)]
metrics_by_group[, lineage_R2_of_prediction := c(NA_real_, NA_real_, R2_TISSUE_PRED)]
metrics_by_group[, lineage_R2_of_observed   := c(NA_real_, NA_real_, R2_TISSUE_OBS)]

## =============================================================================
## 4. THE KEY FRAMING -- place GBM and LGG inside the 30 source LOCO tissues
## =============================================================================
# docs/28 section 6: "place GBM and LGG inside the empirical distribution of the
# 30 source LOCO tissues, and ask 'would CNS look unusual if it had simply been
# the 31st held-out tissue?'" No pass/fail line is drawn here, by design.
#
# The source per-tissue values are RECOMPUTED from the source predictions under
# the same C2 clipping, using the same metrics() function used for CNS above, so
# the two sides of the comparison are produced by identical code.
need_cols(src_raw, c("cancer_type", "actual", PRED_COL), "source LOCO predictions")
src <- data.table::copy(src_raw)
hit <- intersect(CNS_TYPES, unique(as.character(src$cancer_type)))
if (length(hit)) {
  stop("SOURCE GATE: the source LOCO predictions contain locked type(s) ", 
       paste(hit, collapse = ", "), ". The source distribution must be the 30 ",
       "development tissues only, or the placement comparison is circular.",
       call. = FALSE)
}
src[, pred_clip := pmax(as.numeric(get(PRED_COL)), 0)]
src[, actual := as.numeric(actual)]

src_tissue <- data.table::rbindlist(lapply(split(src, src$cancer_type), function(s) {
  mm <- metrics(s$actual, s$pred_clip)
  data.frame(cancer_type = s$cancer_type[1], n = mm$n, Pearson = mm$Pearson,
             Spearman = mm$Spearman, MAE = mm$MAE,
             bias = mm$bias, calibration_slope = mm$calibration_slope,
             stringsAsFactors = FALSE)
}))
N_SRC <- nrow(src_tissue)
message(sprintf("  source LOCO distribution: %d held-out tissues, %d samples",
                N_SRC, nrow(src)))
if (N_SRC < 10L) {
  stop("SOURCE GATE: only ", N_SRC, " source tissues found; the placement ",
       "percentile is not interpretable on so few. Expected 30.", call. = FALSE)
}

# Placement statistics and their orientation.
#   higher_better : Pearson, Spearman
#   lower_better  : MAE, |bias|, |calibration slope - 1|
# For "lower_better" statistics the placement uses the DISTANCE measure named in
# `ranked_on` (e.g. |bias|), because a tissue with bias -8 is no better
# calibrated than one with +8; the signed value is reported alongside it.
PLACEMENT <- list(
  list(stat = "Pearson",           col = "Pearson",           dir = "higher_better", xf = function(v) v),
  list(stat = "Spearman",          col = "Spearman",          dir = "higher_better", xf = function(v) v),
  list(stat = "MAE",               col = "MAE",               dir = "lower_better",  xf = function(v) v),
  list(stat = "bias",              col = "bias",              dir = "lower_better",  xf = function(v) abs(v)),
  list(stat = "calibration_slope", col = "calibration_slope", dir = "lower_better",  xf = function(v) abs(v - 1))
)

#' Empirical percentile of a CNS value within the source tissues.
#' Defined so that HIGH IS ALWAYS GOOD, whichever way the statistic points:
#'   higher_better : pct = 100 * mean(source <= value)  (share of source tissues
#'                   the CNS group matches or beats)
#'   lower_better  : pct = 100 * mean(source >= value)
#' Rank is over the COMBINED set of 30 source tissues + GBM + LGG = 32, ordered
#' best first, ties given the minimum rank, so "GBM MAE is 12th of 32" reads
#' directly.
placement_rows <- list()
for (spec in PLACEMENT) {
  src_signed <- src_tissue[[spec$col]]
  src_score  <- spec$xf(src_signed)
  cns_signed <- vapply(CNS_TYPES, function(g)
    as.numeric(metrics_by_group[group == g][[switch(spec$col,
      Pearson = "Pearson", Spearman = "Spearman", MAE = "MAE",
      bias = "bias_mean_signed_error", calibration_slope = "calibration_slope")]]),
    numeric(1))
  cns_score <- spec$xf(cns_signed)

  combined <- c(src_score, cns_score)
  names(combined) <- c(as.character(src_tissue$cancer_type), CNS_TYPES)
  ord_rank <- if (spec$dir == "higher_better") rank(-combined, ties.method = "min")
              else rank(combined, ties.method = "min")

  for (g in CNS_TYPES) {
    v_signed <- cns_signed[[g]]
    v_score  <- cns_score[[g]]
    pct <- if (!is.finite(v_score)) NA_real_ else
      100 * mean(if (spec$dir == "higher_better") src_score <= v_score
                 else src_score >= v_score, na.rm = TRUE)
    rk <- unname(ord_rank[g])
    placement_rows[[length(placement_rows) + 1L]] <- data.frame(
      statistic = spec$stat, group = g, label = RUN_LABEL,
      value = v_signed, ranked_on = if (identical(spec$stat, "bias")) "abs(bias)"
                                    else if (identical(spec$stat, "calibration_slope")) "abs(slope - 1)"
                                    else spec$stat,
      ranked_value = v_score, direction = spec$dir,
      n_source_tissues = N_SRC,
      source_min = min(src_score, na.rm = TRUE),
      source_q25 = as.numeric(quantile(src_score, 0.25, na.rm = TRUE, names = FALSE)),
      source_median = median(src_score, na.rm = TRUE),
      source_q75 = as.numeric(quantile(src_score, 0.75, na.rm = TRUE, names = FALSE)),
      source_max = max(src_score, na.rm = TRUE),
      source_mean = mean(src_score, na.rm = TRUE),
      source_sd = sd(src_score, na.rm = TRUE),
      percentile_in_source = pct,
      rank_of_combined = rk, n_combined = length(combined),
      rank_statement = sprintf("%s %s is %s of %d (30 source LOCO tissues + GBM + LGG)",
                               g, spec$stat, ordinal(rk), length(combined)),
      inside_source_range = is.finite(v_score) &&
        v_score >= min(src_score, na.rm = TRUE) && v_score <= max(src_score, na.rm = TRUE),
      stringsAsFactors = FALSE)
  }
}
placement <- data.table::rbindlist(placement_rows)

## =============================================================================
## 5. OOD -- frozen flags, and does accuracy differ between the two groups?
## =============================================================================
d[, ood_flag := as.logical(ood)]
d[, reportable_flag := as.logical(reportable)]

ood_subset_row <- function(s, group_label, subset_label) {
  y <- s$observed_HRDsum; p <- s$predicted
  ok <- is.finite(y) & is.finite(p)
  mm <- if (sum(ok) >= 3L) metrics(y[ok], p[ok]) else
    data.frame(n = sum(ok), MAE = NA_real_, RMSE = NA_real_, R2 = NA_real_,
               Pearson = NA_real_, Spearman = NA_real_,
               calibration_intercept = NA_real_, calibration_slope = NA_real_,
               bias = NA_real_)
  os <- s$ood_score[is.finite(s$ood_score)]
  q <- if (length(os)) as.numeric(quantile(os, c(0, .25, .5, .75, 1), names = FALSE)) else rep(NA_real_, 5)
  data.frame(group = group_label, subset = subset_label, label = RUN_LABEL,
             n = nrow(s), n_with_label = sum(ok),
             n_reportable = sum(s$reportable_flag, na.rm = TRUE),
             n_flagged = sum(!s$reportable_flag, na.rm = TRUE),
             n_ood_flagged = sum(s$ood_flag, na.rm = TRUE),
             n_qc_fail = if ("qc_fail" %in% names(s)) sum(as.logical(s$qc_fail), na.rm = TRUE) else NA_integer_,
             frac_flagged = if (nrow(s)) mean(!s$reportable_flag, na.rm = TRUE) else NA_real_,
             ood_min = q[1], ood_q25 = q[2], ood_median = q[3], ood_q75 = q[4], ood_max = q[5],
             Pearson = mm$Pearson, Spearman = mm$Spearman, MAE = mm$MAE,
             RMSE = mm$RMSE, bias = mm$bias,
             calibration_slope = mm$calibration_slope,
             stringsAsFactors = FALSE)
}

message("-- OOD summary (reportable vs flagged)")
ood_rows <- list()
for (g in c(CNS_TYPES, "CNS_pooled")) {
  s <- if (g == "CNS_pooled") d else d[cancer_type == g]
  ood_rows[[length(ood_rows) + 1L]] <- ood_subset_row(s, g, "all")
  ood_rows[[length(ood_rows) + 1L]] <- ood_subset_row(s[reportable_flag %in% TRUE], g, "reportable")
  ood_rows[[length(ood_rows) + 1L]] <- ood_subset_row(s[!(reportable_flag %in% TRUE)], g, "flagged")
}
ood_summary <- data.table::rbindlist(ood_rows)

## ---- purity: secondary residual diagnostic (skipped gracefully if absent) ----
purity_note <- NULL
if (!("purity" %in% names(d)) || !any(is.finite(d$purity))) {
  purity_note <- "purity is not populated for the CNS rows; the residual-vs-purity diagnostic was SKIPPED (docs/28 section 6 makes it conditional)."
  message("  NOTE: ", purity_note)
  metrics_by_group[, purity_available := FALSE]
  metrics_by_group[, residual_vs_purity_Spearman := NA_real_]
  metrics_by_group[, n_with_purity := 0L]
} else {
  pr <- vapply(c(CNS_TYPES, "CNS_pooled"), function(g) {
    s <- if (g == "CNS_pooled") d else d[cancer_type == g]
    ok <- is.finite(s$residual) & is.finite(s$purity)
    if (sum(ok) < 5L || sd(s$purity[ok]) == 0) return(NA_real_)
    suppressWarnings(cor(s$residual[ok], s$purity[ok], method = "spearman"))
  }, numeric(1))
  np <- vapply(c(CNS_TYPES, "CNS_pooled"), function(g) {
    s <- if (g == "CNS_pooled") d else d[cancer_type == g]
    as.numeric(sum(is.finite(s$residual) & is.finite(s$purity)))
  }, numeric(1))
  metrics_by_group[, purity_available := TRUE]
  metrics_by_group[, residual_vs_purity_Spearman := as.numeric(pr)]
  metrics_by_group[, n_with_purity := as.integer(np)]
  purity_note <- sprintf("purity populated for %d of %d CNS rows; Spearman(residual, purity) reported per group.",
                         sum(is.finite(d$purity)), nrow(d))
  message("  ", purity_note)
}

## =============================================================================
## 6. EXPLORATORY THRESHOLD SECTION -- NOT A VALIDATED CUTOFF
## =============================================================================
# docs/28 section 6: "Threshold analysis is exploratory only. HRDsum >= 42 may be
# shown, labelled 'ranking performance at an exploratory cutoff'. It is not a
# primary endpoint and no continuous prediction will be presented as a
# probability."
#
# WHAT THE NUMBERS BELOW ARE AND ARE NOT:
#   * The AUC is the rank-based (Mann-Whitney) AUC of the CONTINUOUS predicted
#     HRDsum against the binary label (observed HRDsum >= 42). It measures
#     ranking only. It is NOT a probability, NOT calibrated, and involves no
#     fitted threshold on the prediction side.
#   * The "crossings" table applies the SAME numeric cut to the prediction. That
#     is a convenience display, not a validated decision rule: the cut of 42
#     describes the TCGA reference distribution (config/analysis_protocol.json)
#     and is not a clinical threshold, least of all a pediatric one.
#   * n near the threshold quantifies how many samples sit within +/- 10 of the
#     cut, i.e. how many decisions this display would get right or wrong by
#     accident given the model's absolute error.
auc_rank <- function(score, positive) {
  ok <- is.finite(score) & !is.na(positive)
  score <- score[ok]; positive <- as.logical(positive[ok])
  n1 <- sum(positive); n0 <- sum(!positive)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(score)                       # ties get the mid-rank: the standard form
  (sum(r[positive]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

thr_rows <- list()
for (g in c(CNS_TYPES, "CNS_pooled")) {
  s <- if (g == "CNS_pooled") d else d[cancer_type == g]
  ok <- is.finite(s$observed_HRDsum) & is.finite(s$predicted)
  y <- s$observed_HRDsum[ok]; p <- s$predicted[ok]
  pos <- y >= EXPLORATORY_THRESHOLD
  phat <- p >= EXPLORATORY_THRESHOLD
  tp <- sum(pos & phat); fp <- sum(!pos & phat)
  fn <- sum(pos & !phat); tn <- sum(!pos & !phat)
  thr_rows[[length(thr_rows) + 1L]] <- data.frame(
    group = g, label = RUN_LABEL, status = "EXPLORATORY_NOT_A_VALIDATED_CUTOFF",
    threshold = EXPLORATORY_THRESHOLD, n = length(y),
    n_observed_above = sum(pos), n_observed_below = sum(!pos),
    n_predicted_above = sum(phat),
    AUC_rank_based_continuous_score = auc_rank(p, pos),
    score_is_probability = FALSE,
    n_observed_near_threshold = sum(abs(y - EXPLORATORY_THRESHOLD) <= NEAR_THRESHOLD_WINDOW),
    n_predicted_near_threshold = sum(abs(p - EXPLORATORY_THRESHOLD) <= NEAR_THRESHOLD_WINDOW),
    near_threshold_window = NEAR_THRESHOLD_WINDOW,
    crossings_TP = tp, crossings_FP = fp, crossings_FN = fn, crossings_TN = tn,
    sensitivity = if ((tp + fn) > 0) tp / (tp + fn) else NA_real_,
    specificity = if ((tn + fp) > 0) tn / (tn + fp) else NA_real_,
    ppv = if ((tp + fp) > 0) tp / (tp + fp) else NA_real_,
    npv = if ((tn + fn) > 0) tn / (tn + fn) else NA_real_,
    n_discordant_crossings = fp + fn,
    stringsAsFactors = FALSE)
}
threshold_tab <- data.table::rbindlist(thr_rows)

## =============================================================================
## 7. PRE-DECLARED VERDICT (docs/28 section 7), computed from the data
## =============================================================================
get_g <- function(g, col) as.numeric(metrics_by_group[group == g][[col]])

rank_ok <- vapply(CNS_TYPES, function(g)
  isTRUE(get_g(g, "Spearman") > 0) && isTRUE(get_g(g, "Spearman_lo95") > 0),
  logical(1))
RANK_TRANSFERS <- all(rank_ok)

bias_cut  <- as.numeric(quantile(abs(src_tissue$bias), BIAS_SOURCE_QUANTILE,
                                 na.rm = TRUE, names = FALSE))
slope_lo  <- min(src_tissue$calibration_slope, na.rm = TRUE)
slope_hi  <- max(src_tissue$calibration_slope, na.rm = TRUE)
cal_ok <- vapply(CNS_TYPES, function(g) {
  isTRUE(get_g(g, "skill_vs_oracle_tissue_mean") > 0) &&
    isTRUE(abs(get_g(g, "bias_mean_signed_error")) <= bias_cut) &&
    isTRUE(get_g(g, "calibration_slope") >= slope_lo) &&
    isTRUE(get_g(g, "calibration_slope") <= slope_hi)
}, logical(1))
CALIBRATION_TRANSFERS <- all(cal_ok)

CLAIMS <- c(
  rank_only = paste("methylation carries HRD-associated signal that transfers in rank",
                    "to an unseen lineage, but absolute genomic-scar calibration remains",
                    "lineage-dependent. Deployment gate stays shut."),
  both      = paste("evidence of ranking and quantitative transfer into adult CNS tumours,",
                    "warranting independent pediatric testing. Still not a clinical claim,",
                    "still not validated, still needs PBTP."),
  neither   = paste("the pan-cancer methylation-HRD association did not generalise to the",
                    "held-out CNS lineage under a strict frozen zero-shot test. This is a",
                    "real and publishable result."))

OUTCOME_KEY <- if (RANK_TRANSFERS && CALIBRATION_TRANSFERS) "both" else
               if (RANK_TRANSFERS) "rank_only" else "neither"
OUTCOME_NAME <- c(rank_only = "RANK TRANSFERS, CALIBRATION DOES NOT",
                  both = "BOTH RANK AND CALIBRATION TRANSFER",
                  neither = "NEITHER TRANSFERS")[[OUTCOME_KEY]]
OUTCOME_CLAIM <- CLAIMS[[OUTCOME_KEY]]

## =============================================================================
## 8. WRITE OUTPUTS
## =============================================================================
wtsv <- function(x, name) {
  p <- file.path(OUT_DIR, name)
  write.table(x, p, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
  message(sprintf("  wrote %-40s %5d rows", name, nrow(x)))
  invisible(p)
}
message("-- writing outputs to ", OUT_DIR)
wtsv(metrics_by_group, "cns_metrics_by_group.tsv")
wtsv(placement,        "cns_placement_in_loco.tsv")
wtsv(ood_summary,      "cns_ood_summary.tsv")
wtsv(threshold_tab,    "cns_threshold_exploratory.tsv")

## ---- human-readable summary -------------------------------------------------
pl <- function(st, g, field) {
  v <- placement[statistic == st & group == g][[field]]
  if (!length(v)) NA else v[1]
}
md_metric_line <- function(g) {
  sprintf("| %s | %d | %s [%s, %s] | %s [%s, %s] | %s | %s | %s | %s | %s | %s |",
          g, get_g(g, "n"),
          fmt(get_g(g, "Pearson")), fmt(get_g(g, "Pearson_lo95")), fmt(get_g(g, "Pearson_hi95")),
          fmt(get_g(g, "Spearman")), fmt(get_g(g, "Spearman_lo95")), fmt(get_g(g, "Spearman_hi95")),
          fmt(get_g(g, "MAE"), 2), fmt(get_g(g, "RMSE"), 2),
          fmt(get_g(g, "median_absolute_error"), 2),
          fmt(get_g(g, "bias_mean_signed_error"), 2),
          fmt(get_g(g, "calibration_slope")),
          fmt(get_g(g, "range_compression_ratio")))
}
md_place_line <- function(st, g) {
  sprintf("| %s | %s | %s | %s | %s | %s | %s%% | %s | %s |",
          st, g, fmt(pl(st, g, "value")), pl(st, g, "ranked_on"),
          fmt(pl(st, g, "source_median")),
          sprintf("[%s, %s]", fmt(pl(st, g, "source_min")), fmt(pl(st, g, "source_max"))),
          fmt(pl(st, g, "percentile_in_source"), 1),
          pl(st, g, "rank_statement"),
          if (isTRUE(pl(st, g, "inside_source_range"))) "yes" else "NO")
}
md_ood_line <- function(g, sub) {
  r <- ood_summary[group == g & subset == sub]
  if (!nrow(r)) return(NULL)
  sprintf("| %s | %s | %d | %s | %s | %s | %s |", g, sub, r$n[1],
          fmt(r$MAE[1], 2), fmt(r$Pearson[1]), fmt(r$Spearman[1]),
          sprintf("%s / %s / %s", fmt(r$ood_min[1], 2), fmt(r$ood_median[1], 2), fmt(r$ood_max[1], 2)))
}
md_thr_line <- function(g) {
  r <- threshold_tab[group == g]
  sprintf("| %s | %d | %d | %s | %d | %d | %d | %d | %d |", g, r$n[1],
          r$n_observed_above[1], fmt(r$AUC_rank_based_continuous_score[1]),
          r$n_observed_near_threshold[1], r$crossings_TP[1], r$crossings_FP[1],
          r$crossings_FN[1], r$crossings_TN[1])
}

md <- c(
  "# Locked CNS evaluation -- results",
  "",
  sprintf("Run label: **%s**.  Generated by `scripts/analyze_cns.R` on %s.",
          RUN_LABEL, format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "Protocol: `docs/28_CANDIDATE_SELECTION_AND_CNS_PROTOCOL.md` section 6, frozen",
  "2026-09-18 11:58 CDT before any CNS label was visible. The metric list below is",
  "exactly that list. Nothing here refits, recalibrates or re-clips: predictions",
  sprintf("arrive from `scripts/predict_frozen.R` under `%s`.", CLIP_RULE),
  "",
  "## Inputs",
  "",
  sprintf("* CNS predictions: `%s`", P_CNS),
  sprintf("* Master samples : `%s`", P_MASTER),
  sprintf("* Source LOCO    : `%s` (%d tissues, %d samples)", P_SRC, N_SRC, nrow(src)),
  sprintf("* Cohort         : %d samples, partition == `%s`; GBM n = %d, LGG n = %d",
          nrow(d), LOCKED_LABEL, sum(d$cancer_type == "GBM"), sum(d$cancer_type == "LGG")),
  sprintf("* Bootstrap      : %d resamples, seed %d, percentile method on resampled pairs",
          BOOT_R, BOOT_SEED),
  "",
  "## 1. Metrics -- GBM and LGG first (primary), pooled CNS second",
  "",
  "| Group | n | Pearson [95% CI] | Spearman [95% CI] | MAE | RMSE | medAE | bias | calib slope | sd(pred)/sd(obs) |",
  "|---|---|---|---|---|---|---|---|---|---|",
  md_metric_line("GBM"), md_metric_line("LGG"), md_metric_line("CNS_pooled"),
  "",
  "Observed vs predicted distribution:",
  "",
  "| Group | obs mean (SD) | obs range | pred mean (SD) | pred range | calib intercept |",
  "|---|---|---|---|---|---|",
  paste0(vapply(c(CNS_TYPES, "CNS_pooled"), function(g)
    sprintf("| %s | %s (%s) | %s - %s | %s (%s) | %s - %s | %s |", g,
            fmt(get_g(g, "observed_mean"), 2), fmt(get_g(g, "observed_sd"), 2),
            fmt(get_g(g, "observed_min"), 1), fmt(get_g(g, "observed_max"), 1),
            fmt(get_g(g, "predicted_mean"), 2), fmt(get_g(g, "predicted_sd"), 2),
            fmt(get_g(g, "predicted_min"), 1), fmt(get_g(g, "predicted_max"), 1),
            fmt(get_g(g, "calibration_intercept"), 2)), character(1)), collapse = "\n"),
  "",
  "### Null comparison -- ORACLE baseline, read the caveat",
  "",
  "The tissue-mean null is computed **within GBM and within LGG**. It uses the CNS",
  "labels themselves and is **not available at deployment** for a genuinely new",
  "lineage, so it is an *oracle* comparator. It answers \"does the model beat",
  "knowing the tissue's mean?\" -- the question that matters -- and it flatters",
  "nothing.",
  "",
  "| Group | model MAE | ORACLE within-group mean-null MAE | skill = 1 - MAE/MAE_null |",
  "|---|---|---|---|",
  paste0(vapply(c(CNS_TYPES, "CNS_pooled"), function(g)
    sprintf("| %s | %s | %s | %s |", g, fmt(get_g(g, "MAE"), 2),
            fmt(get_g(g, "MAE_oracle_tissue_mean_null"), 2),
            fmt(get_g(g, "skill_vs_oracle_tissue_mean"))), character(1)), collapse = "\n"),
  "",
  sprintf("Pooled CNS only: lineage (GBM vs LGG) R2 of the prediction = %s, of the observed HRDsum = %s; within-lineage Pearson = %s, Spearman = %s.",
          fmt(R2_TISSUE_PRED), fmt(R2_TISSUE_OBS),
          fmt(wt$within_Pearson), fmt(wt$within_Spearman)),
  "",
  sprintf("Purity: %s", purity_note),
  "",
  "## 2. Placement in the source LOCO distribution (the primary frame)",
  "",
  sprintf("Would CNS look unusual as the 31st held-out tissue? Each value is placed in the empirical distribution of the %d source tissues, recomputed from `%s` under the same C2 clipping and the same `metrics()` function.",
          N_SRC, basename(P_SRC)),
  "",
  "Percentiles are oriented so **higher is always better**: for MAE, |bias| and",
  "|slope - 1| the percentile is the share of source tissues the CNS group is no",
  "worse than. Ranks are over 30 source tissues + GBM + LGG = 32, best first.",
  "",
  "| Statistic | Group | CNS value | ranked on | source median | source range | percentile | rank | inside source range |",
  "|---|---|---|---|---|---|---|---|---|",
  unlist(lapply(c("Pearson", "Spearman", "MAE", "bias", "calibration_slope"),
                function(st) vapply(CNS_TYPES, function(g) md_place_line(st, g), character(1)))),
  "",
  "## 3. OOD (frozen flags)",
  "",
  "| Group | subset | n | MAE | Pearson | Spearman | ood_score min/median/max |",
  "|---|---|---|---|---|---|---|",
  unlist(lapply(c(CNS_TYPES, "CNS_pooled"),
                function(g) unlist(lapply(c("all", "reportable", "flagged"),
                                          function(s) md_ood_line(g, s))))),
  "",
  sprintf("Pooled: %d reportable, %d flagged (%s%% of the cohort).",
          ood_summary[group == "CNS_pooled" & subset == "all"]$n_reportable[1],
          ood_summary[group == "CNS_pooled" & subset == "all"]$n_flagged[1],
          fmt(100 * ood_summary[group == "CNS_pooled" & subset == "all"]$frac_flagged[1], 1)),
  "",
  "## 4. Exploratory threshold -- NOT A VALIDATED CUTOFF",
  "",
  sprintf("**Exploratory only.** HRDsum >= %g describes the TCGA reference distribution. It is not a clinical threshold, least of all a pediatric one. The AUC below is the rank-based (Mann-Whitney) AUC of the *continuous* predicted HRDsum; the continuous score is **not** a probability and no probability mapping has been fitted or validated.",
          EXPLORATORY_THRESHOLD),
  "",
  sprintf("| Group | n | n obs >= %g | AUC (rank-based) | n obs within +/-%g | TP | FP | FN | TN |",
          EXPLORATORY_THRESHOLD, NEAR_THRESHOLD_WINDOW),
  "|---|---|---|---|---|---|---|---|---|",
  vapply(c(CNS_TYPES, "CNS_pooled"), md_thr_line, character(1)),
  "",
  "## 5. Pre-declared outcome (docs/28 section 7)",
  "",
  "The three outcomes were fixed in advance. The numeric readings used here were",
  "declared in the header of `scripts/analyze_cns.R`, written before unlock:",
  "",
  sprintf("* rank transfers: Spearman > 0 and bootstrap 95%% CI lower bound > 0 in BOTH lineages -- GBM %s, LGG %s -> **%s**",
          ifelse(rank_ok[["GBM"]], "yes", "no"), ifelse(rank_ok[["LGG"]], "yes", "no"),
          ifelse(RANK_TRANSFERS, "YES", "NO")),
  sprintf("* calibration transfers: beats the ORACLE within-group mean null, |bias| <= source %g-quantile (%s), and calibration slope inside the source range [%s, %s], in BOTH lineages -- GBM %s, LGG %s -> **%s**",
          BIAS_SOURCE_QUANTILE, fmt(bias_cut, 2), fmt(slope_lo), fmt(slope_hi),
          ifelse(cal_ok[["GBM"]], "yes", "no"), ifelse(cal_ok[["LGG"]], "yes", "no"),
          ifelse(CALIBRATION_TRANSFERS, "YES", "NO")),
  "",
  sprintf("### Outcome: %s", OUTCOME_NAME),
  "",
  sprintf("Pre-declared claim (docs/28 section 7): *%s*", OUTCOME_CLAIM),
  "",
  "In all three cases the deployment gate (`docs/27` section 2A) remains shut.",
  "Nothing here produces a validated clinical HRD assay. The CNS result cannot",
  "change the model (docs/28 section 6).",
  "")
writeLines(md, file.path(OUT_DIR, "cns_analysis_summary.md"))
message("  wrote cns_analysis_summary.md")

## ---- provenance --------------------------------------------------------------
file_hash <- function(p) {
  if (!file.exists(p)) return("file_absent")
  if (requireNamespace("digest", quietly = TRUE))
    return(digest::digest(file = p, algo = "sha256"))
  unname(tools::md5sum(p))
}
si <- file.path(OUT_DIR, "sessionInfo.txt")
con <- file(si, open = "wt")
writeLines(c(
  "KIDS26-Team12 -- locked CNS evaluation provenance",
  paste0("script                 : scripts/analyze_cns.R"),
  paste0("run label              : ", RUN_LABEL),
  paste0("timestamp              : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("protocol               : docs/28_CANDIDATE_SELECTION_AND_CNS_PROTOCOL.md section 6 (frozen)"),
  paste0("bootstrap seed         : ", BOOT_SEED),
  paste0("bootstrap resamples    : ", BOOT_R),
  paste0("bootstrap method       : percentile CI on resampled (observed, predicted) pairs"),
  paste0("exploratory threshold  : ", EXPLORATORY_THRESHOLD,
         " (NOT validated; near-threshold window +/-", NEAR_THRESHOLD_WINDOW, ")"),
  paste0("clipping rule (asserted, not applied here): ", CLIP_RULE),
  paste0("null baseline          : within-CNS-group tissue mean -- ORACLE, uses CNS labels"),
  paste0("cns predictions        : ", P_CNS, "  sha256=", file_hash(P_CNS)),
  paste0("master samples         : ", P_MASTER, "  sha256=", file_hash(P_MASTER)),
  paste0("source LOCO            : ", P_SRC, "  sha256=", file_hash(P_SRC)),
  paste0("n CNS samples          : ", nrow(d),
         "  (GBM ", sum(d$cancer_type == "GBM"), ", LGG ", sum(d$cancer_type == "LGG"), ")"),
  paste0("n source tissues       : ", N_SRC),
  paste0("outcome (docs/28 s7)   : ", OUTCOME_NAME),
  paste0("R                      : ", R.version.string),
  "",
  "---- sessionInfo() ----"), con)
capture.output(sessionInfo(), file = con, append = TRUE)
close(con)
message("  wrote ", si)

## =============================================================================
## 9. CONSOLE SUMMARY
## =============================================================================
line <- function(...) message(sprintf(...))
message("")
message("================= LOCKED CNS EVALUATION -- ", RUN_LABEL, " =================")
line("cohort            : %d samples (GBM %d, LGG %d), partition = %s",
     nrow(d), sum(d$cancer_type == "GBM"), sum(d$cancer_type == "LGG"), LOCKED_LABEL)
line("clipping (checked): %s", CLIP_RULE)
message("")
message("-- PRIMARY, per lineage ------------------------------------------------")
for (g in c(CNS_TYPES, "CNS_pooled")) {
  line("%-11s n=%4d  r=%s [%s,%s]  rho=%s [%s,%s]  MAE=%s  bias=%s  slope=%s  sd(p)/sd(o)=%s",
       g, get_g(g, "n"),
       fmt(get_g(g, "Pearson")), fmt(get_g(g, "Pearson_lo95")), fmt(get_g(g, "Pearson_hi95")),
       fmt(get_g(g, "Spearman")), fmt(get_g(g, "Spearman_lo95")), fmt(get_g(g, "Spearman_hi95")),
       fmt(get_g(g, "MAE"), 2), fmt(get_g(g, "bias_mean_signed_error"), 2),
       fmt(get_g(g, "calibration_slope")), fmt(get_g(g, "range_compression_ratio")))
}
message("")
message("-- ORACLE null (within-group tissue mean; uses CNS labels, unavailable at deployment)")
for (g in c(CNS_TYPES, "CNS_pooled")) {
  line("%-11s model MAE=%s  oracle-null MAE=%s  skill=%s", g,
       fmt(get_g(g, "MAE"), 2), fmt(get_g(g, "MAE_oracle_tissue_mean_null"), 2),
       fmt(get_g(g, "skill_vs_oracle_tissue_mean")))
}
message("")
message("-- PLACEMENT in the ", N_SRC, " source LOCO tissues (would CNS look unusual as the 31st?)")
for (i in seq_len(nrow(placement))) {
  r <- placement[i]
  line("%-18s %-4s value=%8s  source median=%8s  pct=%5s%%  %s%s",
       r$statistic, r$group, fmt(r$value), fmt(r$source_median),
       fmt(r$percentile_in_source, 1), r$rank_statement,
       if (isTRUE(r$inside_source_range)) "" else "  [OUTSIDE source range]")
}
message("")
message("-- OOD")
for (g in c(CNS_TYPES, "CNS_pooled")) {
  r <- ood_summary[group == g & subset == "all"]
  line("%-11s reportable=%d  flagged=%d  ood_score min/q25/med/q75/max = %s/%s/%s/%s/%s",
       g, r$n_reportable[1], r$n_flagged[1], fmt(r$ood_min[1], 2), fmt(r$ood_q25[1], 2),
       fmt(r$ood_median[1], 2), fmt(r$ood_q75[1], 2), fmt(r$ood_max[1], 2))
  for (s in c("reportable", "flagged")) {
    rr <- ood_summary[group == g & subset == s]
    line("   %-10s n=%4d  MAE=%s  r=%s  rho=%s", s, rr$n[1], fmt(rr$MAE[1], 2),
         fmt(rr$Pearson[1]), fmt(rr$Spearman[1]))
  }
}
message("")
message("!! EXPLORATORY THRESHOLD SECTION -- NOT A VALIDATED CUTOFF !!")
message("!! HRDsum >= ", EXPLORATORY_THRESHOLD, " describes the TCGA reference distribution. The AUC below is")
message("!! rank-based on the CONTINUOUS prediction. The score is NOT a probability and")
message("!! must never be reported as one. No decision rule is being proposed.")
for (g in c(CNS_TYPES, "CNS_pooled")) {
  r <- threshold_tab[group == g]
  line("%-11s AUC=%s  n>=%g: %d  near-threshold(+/-%g): %d  crossings TP/FP/FN/TN = %d/%d/%d/%d",
       g, fmt(r$AUC_rank_based_continuous_score[1]), EXPLORATORY_THRESHOLD,
       r$n_observed_above[1], NEAR_THRESHOLD_WINDOW, r$n_observed_near_threshold[1],
       r$crossings_TP[1], r$crossings_FP[1], r$crossings_FN[1], r$crossings_TN[1])
}
message("")
message("-- PRE-DECLARED OUTCOME (docs/28 section 7), computed from the numbers above")
line("rank transfers        : %s  (GBM %s, LGG %s)", ifelse(RANK_TRANSFERS, "YES", "NO"),
     ifelse(rank_ok[["GBM"]], "yes", "no"), ifelse(rank_ok[["LGG"]], "yes", "no"))
line("calibration transfers : %s  (GBM %s, LGG %s)", ifelse(CALIBRATION_TRANSFERS, "YES", "NO"),
     ifelse(cal_ok[["GBM"]], "yes", "no"), ifelse(cal_ok[["LGG"]], "yes", "no"))
message("")
message("PLAIN LANGUAGE: ", OUTCOME_NAME, ".")
message("  ", OUTCOME_CLAIM)
message("  The deployment gate (docs/27 section 2A) remains shut either way, and the")
message("  CNS result cannot change the model (docs/28 section 6).")
message("")
message("== analyze_cns.R complete ==  outputs in ", OUT_DIR)

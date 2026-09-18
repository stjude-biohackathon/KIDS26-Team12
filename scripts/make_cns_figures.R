# =============================================================================
# scripts/make_cns_figures.R -- figures for the locked CNS evaluation
# -----------------------------------------------------------------------------
# USAGE (from the repository root)
#   Rscript scripts/make_cns_figures.R \
#     <cns_predictions.tsv> <master_samples.tsv> <source_loco_predictions.tsv> \
#     [out_dir] [label]
#
#   e.g. Rscript scripts/make_cns_figures.R \
#          results/locked_cns/cns_predictions.tsv \
#          data/processed/master_samples.tsv \
#          results/loco_run01/loco_predictions.tsv \
#          results/figures primary_run01
#
# PURPOSE
#   The visual half of the frozen CNS test (docs/28 section 6). It recomputes the
#   same quantities scripts/analyze_cns.R writes to TSV, from the same inputs and
#   the same shared metrics() in R/model.R, so the figures and the tables cannot
#   disagree. It fits nothing, recalibrates nothing and re-clips nothing:
#   predictions arrive already C2-clipped and that is asserted, not assumed.
#
# INPUTS (all required; the script fails loudly if any is missing)
#   1  cns_predictions.tsv          output of scripts/predict_frozen.R on the
#                                   locked partition (authoritative clipped
#                                   prediction + clipping_rule + reportable).
#   2  master_samples.tsv           supplies cancer_type, partition, HRDsum.
#   3  source_loco_predictions.tsv  results/loco_run01/loco_predictions.tsv,
#                                   the 30 development tissues.
#   4  out_dir                      default results/figures
#   5  label                        optional run label, printed in the captions.
#   R/figures.R                     shared theme, palette, save_fig(), fmt()
#   R/model.R                       metrics() -- identical code path to the tables
#
# OUTPUTS (out_dir, default results/figures/)
#   fig07_cns_pred_vs_obs.png   predicted vs observed HRDsum, GBM and LGG as
#                               separate panels, identity line, per-panel
#                               r / rho / MAE annotation.
#   fig08_cns_placement.png     THE MOST IMPORTANT CNS FIGURE. The source-LOCO
#                               distributions of Pearson, MAE and bias shown as
#                               histograms over the 30 held-out tissues, with GBM
#                               and LGG overlaid as marked vertical lines. Answers
#                               "would CNS look unusual as the 31st held-out
#                               tissue?" (docs/28 section 6).
#   fig07_08_sessionInfo.txt    seed, inputs and sessionInfo() for these figures.
#
# NOTE ON THE CNS LOCK
#   R/figures.R::assert_development() exists to keep GBM/LGG out of the
#   development figure layer. This script is the ONE deliberate exception,
#   executed only after the locked cohort has been opened under docs/28. It
#   therefore applies assert_development() to the SOURCE LOCO data (which must
#   still be CNS-free, or the placement comparison is circular) and never to the
#   CNS frame itself.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

## ---- setup -------------------------------------------------------------------
if (!file.exists("R/figures.R") || !file.exists("R/model.R")) {
  stop("Run this script from the repository root (R/figures.R and R/model.R not ",
       "found at ", normalizePath(".", mustWork = FALSE), ")", call. = FALSE)
}
source("R/figures.R")   # theme_kids26(), PAL, read_required(), save_fig(), fmt(), FIG_SEED
source("R/model.R")     # metrics() -- the same function the tables use
set.seed(FIG_SEED)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L || length(args) > 5L) {
  stop("Expected <cns_predictions.tsv> <master_samples.tsv> ",
       "<source_loco_predictions.tsv> [out_dir] [label]", call. = FALSE)
}
P_CNS    <- args[1]
P_MASTER <- args[2]
P_SRC    <- args[3]
FIG_DIR  <- if (length(args) >= 4L) args[4] else file.path("results", "figures")
RUN_LABEL <- if (length(args) >= 5L) args[5] else "locked CNS evaluation"

CNS_TYPES    <- c("GBM", "LGG")
LOCKED_LABEL <- "locked_CNS"
PRED_COL     <- "predicted_reference_HRDsum"

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
message("== make_cns_figures.R ==  seed = ", FIG_SEED, "  label = ", RUN_LABEL)
message("-- reading inputs")

cns    <- read_required(P_CNS)
master <- read_required(P_MASTER)
src    <- read_required(P_SRC)

need_cols <- function(dt, cols, what) {
  miss <- setdiff(cols, names(dt))
  if (length(miss)) stop(what, " is missing required column(s): ",
                         paste(miss, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}
need_cols(cns, c("sample_id", PRED_COL, "clipping_rule"), "CNS predictions file")
need_cols(master, c("sample_id", "cancer_type", "partition", "HRDsum"), "master sample table")
need_cols(src, c("cancer_type", "actual", PRED_COL), "source LOCO predictions")

## ---- join, cohort gate, C2 gate ----------------------------------------------
idx <- match(cns$sample_id, master$sample_id)
if (anyNA(idx)) {
  stop(sprintf("%d of %d predicted samples are absent from the master table; ",
               sum(is.na(idx)), nrow(cns)),
       "refusing to plot samples whose lock status and label cannot be established.",
       call. = FALSE)
}
d <- data.table::data.table(
  sample_id   = cns$sample_id,
  predicted   = as.numeric(cns[[PRED_COL]]),
  cancer_type = as.character(master$cancer_type[idx]),
  partition   = as.character(master$partition[idx]),
  observed    = as.numeric(master$HRDsum[idx]),
  clipping_rule = as.character(cns$clipping_rule))

if (!all(d$partition == LOCKED_LABEL)) {
  stop("COHORT GATE: not every row has partition == '", LOCKED_LABEL,
       "'. This script plots the locked CNS partition only.", call. = FALSE)
}
if (!all(d$cancer_type %chin% CNS_TYPES) || !all(CNS_TYPES %chin% unique(d$cancer_type))) {
  stop("COHORT GATE: expected exactly {GBM, LGG}; found: ",
       paste(sort(unique(d$cancer_type)), collapse = ", "), call. = FALSE)
}
cr <- unique(d$clipping_rule)
if (length(cr) != 1L || is.na(cr) || !grepl("^C2_ADOPTED:", cr)) {
  stop("C2 GATE: clipping_rule is absent, mixed, or does not record the adopted ",
       "C2 rule (found: ", paste(cr, collapse = " | "), "). These predictions were ",
       "not produced by a C2-compliant inference run and must not be plotted.",
       call. = FALSE)
}
if (any(is.finite(d$predicted) & d$predicted < 0)) {
  stop("C2 GATE: negative authoritative predictions despite clipping_rule = ", cr,
       call. = FALSE)
}
CLIP_RULE <- cr
d <- d[is.finite(observed) & is.finite(predicted)]
if (!nrow(d)) stop("No CNS rows with both a finite prediction and a finite HRDsum.",
                   call. = FALSE)
message(sprintf("  CNS cohort: %d rows (GBM %d, LGG %d); clipping rule = %s",
                nrow(d), sum(d$cancer_type == "GBM"), sum(d$cancer_type == "LGG"), CLIP_RULE))

# The source distribution must be CNS-free or the placement comparison is circular.
assert_development(src, "cancer_type", "source LOCO predictions")
src <- data.table::copy(src)
src[, pred_clip := pmax(as.numeric(get(PRED_COL)), 0)]   # same C2 rule as the CNS side
src[, actual := as.numeric(actual)]

## ---- per-group and per-tissue metrics (same metrics() as the tables) ---------
cns_stats <- data.table::rbindlist(lapply(CNS_TYPES, function(g) {
  s <- d[cancer_type == g]
  mm <- metrics(s$observed, s$predicted)
  data.frame(group = g, n = mm$n, Pearson = mm$Pearson, Spearman = mm$Spearman,
             MAE = mm$MAE, bias = mm$bias, stringsAsFactors = FALSE)
}))
src_tissue <- data.table::rbindlist(lapply(split(src, src$cancer_type), function(s) {
  mm <- metrics(s$actual, s$pred_clip)
  data.frame(cancer_type = s$cancer_type[1], Pearson = mm$Pearson,
             MAE = mm$MAE, bias = mm$bias, stringsAsFactors = FALSE)
}))
N_SRC <- nrow(src_tissue)
message(sprintf("  source LOCO distribution: %d held-out tissues", N_SRC))

## =============================================================================
## fig07 -- predicted vs observed, GBM and LGG panels
## =============================================================================
message("-- fig07_cns_pred_vs_obs")
lim <- range(c(d$observed, d$predicted), finite = TRUE)
pad <- 0.04 * diff(lim)
lim <- c(min(0, lim[1] - pad), lim[2] + pad)

ann <- data.table::copy(cns_stats)
ann[, panel := sprintf("%s  (n = %d)", group, n)]
ann[, txt := sprintf("r = %s\nrho = %s\nMAE = %s", fmt(Pearson, 3), fmt(Spearman, 3), fmt(MAE, 2))]
d[, panel := sprintf("%s  (n = %d)", cancer_type,
                     ann$n[match(cancer_type, ann$group)])]

fig07 <- ggplot(d, aes(x = predicted, y = observed)) +
  geom_abline(slope = 1, intercept = 0, colour = PAL[["neutral"]],
              linetype = "dashed", linewidth = 0.5) +
  geom_point(colour = PAL[["accent"]], alpha = 0.45, size = 1.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = PAL[["warn"]], linewidth = 0.6) +
  geom_text(data = ann, aes(x = lim[1] + 0.03 * diff(lim), y = lim[2] - 0.03 * diff(lim),
                            label = txt), inherit.aes = FALSE,
            hjust = 0, vjust = 1, size = 3.6, colour = "grey20", lineheight = 1.15) +
  facet_wrap(~ panel, nrow = 1) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(
    title = "Zero-shot prediction into the locked CNS lineage",
    subtitle = paste0("Frozen model, one-shot locked evaluation (", RUN_LABEL,
                      "). Dashed line = identity; red line = lm(observed ~ predicted)."),
    x = "Predicted HRDsum (frozen model, C2-clipped)",
    y = "Observed HRDsum (reference, SNP6 + ABSOLUTE)",
    caption = paste0(
      "Source: ", basename(P_CNS), " scored by scripts/predict_frozen.R; labels from ",
      basename(P_MASTER), ".\n",
      "Clipping rule applied at inference and asserted here: ", CLIP_RULE,
      ". No refitting, recalibration or re-clipping. Protocol: docs/28 section 6.")) +
  theme_kids26(base_size = 13)
save_fig(fig07, file.path(FIG_DIR, "fig07_cns_pred_vs_obs.png"), width = 10.0, height = 5.6)

## =============================================================================
## fig08 -- placement of GBM and LGG in the source LOCO distribution
## =============================================================================
# The most important CNS figure. Each panel is the empirical distribution of one
# statistic across the N_SRC source tissues, each of which was genuinely held
# out. GBM and LGG are overlaid as vertical lines. The question the figure is
# built to answer is docs/28 section 6's: would CNS look unusual if it had simply
# been the 31st held-out tissue? No pass/fail line is drawn, by design.
message("-- fig08_cns_placement")

STATS <- c(Pearson = "Pearson r (within held-out tissue)",
           MAE     = "MAE (HRDsum units)",
           bias    = "Signed bias (predicted - observed)")

# Percentile is computed exactly as in scripts/analyze_cns.R, oriented so that
# higher is always better: share of source tissues the CNS group is no worse than.
pct_of <- function(stat, v) {
  sv <- src_tissue[[stat]]
  if (!is.finite(v)) return(NA_real_)
  if (stat == "Pearson") 100 * mean(sv <= v, na.rm = TRUE)
  else if (stat == "MAE") 100 * mean(sv >= v, na.rm = TRUE)
  else 100 * mean(abs(sv) >= abs(v), na.rm = TRUE)   # bias ranked on |bias|
}

src_long <- data.table::rbindlist(lapply(names(STATS), function(s)
  data.table::data.table(statistic = STATS[[s]], value = src_tissue[[s]],
                         cancer_type = src_tissue$cancer_type)))
cns_long <- data.table::rbindlist(lapply(names(STATS), function(s)
  data.table::data.table(
    statistic = STATS[[s]], group = cns_stats$group, value = cns_stats[[s]],
    pct = vapply(cns_stats[[s]], function(v) pct_of(s, v), numeric(1)))))
cns_long[, lab := sprintf("%s: %s  (%sth pct)", group, fmt(value, 3), fmt(pct, 0))]
src_long[, statistic := factor(statistic, levels = unname(STATS))]
cns_long[, statistic := factor(statistic, levels = unname(STATS))]

# A zero reference only means something for the two signed statistics.
zero_ref <- data.table::data.table(
  statistic = factor(c(STATS[["Pearson"]], STATS[["bias"]]), levels = unname(STATS)),
  x = c(0, 0))

fig08 <- ggplot(src_long, aes(x = value)) +
  geom_histogram(bins = 12, fill = PAL[["light"]], colour = "grey45", linewidth = 0.25) +
  geom_rug(colour = PAL[["neutral"]], alpha = 0.7, length = grid::unit(0.03, "npc")) +
  geom_vline(data = zero_ref, aes(xintercept = x), inherit.aes = FALSE,
             colour = "grey55", linetype = "dotted", linewidth = 0.4) +
  geom_vline(data = cns_long, aes(xintercept = value, colour = group),
             inherit.aes = FALSE, linewidth = 1.0) +
  geom_text(data = cns_long, aes(x = value, y = Inf, label = lab, colour = group),
            inherit.aes = FALSE, angle = 90, hjust = 1.05, vjust = -0.35,
            size = 3.2, show.legend = FALSE) +
  facet_wrap(~ statistic, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = c(GBM = unname(PAL[["warn"]]),
                                 LGG = unname(PAL[["amber"]])), name = NULL) +
  labs(
    title = "Would CNS look unusual as the 31st held-out tissue?",
    subtitle = paste0("Grey histogram = the ", N_SRC,
                      " source tissues, each genuinely held out in the LOCO run.\n",
                      "Vertical lines = the locked CNS lineages (", RUN_LABEL, ")."),
    x = NULL, y = "Number of source tissues",
    caption = paste0(
      "Per-tissue values recomputed from ", basename(P_SRC),
      " under the same C2 clipping and the same metrics() in R/model.R as the CNS values.\n",
      "Percentiles are oriented so higher is always better (for MAE and |bias|, the share of source tissues the CNS group is no worse than).\n",
      "No pass/fail line is drawn: docs/28 section 6 fixes placement inside this distribution as the interpretive frame, not a threshold.")) +
  theme_kids26(base_size = 13) +
  theme(panel.grid.major.x = element_blank())
save_fig(fig08, file.path(FIG_DIR, "fig08_cns_placement.png"), width = 12.0, height = 5.4)

## ---- provenance ---------------------------------------------------------------
si <- file.path(FIG_DIR, "fig07_08_sessionInfo.txt")
con <- file(si, open = "wt")
writeLines(c(
  "KIDS26-Team12 -- locked CNS figures provenance",
  paste0("script            : scripts/make_cns_figures.R"),
  paste0("run label         : ", RUN_LABEL),
  paste0("timestamp         : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("seed (FIG_SEED)   : ", FIG_SEED),
  paste0("cns predictions   : ", P_CNS),
  paste0("master samples    : ", P_MASTER),
  paste0("source LOCO       : ", P_SRC, "  (", N_SRC, " tissues)"),
  paste0("clipping rule     : ", CLIP_RULE, "  (applied at inference, asserted here)"),
  paste0("cohort            : ", nrow(d), " locked_CNS rows (GBM ",
         sum(d$cancer_type == "GBM"), ", LGG ", sum(d$cancer_type == "LGG"), ")"),
  paste0("R                 : ", R.version.string),
  "",
  "---- CNS per-group statistics plotted ----",
  capture.output(print(as.data.frame(cns_stats))),
  "",
  "---- sessionInfo() ----",
  capture.output(sessionInfo())), con)
close(con)
message("  wrote ", si)

message("== make_cns_figures.R complete ==  figures in ", FIG_DIR)

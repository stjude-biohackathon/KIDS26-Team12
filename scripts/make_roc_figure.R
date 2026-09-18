# =============================================================================
# scripts/make_roc_figure.R - ROC both ways for the frozen Candidate A model
# =============================================================================
#
# USAGE
#   Rscript scripts/make_roc_figure.R
#
# WHAT THIS IS
#   The flash-talk ROC figure, built from Candidate A - the locked-down
#   pooled/run-01 model (results/frozen_2026-09-18/frozen_nonCNS.rds,
#   sha256 0e975cd6..., alpha 0.1, lambda 1.4201, 774 non-zero coefficients).
#   NOT V3-abs.
#
#   It scores the DEVELOPMENT LOCO predictions - 7,065 adult TCGA samples across
#   30 non-CNS cancer types, each predicted by a model that never saw its cancer
#   type during training. The locked CNS cohort (642 GBM+LGG) is NOT touched by
#   this script and remains unscored.
#
# WHY THREE CURVES AND NOT ONE
#   Panel A pools all tissues at the exploratory HRDsum >= 42 cutoff. That is
#   the number that looks best, and on its own it is misleading: HRDsum varies
#   strongly BETWEEN tumour types, so a model that only recognised the lineage
#   would still score well. The dashed curve is exactly that control - every
#   sample replaced by its own tissue's MEAN prediction, which contains no
#   within-patient information whatsoever. The gap between solid and dashed is
#   what methylation contributes beyond knowing the tissue.
#
#   Panel B asks the question the project can actually answer: within a single
#   cancer type, does the score rank the higher-HRD tumours above the lower
#   ones? Predictions are centred on their own tissue mean, and the label is
#   that tissue's own top quartile. This removes the between-tissue axis
#   entirely, so the lineage shortcut is unavailable by construction.
#
# INPUTS
#   results/loco_run01/loco_predictions.tsv
# OUTPUTS
#   results/figures/fig09_roc_both_ways.png
#   results/tables/roc_summary.tsv
#
# CAVEATS THAT MUST TRAVEL WITH THIS FIGURE
#   * HRDsum >= 42 is an EXPLORATORY cutoff, not a validated clinical
#     threshold, and not one this project prospectively justified for
#     methylation. Only 10.6% of the development cohort sits above it.
#   * The score is a continuous regression output. It is NOT a probability and
#     must never be presented as one.
#   * This is a DEVELOPMENT result. Clipping, the log1p rejection and the whole
#     C1 investigation were all developed against these same folds.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})
source("R/figures.R")

SEED <- 4827
set.seed(SEED)

IN   <- "results/loco_run01/loco_predictions.tsv"
OUTF <- "results/figures/fig09_roc_both_ways.png"
OUTT <- "results/tables/roc_summary.tsv"
if (!file.exists(IN)) stop("Missing input: ", IN)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables",  showWarnings = FALSE, recursive = TRUE)

d <- fread(IN)
if (any(d$cancer_type %in% c("GBM", "LGG")))
  stop("CNS rows present in a development figure - refusing to plot.")

# Adopted C2 rule: negative HRDsum is physically meaningless.
d[, pred_clip := pmax(predicted_reference_HRDsum, 0)]

# --- ROC/AUC helpers --------------------------------------------------------
auc_mw <- function(score, pos) {
  pos <- as.logical(pos)
  r <- rank(score)
  (sum(r[pos]) - sum(pos) * (sum(pos) + 1) / 2) / (sum(pos) * sum(!pos))
}
roc_pts <- function(score, pos, label) {
  pos <- as.logical(pos)
  o <- order(score, decreasing = TRUE)
  tp <- cumsum(pos[o]); fp <- cumsum(!pos[o])
  data.table(fpr = c(0, fp / sum(!pos)), tpr = c(0, tp / sum(pos)), curve = label)
}
# Stratified bootstrap CI for AUC: resample within class so prevalence is fixed.
boot_auc <- function(score, pos, B = 2000) {
  pos <- as.logical(pos); ip <- which(pos); in_ <- which(!pos)
  v <- replicate(B, {
    s <- c(sample(ip, length(ip), TRUE), sample(in_, length(in_), TRUE))
    auc_mw(score[s], pos[s])
  })
  quantile(v, c(0.025, 0.975), names = FALSE)
}

# --- Panel A: absolute, pooled ---------------------------------------------
d[, hi_abs := actual >= 42]
d[, tissue_mean_pred := mean(pred_clip), by = cancer_type]

a_model  <- auc_mw(d$pred_clip,        d$hi_abs)
a_tissue <- auc_mw(d$tissue_mean_pred, d$hi_abs)
ci_a     <- boot_auc(d$pred_clip, d$hi_abs)

pa <- rbind(
  roc_pts(d$pred_clip,        d$hi_abs, sprintf("Methylation score (AUC %.3f)", a_model)),
  roc_pts(d$tissue_mean_pred, d$hi_abs, sprintf("Tissue identity only (AUC %.3f)", a_tissue))
)

# --- Panel B: within-tissue top quartile ------------------------------------
d[, q75 := quantile(actual, 0.75), by = cancer_type]
d[, hi_within := actual > q75]
d[, pred_centred := pred_clip - mean(pred_clip), by = cancer_type]

b_model <- auc_mw(d$pred_centred, d$hi_within)
ci_b    <- boot_auc(d$pred_centred, d$hi_within)
pb <- roc_pts(d$pred_centred, d$hi_within,
              sprintf("Tissue-centred score (AUC %.3f)", b_model))

# --- Plot -------------------------------------------------------------------
mk <- function(dat, title, sub) {
  ggplot(dat, aes(fpr, tpr, colour = curve, linetype = curve)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_line(linewidth = 0.9) +
    scale_colour_manual(values = c("#B2182B", "grey45")) +
    scale_linetype_manual(values = c("solid", "22")) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(title = title, subtitle = sub,
         x = "False positive rate", y = "True positive rate") +
    theme_kids26() +
    theme(legend.position = c(0.97, 0.03), legend.justification = c(1, 0),
          legend.title = element_blank(),
          legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
          legend.key.width = unit(1.6, "lines"), plot.title.position = "plot")
}

gA <- mk(pa, "A. Absolute cutoff, pooled across tissues",
         sprintf("HRDsum >= 42 (exploratory) | %d of %d positive (%.1f%%)\nMost of this is lineage: tissue identity alone reaches %.3f",
                 sum(d$hi_abs), nrow(d), 100 * mean(d$hi_abs), a_tissue)) +
  theme(plot.margin = margin(5.5, 22, 5.5, 5.5))
gB <- mk(pb, "B. Within-tissue top quartile",
         sprintf("Each tissue's own upper quartile | %d of %d positive (%.1f%%)\nLineage shortcut removed by construction",
                 sum(d$hi_within), nrow(d), 100 * mean(d$hi_within))) +
  theme(plot.margin = margin(5.5, 5.5, 5.5, 22))

if (!requireNamespace("patchwork", quietly = TRUE))
  stop("patchwork is required for the two-panel layout; install via scripts/setup.R")

g <- patchwork::wrap_plots(gA, gB, nrow = 1) +
  patchwork::plot_annotation(
    title = "Methylation-derived HRD score: ranking performance two ways",
    subtitle = paste0("Frozen Candidate A: pooled features, alpha 0.1, 774 non-zero coefficients.\n",
                      "7,065 adult TCGA samples, 30 non-CNS types, each held out by cancer type.\n",
                      "Development result - not a validated assay, not a probability. Locked CNS cohort not scored."),
    theme = theme_kids26()
  )
save_fig(g, OUTF, width = 13.5, height = 6.6)

# --- Table ------------------------------------------------------------------
tab <- data.table(
  framing = c("Absolute HRDsum >= 42 (pooled)",
              "Tissue-identity-only control (pooled)",
              "Within-tissue top quartile (tissue-centred)"),
  n_positive = c(sum(d$hi_abs), sum(d$hi_abs), sum(d$hi_within)),
  prevalence = sprintf("%.1f%%", 100 * c(mean(d$hi_abs), mean(d$hi_abs), mean(d$hi_within))),
  AUC = round(c(a_model, a_tissue, b_model), 4),
  ci_low  = round(c(ci_a[1], NA, ci_b[1]), 4),
  ci_high = round(c(ci_a[2], NA, ci_b[2]), 4),
  note = c("Exploratory cutoff; inflated by between-tissue HRD differences",
           "No within-patient information; the shortcut a pooled AUC rewards",
           "The defensible ranking claim")
)
fwrite(tab, OUTT, sep = "\t")

cat("\n"); print(tab)
cat(sprintf("\nPooled AUC %.3f, of which tissue identity alone supplies %.3f.\n", a_model, a_tissue))
cat(sprintf("Within-tissue AUC %.3f [%.3f, %.3f] - the claim that survives.\n", b_model, ci_b[1], ci_b[2]))
cat("Written:\n  ", OUTF, "\n  ", OUTT, "\n", sep = "")

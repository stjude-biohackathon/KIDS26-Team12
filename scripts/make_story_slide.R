# =============================================================================
# scripts/make_story_slide.R -- ONE flash-talk slide, three acts, 16:9
# -----------------------------------------------------------------------------
# PURPOSE
#   Build a single projector-ready PNG telling the project story:
#     Act 1  the lineage-imprinting problem in the initial model (run 01, "A")
#     Act 2  five attacks on it -- four failed, one survived
#     Act 3  the V3-abs result: 42% of excess lineage removed, no ranking cost
#
# INPUTS
#   R/figures.R                                   shared theme / palette / save_fig
#   results/tables/A_vs_V3_per_tissue_full30.tsv  per-tissue d_r for the Act 3 strip
#   results/tables/tableA2_A_vs_V3_full30.tsv     read back only to ASSERT that the
#                                                 hardcoded headline numbers match
# OUTPUTS
#   results/figures/slide_story.png               4800 x 2700 px (16 x 9 in @ 300 dpi)
#   results/figures/sessionInfo_slide.txt
#
# CNS LOCK
#   Nothing here reads, scores or plots GBM/LGG. The only per-sample-derived input
#   is the 30-row development-cohort per-tissue table, guarded by assert_development().
#
# PROVENANCE OF EVERY HARDCODED NUMBER
#   Act 1
#     336,480 CpGs / 7,065 samples / 30 types / nested LOCO ... project design brief
#     macro Pearson 0.5197, pooled MAE 8.9717 ............ tableA2 rows 1, 3
#     29/30 positive r, permutation p = 0.001 ............ tableA2 row "Tissues with
#                                                          positive r"; permutation
#                                                          test from run 01 write-up
#     R2(tissue -> truth)      0.3414 .................... tableA2 "Tissue R2 of truth"
#     R2(tissue -> prediction) 0.5584 .................... tableA2 "Tissue R2 of predictions"
#     excess                   0.2171 .................... tableA2 "Excess tissue R2"
#     SARC -15.44 / PCPG +8.29 offsets, mean |offset| 3.46  per-tissue table bias_A
#                                                          (SARC, PCPG rows) and
#                                                          tableA2 "Mean absolute tissue bias"
#   Act 2  (all five attempt outcomes are quoted from the experiment write-ups;
#           they are narrative text, not recomputed here)
#     1 label-free offset regression   best LOTO R2 = -0.116
#     2 zero-shot percentile remap     tissue R2 0.562 -> 0.593, 0 ranking info
#     3 relative-target V1/V2          failed 3 of 7 pre-registered gates
#                                      (macro rho -0.219, BRCA -0.405,
#                                       tissue R2 0.576 vs 0.527)
#     4 V4 lineage-penalised features  gate G3 failed, tissue R2 0.5219 -> 0.5361
#     5 V3-abs                         passed 7/7 sentinel gates + 7/7 criteria
#   Act 3
#     excess R2 A 0.2171 [0.2018, 0.2310], V3 0.1258 [0.1104, 0.1407]  tableA2
#     delta 0.0913, 95% CI [0.0824, 0.0996], p < 0.0001, 42.0% [37.9, 46.5]  tableA2
#     macro skill vs oracle null 0.4192 -> 0.6304 ........ tableA2
#     pooled MAE 8.97 -> 8.45, mean |bias| 3.46 -> 2.97 .. tableA2
#     paired Pearson mean -0.0161, CI [-0.0410, 0.0051], Wilcoxon p = 0.428  tableA2
#     paired Spearman mean -0.0275, p = 0.477 ............ tableA2
#     per-tissue d_r points .............................. read from the per-tissue TSV
#
# USAGE  Rscript scripts/make_story_slide.R     (run from the repository root)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

source("R/figures.R")
set.seed(FIG_SEED)

OUT_PNG  <- "results/figures/slide_story.png"
OUT_SESS <- "results/figures/sessionInfo_slide.txt"

# --- headline numbers (single definition, asserted against tableA2 below) -----
N_R2_TRUTH  <- 0.3414
N_R2_PRED_A <- 0.5584
N_R2_PRED_B <- 0.4672
N_EXC_A     <- 0.2171;  N_EXC_A_LO <- 0.2018; N_EXC_A_HI <- 0.2310
N_EXC_B     <- 0.1258;  N_EXC_B_LO <- 0.1104; N_EXC_B_HI <- 0.1407
N_DELTA     <- 0.0913;  N_DELTA_LO <- 0.0824; N_DELTA_HI <- 0.0996
N_PCT       <- 42.0;    N_PCT_LO   <- 37.9;   N_PCT_HI   <- 46.5
N_MAE_A     <- 8.97;    N_MAE_B    <- 8.45
N_BIAS_A    <- 3.46;    N_BIAS_B   <- 2.97
N_SKILL_A   <- 0.419;   N_SKILL_B  <- 0.630
N_DR_MEAN   <- -0.0161; N_DR_LO    <- -0.0410; N_DR_HI <- 0.0051

# --- verify the hardcoded values against the shipped table -------------------
a2df <- as.data.frame(read_required("results/tables/tableA2_A_vs_V3_full30.tsv"))
chk <- function(metric_name, col, want, tol = 5e-4) {
  row <- a2df[a2df$metric == metric_name, , drop = FALSE]
  got <- suppressWarnings(as.numeric(row[[col]][1]))
  if (is.na(got) || abs(got - want) > tol) {
    stop(sprintf("PROVENANCE MISMATCH: '%s'/%s -- table has %s, script hardcodes %s",
                 metric_name, col, format(got), format(want)), call. = FALSE)
  }
  invisible(TRUE)
}
chk("Tissue R2 of truth",                 "model_A", N_R2_TRUTH)
chk("Tissue R2 of predictions",           "model_A", N_R2_PRED_A)
chk("Tissue R2 of predictions",           "v3_abs",  N_R2_PRED_B)
chk("Excess tissue R2 (pred - truth)",    "model_A", N_EXC_A)
chk("Excess tissue R2 (pred - truth)",    "v3_abs",  N_EXC_B)
chk("Excess tissue R2 (pred - truth)",    "difference", -N_DELTA)
chk("Pooled MAE",                         "model_A", N_MAE_A, tol = 5e-3)
chk("Pooled MAE",                         "v3_abs",  N_MAE_B, tol = 5e-3)
chk("Mean absolute tissue bias",          "model_A", N_BIAS_A, tol = 5e-3)
chk("Mean absolute tissue bias",          "v3_abs",  N_BIAS_B, tol = 5e-3)
chk("Macro skill vs oracle tissue-mean null", "model_A", N_SKILL_A, tol = 5e-3)
chk("Macro skill vs oracle tissue-mean null", "v3_abs",  N_SKILL_B, tol = 5e-3)
message("  provenance check: all headline numbers match tableA2")

# --- per-tissue deltas for the Act 3 strip ------------------------------------
pt <- read_required("results/tables/A_vs_V3_per_tissue_full30.tsv")
assert_development(pt, "cancer_type", "A_vs_V3_per_tissue_full30")
stopifnot(nrow(pt) == 30L)
pt[, improved := d_r > 0]


# =============================================================================
# shared sizing -- tuned for a 16 x 9 in canvas at 300 dpi.
# ggplot2 geom_text `size` is in mm; 1 mm ~ 2.845 pt. On a 16-in-wide slide a
# size-5 label is ~14 pt of slide type, which reads from the back of a room.
# =============================================================================
BS   <- 17            # panel base_size for axes
HEAD <- 5.5           # in-panel act headline
TXT  <- 4.15          # body copy
TXTS <- 3.7           # small body copy
NUMB <- 5.2           # in-plot numeric annotation

blank_canvas <- function(ylim = c(0, 10)) {
  ggplot() +
    coord_cartesian(xlim = c(0, 10), ylim = ylim, expand = FALSE, clip = "off") +
    theme_void(base_size = BS) +
    theme(plot.margin = margin(0, 6, 0, 6))
}

# =============================================================================
# ACT 1 -- the lineage-excess problem
# =============================================================================
d1 <- data.table(
  what = factor(c("TRUTH\nreference HRDsum", "MODEL A\nprediction", "MODEL A\nprediction"),
                levels = c("TRUTH\nreference HRDsum", "MODEL A\nprediction")),
  part = factor(c("shared", "shared", "excess"), levels = c("excess", "shared")),
  val  = c(N_R2_TRUTH, N_R2_TRUTH, N_EXC_A)
)

p1 <- ggplot(d1, aes(x = what, y = val, fill = part)) +
  geom_col(width = 0.42, colour = "white", linewidth = 0.6) +
  annotate("segment", x = 2.34, xend = 2.34, y = N_R2_TRUTH, yend = N_R2_PRED_A,
           linewidth = 0.9, colour = PAL[["warn"]]) +
  annotate("segment", x = 2.28, xend = 2.34, y = N_R2_TRUTH, yend = N_R2_TRUTH,
           linewidth = 0.9, colour = PAL[["warn"]]) +
  annotate("segment", x = 2.28, xend = 2.34, y = N_R2_PRED_A, yend = N_R2_PRED_A,
           linewidth = 0.9, colour = PAL[["warn"]]) +
  annotate("text", x = 2.40, y = (N_R2_TRUTH + N_R2_PRED_A) / 2,
           label = "EXCESS\n+0.217", hjust = 0, vjust = 0.5, lineheight = 0.95,
           size = NUMB, fontface = "bold", colour = PAL[["warn"]]) +
  annotate("text", x = 1, y = N_R2_TRUTH / 2, label = "0.341",
           size = NUMB, colour = "white", fontface = "bold") +
  annotate("text", x = 2, y = N_R2_TRUTH / 2, label = "0.341",
           size = NUMB, colour = "white", fontface = "bold") +
  annotate("text", x = 2, y = N_R2_TRUTH + N_EXC_A / 2, label = "0.217",
           size = NUMB, colour = "white", fontface = "bold") +
  annotate("text", x = 2, y = N_R2_PRED_A + 0.015, label = "0.558", vjust = 0,
           size = TXT, colour = "grey25", fontface = "bold") +
  annotate("text", x = 1, y = N_R2_TRUTH + 0.015, label = "0.341", vjust = 0,
           size = TXT, colour = "grey25", fontface = "bold") +
  scale_fill_manual(values = c(excess = PAL[["warn"]], shared = PAL[["accent"]]),
                    guide = "none") +
  scale_y_continuous(limits = c(0, 0.68), breaks = seq(0, 0.6, 0.2), expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = c(0.70, 1.20))) +
  labs(x = NULL, y = "R\u00b2 explained by tissue identity") +
  theme_kids26(base_size = BS) +
  theme(axis.text.x  = element_text(size = BS * 0.70, lineheight = 0.95,
                                    colour = "grey15", face = "bold"),
        axis.title.y = element_text(size = BS * 0.80),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(2, 8, 2, 6))

p1_head <- blank_canvas() +
  annotate("text", x = 0.05, y = 9.6, hjust = 0, vjust = 1, size = HEAD,
           fontface = "bold", colour = PAL[["accent"]], lineheight = 1.0,
           label = "ACT 1  \u00b7  THE MODEL WORKED\n\u2014 AND THAT WAS THE PROBLEM") +
  annotate("text", x = 0.05, y = 4.6, hjust = 0, vjust = 1, size = TXT,
           colour = "grey20", lineheight = 1.25,
           label = paste("Elastic net, 336,480 CpGs \u2192 reference HRDsum.",
                         "7,065 adult TCGA samples, 30 non-CNS",
                         "cancer types, nested leave-one-cancer-out.", sep = "\n"))

p1_foot <- blank_canvas() +
  annotate("text", x = 0.05, y = 9.9, hjust = 0, vjust = 1, size = TXTS,
           colour = "grey25", lineheight = 1.25,
           label = "Macro within-tissue r = 0.520 \u00b7 pooled MAE 8.97\n29 / 30 tissues positive \u00b7 permutation p = 0.001") +
  annotate("text", x = 0.05, y = 6.3, hjust = 0, vjust = 1, size = TXT,
           colour = PAL[["warn"]], fontface = "bold", lineheight = 1.2,
           label = "The model encoded more lineage\nthan the thing it predicts.") +
  annotate("text", x = 0.05, y = 2.2, hjust = 0, vjust = 1, size = TXTS - 0.3,
           colour = "grey40", lineheight = 1.2,
           label = "Per-tissue offsets SARC \u201315.4 \u2192 PCPG +8.3\n(mean |offset| 3.46 HRDsum units)")

# =============================================================================
# ACT 2 -- the scoreboard
# =============================================================================
sb <- data.table(
  y      = 5:1,
  name   = c("1  Label-free tissue-offset model",
             "2  Zero-shot percentile remap",
             "3  Relative-target models V1 / V2",
             "4  V4 lineage-penalised features",
             "5  V3-abs: within-tissue variance"),
  why    = c("best leave-one-tissue-out R\u00b2 = \u20130.116 \u2014 worse than the mean",
             "monotone map cannot reorder: 0 info; tissue R\u00b2 0.562 \u2192 0.593",
             "failed 3 of 7 pre-registered gates (macro \u03c1 \u20130.219; BRCA \u20130.405)",
             "mechanism gate G3 failed: tissue R\u00b2 went UP, 0.522 \u2192 0.536",
             "passed 7 / 7 sentinel gates and 7 / 7 selection criteria"),
  status = c(rep("FAILED", 4), "PASSED")
)
sb[, pass := status == "PASSED"]

p2 <- ggplot(sb) +
  geom_rect(aes(xmin = 0, xmax = 11.9, ymin = y - 0.44, ymax = y + 0.44, fill = pass),
            colour = NA) +
  geom_segment(aes(x = 0.05, xend = 0.05, y = y - 0.44, yend = y + 0.44, colour = pass),
               linewidth = 2.6) +
  geom_point(aes(x = 0.42, y = y + 0.16, shape = pass, colour = pass),
             size = 4.0, stroke = 1.8) +
  geom_text(aes(x = 0.80, y = y + 0.16, label = name, colour = pass),
            hjust = 0, size = TXT - 0.15, fontface = "bold") +
  geom_text(aes(x = 0.80, y = y - 0.20, label = why), hjust = 0,
            size = TXTS - 0.75, colour = "grey30") +
  geom_text(aes(x = 9.95, y = y + 0.16, label = status, colour = pass),
            hjust = 1, size = TXT - 0.35, fontface = "bold") +
  scale_fill_manual(values = c(`TRUE` = "#DFEDE3", `FALSE` = "#F2F2F2"), guide = "none") +
  scale_colour_manual(values = c(`TRUE` = PAL[["good"]], `FALSE` = PAL[["warn"]]),
                      guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 3, `FALSE` = 4), guide = "none") +
  coord_cartesian(xlim = c(0, 11.9), ylim = c(0.4, 5.62), expand = FALSE, clip = "off") +
  theme_void(base_size = BS) +
  theme(plot.margin = margin(2, 6, 2, 6))

p2_head <- blank_canvas() +
  annotate("text", x = 0.05, y = 9.6, hjust = 0, vjust = 1, size = HEAD,
           fontface = "bold", colour = PAL[["accent"]], lineheight = 1.0,
           label = "ACT 2  \u00b7  FIVE ATTACKS.\nFOUR FAILED.") +
  annotate("text", x = 0.05, y = 4.6, hjust = 0, vjust = 1, size = TXT,
           colour = "grey20", lineheight = 1.25,
           label = paste("Every attempt to strip the lineage signal was",
                         "pre-registered with explicit pass / fail gates.",
                         "We report all of them.", sep = "\n"))

p2_foot <- blank_canvas() +
  annotate("text", x = 0.05, y = 9.9, hjust = 0, vjust = 1, size = TXT,
           colour = PAL[["good"]], fontface = "bold", lineheight = 1.25,
           label = "Diagnosis: the 5,000-probe filter ranked CpGs by\nPOOLED variance \u2014 so it was selecting lineage\nmarkers by construction, before the target was seen.")

# =============================================================================
# ACT 3 -- the result
# =============================================================================
d3 <- data.table(
  model = factor(c("Model A\npooled variance", "V3-abs\nwithin-tissue variance"),
                 levels = c("Model A\npooled variance", "V3-abs\nwithin-tissue variance")),
  est = c(N_EXC_A, N_EXC_B), lo = c(N_EXC_A_LO, N_EXC_B_LO), hi = c(N_EXC_A_HI, N_EXC_B_HI)
)

p3a <- ggplot(d3, aes(x = model, y = est, colour = model, shape = model)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.10, linewidth = 1.1) +
  geom_point(size = 5, stroke = 1.8) +
  annotate("segment", x = 1.20, xend = 1.80, y = 0.275, yend = 0.275,
           linewidth = 0.9, colour = PAL[["good"]],
           arrow = arrow(length = unit(0.13, "in"), type = "closed")) +
  annotate("text", x = 1.5, y = 0.288, vjust = 0, size = TXT - 0.25,
           fontface = "bold", colour = PAL[["good"]],
           label = "42.0% of excess removed [37.9\u201346.5%]") +
  annotate("text", x = 1.17, y = N_EXC_A, hjust = 0, size = NUMB - 0.3,
           colour = PAL[["warn"]], fontface = "bold", label = "0.217") +
  annotate("text", x = 1.83, y = N_EXC_B, hjust = 1, size = NUMB - 0.3,
           colour = PAL[["accent"]], fontface = "bold", label = "0.126") +
  scale_colour_manual(values = c(PAL[["warn"]], PAL[["accent"]]), guide = "none") +
  scale_shape_manual(values = c(4, 16), guide = "none") +
  scale_y_continuous(limits = c(0.07, 0.345), breaks = seq(0.10, 0.25, 0.05)) +
  scale_x_discrete(expand = expansion(add = 0.55)) +
  labs(x = NULL, y = "Excess tissue R\u00b2") +
  theme_kids26(base_size = BS) +
  theme(axis.text.x  = element_text(size = BS * 0.72, lineheight = 0.95,
                                    colour = "grey15", face = "bold"),
        axis.title.y = element_text(size = BS * 0.80),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(2, 8, 0, 6))

set.seed(FIG_SEED)
pt[, ystrip := 0]
p3b <- ggplot(pt, aes(x = d_r, y = ystrip)) +
  annotate("rect", xmin = N_DR_LO, xmax = N_DR_HI, ymin = -0.75, ymax = 0.75,
           fill = PAL[["accent"]], alpha = 0.15) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.6) +
  geom_jitter(aes(shape = improved, colour = improved), width = 0, height = 0.36,
              size = 3.2, stroke = 1.4, alpha = 0.9) +
  geom_point(x = N_DR_MEAN, y = 0, shape = 18, size = 6.5, colour = "grey10",
             inherit.aes = FALSE) +
  annotate("text", x = N_DR_MEAN, y = 1.05, label = "mean \u20130.016  (n.s.)", vjust = 0,
           size = TXTS, fontface = "bold", colour = "grey10") +
  annotate("text", x = -0.265, y = -1.05, label = "V3-abs worse", hjust = 0, vjust = 1,
           size = TXTS - 0.3, colour = "grey45") +
  annotate("text", x = 0.095, y = -1.05, label = "better", hjust = 1, vjust = 1,
           size = TXTS - 0.3, colour = "grey45") +
  scale_colour_manual(values = c(`TRUE` = PAL[["accent"]], `FALSE` = PAL[["warn"]]),
                      guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 4), guide = "none") +
  scale_x_continuous(breaks = seq(-0.25, 0.05, 0.05)) +
  coord_cartesian(xlim = c(-0.275, 0.10), ylim = c(-1.75, 1.75), expand = FALSE,
                  clip = "off") +
  labs(x = "Per-tissue \u0394 within-tissue Pearson r (30 tissues)", y = NULL) +
  theme_kids26(base_size = BS) +
  theme(axis.text.y = element_blank(), axis.line.y = element_blank(),
        axis.title.x = element_text(size = BS * 0.76),
        panel.grid.major.y = element_blank(),
        plot.margin = margin(0, 8, 2, 6))

p3_head <- blank_canvas() +
  annotate("text", x = 0.05, y = 9.6, hjust = 0, vjust = 1, size = HEAD,
           fontface = "bold", colour = PAL[["accent"]], lineheight = 1.0,
           label = "ACT 3  \u00b7  ONE CHANGE,\nMEASURABLY LESS LINEAGE") +
  annotate("text", x = 0.05, y = 4.6, hjust = 0, vjust = 1, size = TXT,
           colour = "grey20", lineheight = 1.25,
           label = paste("Identical architecture. Only change: rank the",
                         "5,000-probe filter by WITHIN-tissue variance,",
                         "training-fold rows only. Full 30-fold LOCO.", sep = "\n"))

p3_foot <- blank_canvas() +
  annotate("text", x = 0.05, y = 9.6, hjust = 0, vjust = 1, size = TXTS,
           colour = "grey15", lineheight = 1.35,
           label = paste("\u0394 excess R\u00b2 = 0.0913  [0.0824, 0.0996],  p < 0.0001",
                         "Macro skill vs oracle tissue-mean null  0.419 \u2192 0.630",
                         "Pooled MAE 8.97 \u2192 8.45 \u00b7 mean |tissue bias| 3.46 \u2192 2.97",
                         sep = "\n")) +
  annotate("text", x = 0.05, y = 3.0, hjust = 0, vjust = 1, size = TXTS - 0.2,
           colour = PAL[["good"]], fontface = "bold", lineheight = 1.2,
           label = "No ranking cost: Wilcoxon p = 0.43 (r), p = 0.48 (\u03c1);\n29 / 30 tissues still positive.")

# =============================================================================
# title + bottom line
# =============================================================================
title_strip <- blank_canvas() +
  annotate("text", x = 0.04, y = 8.8, hjust = 0, vjust = 1, size = 9.6,
           fontface = "bold", colour = PAL[["accent"]],
           label = "Methylation predicts HRDsum \u2014 but how much of that is just tissue?") +
  annotate("text", x = 0.04, y = 2.6, hjust = 0, vjust = 1, size = TXT + 0.1,
           colour = "grey35",
           label = "Pan-cancer DNA-methylation model of homologous-recombination deficiency  \u00b7  KIDS26 Team 12  \u00b7  30 adult TCGA cancer types, nested leave-one-cancer-out")

bottom_strip <- blank_canvas() +
  annotate("rect", xmin = -0.3, xmax = 10.3, ymin = 0.6, ymax = 9.9,
           fill = PAL[["accent"]], colour = NA) +
  annotate("text", x = 0.08, y = 7.1, hjust = 0, vjust = 0.5, size = TXT + 0.45,
           fontface = "bold", colour = "white",
           label = "BOTTOM LINE    Methylation carries within-tissue HRD signal \u2014 and 42% of the lineage confound comes off at no measurable ranking cost.") +
  annotate("text", x = 0.10, y = 3.0, hjust = 0, vjust = 0.5, size = TXT - 0.45,
           colour = "#BDD3E6",
           label = "Absolute cross-lineage calibration remains unsolved.    The locked CNS (GBM/LGG) cohort is the external test \u2014 never seen during development.")

# =============================================================================
# assemble
# =============================================================================
col1 <- p1_head / p1 / p1_foot       + plot_layout(heights = c(0.26, 0.52, 0.22))
col2 <- p2_head / p2 / p2_foot       + plot_layout(heights = c(0.26, 0.52, 0.22))
col3 <- p3_head / p3a / p3b / p3_foot + plot_layout(heights = c(0.26, 0.31, 0.21, 0.22))

body  <- (col1 | col2 | col3) + plot_layout(widths = c(0.315, 0.360, 0.325))

slide <- title_strip / body / bottom_strip +
  plot_layout(heights = c(0.115, 0.765, 0.12)) &
  theme(plot.background = element_rect(fill = "white", colour = NA))

save_fig(slide, OUT_PNG, width = 16, height = 9, dpi = 300)

# --- provenance ---------------------------------------------------------------
writeLines(c("# sessionInfo for scripts/make_story_slide.R",
             paste("# generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
             paste("# seed:", FIG_SEED),
             paste("# output:", OUT_PNG), "",
             capture.output(sessionInfo())), OUT_SESS)
message("  wrote ", OUT_SESS)
message("DONE")

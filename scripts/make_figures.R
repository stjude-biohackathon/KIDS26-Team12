# =============================================================================
# scripts/make_figures.R -- presentation figures for KIDS26-Team12
# -----------------------------------------------------------------------------
# Run from the repository root:   Rscript scripts/make_figures.R
#
# INPUTS (all required; the script fails loudly if any is missing)
#   results/loco_run01/loco_predictions.tsv             7,065 x 12, source LOCO
#   results/loco_run01/loco_metrics.tsv                 per-tissue LOCO metrics
#   results/loco_run01/per_tissue_within_correlation.tsv  cancer_type,n,r,rho,bias
#   results/loco_run01/pooled_within_tissue_metrics.tsv
#   results/c1_shrinkage/shrinkage_macro_by_k.tsv       EB few-shot calibration
#   results/components_2026-09-18/04_pred_component_macro.tsv
#   results/components_2026-09-18/04b_pred_component_paired.tsv
#   R/figures.R                                         shared theme + helpers
#
# OUTPUTS (results/figures/)
#   fig01_workflow.mmd              mermaid source (ggplot is the wrong tool here)
#   fig02_pred_vs_obs.png           predicted vs observed HRDsum, hexbin
#   fig03_per_tissue_correlation.png  within-tissue Pearson + Spearman, 30 tissues
#   fig04_tissue_bias.png           per-tissue signed offset (the C1 failure mode)
#   fig05_lineage_imprinting.png    tissue R2 of prediction vs of truth
#   fig06_fewshot_calibration.png   fraction of oracle gain recovered vs k
#   fig07_components.png            LOH / LST / TAI / HRDsum tracking
#   sessionInfo.txt                 session provenance
#
# CONVENTIONS
#   * Adopted C2 rule applied before ANY computation: pred_clip = pmax(pred, 0).
#   * CNS LOCK: development cancers only. GBM/LGG are filtered explicitly and the
#     script aborts if either appears. No CNS HRDsum value is read or plotted.
#   * No LSF submission, no access to data/processed/beta.tsv.
#   * Deterministic: fixed seed (FIG_SEED), no stochastic steps in the figures.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

## ---- setup ------------------------------------------------------------------
ROOT <- normalizePath(".", mustWork = TRUE)
if (!file.exists(file.path(ROOT, "R", "figures.R"))) {
  stop("Run this script from the repository root (R/figures.R not found at ", ROOT, ")",
       call. = FALSE)
}
source(file.path(ROOT, "R", "figures.R"))
set.seed(FIG_SEED)

FIG_DIR <- file.path(ROOT, "results", "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

message("== make_figures.R ==  seed = ", FIG_SEED)
message("-- reading inputs")

P_PRED    <- file.path(ROOT, "results/loco_run01/loco_predictions.tsv")
P_MET     <- file.path(ROOT, "results/loco_run01/loco_metrics.tsv")
P_WITHIN  <- file.path(ROOT, "results/loco_run01/per_tissue_within_correlation.tsv")
P_POOLED  <- file.path(ROOT, "results/loco_run01/pooled_within_tissue_metrics.tsv")
P_SHRINK  <- file.path(ROOT, "results/c1_shrinkage/shrinkage_macro_by_k.tsv")
P_CMACRO  <- file.path(ROOT, "results/components_2026-09-18/04_pred_component_macro.tsv")
P_CPAIR   <- file.path(ROOT, "results/components_2026-09-18/04b_pred_component_paired.tsv")

pred     <- read_required(P_PRED)
met      <- read_required(P_MET)
within   <- read_required(P_WITHIN)
pooled   <- read_required(P_POOLED)
shrink   <- read_required(P_SHRINK)
cmacro   <- read_required(P_CMACRO)
cpair    <- read_required(P_CPAIR)

## ---- cohort lock + C2 rule ---------------------------------------------------
pred   <- keep_development(pred,   "cancer_type", "loco_predictions")
met    <- keep_development(met,    "cancer_type", "loco_metrics")
within <- keep_development(within, "cancer_type", "per_tissue_within_correlation")
pred   <- apply_c2_clip(pred)

stopifnot(nrow(pred) == 7065L, uniqueN(pred$cancer_type) == 30L,
          nrow(within) == 30L, nrow(met) == 30L)

# The shipped per_tissue_within_correlation.tsv is computed on the UNCLIPPED
# prediction. The adopted C2 rule must be applied before anything is computed, so
# the per-tissue statistics used in fig03/fig04 are recomputed here from pred_clip.
# Clipping ties all negative predictions at zero, which moves both the correlation
# and the offset in tissues with many negatives (KICH: r 0.803 -> 0.888, offset
# -0.787 -> +1.088). The shipped file is retained only as an integrity check on the
# unclipped statistics.
within_shipped <- copy(within)
within <- pred[, .(n    = .N,
                   r    = cor(pred_clip, actual),
                   rho  = cor(pred_clip, actual, method = "spearman"),
                   bias = mean(pred_clip - actual)), by = cancer_type]
raw_chk <- pred[, .(r_raw = cor(predicted_reference_HRDsum, actual)), by = cancer_type]
int_chk <- merge(raw_chk, within_shipped[, .(cancer_type, r)], by = "cancer_type")
if (max(abs(int_chk$r_raw - int_chk$r)) > 1e-10) {
  stop("Integrity check failed: unclipped per-tissue Pearson does not reproduce ",
       "results/loco_run01/per_tissue_within_correlation.tsv", call. = FALSE)
}
message(sprintf("-- per-tissue stats recomputed under C2 clipping (max |dr| vs shipped = %.4f)",
                max(abs(within$r - within_shipped$r[match(within$cancer_type,
                                                          within_shipped$cancer_type)]))))

# Within-tissue pooled correlation, also recomputed under clipping.
pred[, `:=`(pc_c = pred_clip - mean(pred_clip), ac_c = actual - mean(actual)),
     by = cancer_type]
w_r_clip   <- cor(pred$pc_c, pred$ac_c)
w_rho_clip <- cor(rank(pred$pc_c), rank(pred$ac_c))
macro_r    <- mean(within$r)
message(sprintf("-- cohort: %d specimens, %d development cancer types (GBM/LGG absent)",
                nrow(pred), uniqueN(pred$cancer_type)))

## =============================================================================
## fig01 -- workflow (mermaid; a node/edge diagram is awkward and fragile in ggplot)
## =============================================================================
mmd <- c(
  "%% fig01_workflow.mmd -- KIDS26-Team12 evaluation design",
  "%% Render: mermaid-cli  ->  mmdc -i fig01_workflow.mmd -o fig01_workflow.png -s 3",
  "flowchart TD",
  "    A[\"TCGA non-CNS development cohort<br/>7,065 specimens &middot; 30 cancer types\"]",
  "    B[\"Nested leave-one-cancer-out (LOCO)<br/>outer loop: hold out 1 cancer type<br/>inner folds: one per training cancer type\"]",
  "    C[\"Candidate selection<br/>inner-fold metrics only<br/>held-out cancer never seen\"]",
  "    D[\"FREEZE<br/>single candidate, hyperparameters<br/>and decision rules locked\"]",
  "    E[\"LOCKED CNS cohort<br/>642 GBM + LGG<br/>external test, opened ONCE\"]",
  "    A --> B --> C --> D --> E",
  "    D -. \"no refitting, no re-selection<br/>after the freeze\" .-> E",
  "    classDef dev  fill:#e7eff7,stroke:#1F4E79,stroke-width:2px,color:#102a43;",
  "    classDef frz  fill:#d1e7dd,stroke:#0F5132,stroke-width:2px,color:#0F5132;",
  "    classDef lock fill:#f8d7da,stroke:#B02A37,stroke-width:3px,color:#7a1420;",
  "    class A,B,C dev;",
  "    class D frz;",
  "    class E lock;"
)
writeLines(mmd, file.path(FIG_DIR, "fig01_workflow.mmd"))
message("  wrote results/figures/fig01_workflow.mmd (ggplot skipped: node/edge diagram)")

## =============================================================================
## fig02 -- predicted vs observed HRDsum (source LOCO, 7,065 specimens)
## =============================================================================
pool_r   <- cor(pred$pred_clip, pred$actual)
pool_rho <- cor(pred$pred_clip, pred$actual, method = "spearman")
w_r      <- w_r_clip
w_rho    <- w_rho_clip
ax_max   <- ceiling(max(pred$actual, pred$pred_clip) / 10) * 10

fig02 <- ggplot(pred, aes(actual, pred_clip)) +
  geom_abline(slope = 1, intercept = 0, colour = PAL[["warn"]],
              linewidth = 0.7, linetype = "22") +
  geom_hex(bins = 55) +
  geom_point(alpha = 0.06, size = 0.35, colour = "grey15") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = PAL[["accent"]], linewidth = 0.8) +
  scale_fill_gradient(low = "#DCE6F1", high = PAL[["accent"]],
                      trans = "log10", name = "specimens per hex",
                      breaks = c(1, 10, 100)) +
  coord_fixed(xlim = c(0, ax_max), ylim = c(0, ax_max), expand = FALSE) +
  labs(
    title = "Predicted vs observed HRDsum on held-out cancer types",
    subtitle = sprintf(
      paste0("Leave-one-cancer-out, 7,065 TCGA specimens, 30 development cancer types.\n",
             "Pooled r = %s, rho = %s;  within-tissue r = %s;  macro mean of the 30 per-tissue r = %s."),
      fmt(pool_r, 3), fmt(pool_rho, 3), fmt(w_r, 3), fmt(macro_r, 3)),
    x = "Observed HRDsum (LOH + LST + TAI)",
    y = "Predicted HRDsum (clipped at 0)",
    caption = paste0(
      "Red dashed = identity; blue = OLS fit. Every point is a specimen whose cancer type was absent from training.\n",
      "Pooled correlation is inflated by between-tissue mean differences; the within-tissue value is the honest one.")
  ) +
  theme_kids26(base_size = 12) +
  theme(legend.key.width = grid::unit(1.6, "lines"))
save_fig(fig02, file.path(FIG_DIR, "fig02_pred_vs_obs.png"), width = 8.2, height = 8.6)

## =============================================================================
## fig03 -- per-tissue within-tissue Pearson and Spearman
## =============================================================================
w <- copy(within)
setorder(w, r)
w[, flagged := cancer_type %chin% NEAR_ZERO_TISSUES]
w[, lab := sprintf("%s  (n = %d)", cancer_type, n)]
w[, lab := factor(lab, levels = lab)]
wl <- melt(w, id.vars = c("lab", "n", "flagged"),
           measure.vars = c("r", "rho"), variable.name = "stat", value.name = "value")
wl[, stat := factor(stat, levels = c("r", "rho"),
                    labels = c("Pearson r", "Spearman rho"))]
seg <- w[, .(lab, flagged, lo = pmin(r, rho), hi = pmax(r, rho))]
flag_idx <- which(w$flagged)

fig03 <- ggplot() +
  geom_rect(data = w[flagged == TRUE], aes(xmin = as.numeric(lab) - 0.5,
                                           xmax = as.numeric(lab) + 0.5),
            ymin = -0.25, ymax = 0.92, fill = PAL[["warn"]], alpha = 0.09,
            inherit.aes = FALSE) +
  geom_text(data = w[flagged == TRUE][which.max(as.numeric(lab))],
            aes(x = lab, y = 0.90), hjust = 1, vjust = -0.6,
            label = "near-zero within-tissue signal", colour = PAL[["warn"]],
            size = 3.3, fontface = "bold", inherit.aes = FALSE) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.5) +
  geom_segment(data = seg, aes(x = lab, xend = lab, y = lo, yend = hi),
               colour = "grey75", linewidth = 0.6) +
  geom_point(data = wl, aes(lab, value, shape = stat, colour = stat),
             size = 2.7, stroke = 0.9) +
  scale_colour_manual(values = c("Pearson r" = PAL[["accent"]],
                                 "Spearman rho" = PAL[["amber"]]), name = NULL) +
  scale_shape_manual(values = c("Pearson r" = 16, "Spearman rho" = 17), name = NULL) +
  scale_y_continuous(limits = c(-0.25, 0.92), breaks = seq(-0.2, 0.8, 0.2),
                     expand = expansion(mult = 0)) +
  coord_flip() +
  labs(
    title = "Within-tissue correlation varies enormously by cancer type",
    subtitle = sprintf(
      paste0("Leave-one-cancer-out run 01, 30 held-out cancer types.\n",
             "Four tissues (%s) carry essentially no within-tissue signal (shaded)."),
      paste(NEAR_ZERO_TISSUES, collapse = ", ")),
    x = NULL, y = "Correlation of predicted with observed HRDsum, within the held-out tissue",
    caption = paste0("Grey bar spans Pearson to Spearman. Sample size shown beside each cancer type.\n",
                     "Computed under the adopted C2 rule pred_clip = pmax(pred, 0).")
  ) +
  theme_kids26(base_size = 12) +
  theme(panel.grid.major.y = element_blank())
save_fig(fig03, file.path(FIG_DIR, "fig03_per_tissue_correlation.png"),
         width = 9.0, height = 8.4)

## =============================================================================
## fig04 -- per-tissue mean signed bias (the C1 offset) -- central failure mode
## =============================================================================
b <- copy(within)
setorder(b, bias)
b[, cancer_type := factor(cancer_type, levels = cancer_type)]
b[, sign_lab := factor(ifelse(bias < 0, "Under-predicted (offset < 0)",
                              "Over-predicted (offset > 0)"),
                       levels = c("Under-predicted (offset < 0)",
                                  "Over-predicted (offset > 0)"))]
mean_abs_offset <- mean(abs(b$bias))
sd_offset       <- sd(b$bias)

fig04 <- ggplot(b, aes(cancer_type, bias, colour = sign_lab)) +
  geom_hline(yintercept = 0, colour = "grey25", linewidth = 0.6) +
  geom_hline(yintercept = c(-mean_abs_offset, mean_abs_offset),
             colour = PAL[["neutral"]], linetype = "22", linewidth = 0.4) +
  geom_segment(aes(xend = cancer_type, y = 0, yend = bias), linewidth = 0.9) +
  geom_point(size = 3) +
  geom_text(aes(label = fmt(bias, 1),
                hjust = ifelse(bias < 0, 1.35, -0.35)),
            size = 2.9, show.legend = FALSE) +
  annotate("label", x = 4.0, y = -17.8,
           label = sprintf("mean |offset| = %s HRD units\nSD %s   range %s to %s",
                           fmt(mean_abs_offset, 2), fmt(sd_offset, 2),
                           fmt(min(b$bias), 2), fmt(max(b$bias), 2)),
           hjust = 0, size = 3.4, colour = "grey15",
           fill = "white", lineheight = 1.1) +
  scale_colour_manual(values = c("Under-predicted (offset < 0)" = PAL[["accent"]],
                                 "Over-predicted (offset > 0)"  = PAL[["warn"]]),
                      name = NULL) +
  scale_y_continuous(limits = c(-18.5, 11), breaks = seq(-15, 10, 5)) +
  coord_flip() +
  labs(
    title = "C1: every cancer type carries its own calibration offset",
    subtitle = paste0("Mean signed error (predicted - observed HRDsum) for each cancer type,\n",
                      "computed when that type was entirely absent from training."),
    x = NULL, y = "Per-tissue mean signed bias (HRD units)",
    caption = paste0(
      "Dashed lines mark +/- the mean absolute offset. This offset is the project's central failure mode: it cannot be\n",
      "estimated without labels from the new tissue, so an absolute HRDsum for an unseen cancer type is not yet trustworthy.\n",
      "Computed under the adopted C2 rule pred_clip = pmax(pred, 0).")
  ) +
  theme_kids26(base_size = 12) +
  theme(panel.grid.major.y = element_blank())
save_fig(fig04, file.path(FIG_DIR, "fig04_tissue_bias.png"), width = 9.0, height = 8.4)

## =============================================================================
## fig05 -- LINEAGE IMPRINTING (key conceptual figure)
##   R2(cancer_type -> pred_clip)  vs  R2(cancer_type -> actual), 30 tissues
## =============================================================================
r2_of <- function(y, g) {
  gm  <- mean(y)
  ssr <- sum((ave_by(y, g) - gm)^2)
  sst <- sum((y - gm)^2)
  ssr / sst
}
ave_by <- function(y, g) {
  m <- tapply(y, g, mean)
  as.numeric(m[as.character(g)])
}

g <- pred$cancer_type
R2_PRED   <- r2_of(pred$pred_clip, g)
R2_ACTUAL <- r2_of(pred$actual,    g)
R2_GAP    <- R2_PRED - R2_ACTUAL

message(sprintf("-- fig05: R2(tissue -> pred_clip) = %.4f ; R2(tissue -> actual) = %.4f ; gap = %.4f",
                R2_PRED, R2_ACTUAL, R2_GAP))

tm <- pred[, .(mean_pred = mean(pred_clip), mean_actual = mean(actual), n = .N),
           by = cancer_type]
grand_actual <- mean(pred$actual)
grand_pred   <- mean(pred$pred_clip)
# compression: how far tissue means sit from their own grand mean, pred vs truth
sd_tm_pred   <- sd(tm$mean_pred)
sd_tm_actual <- sd(tm$mean_actual)

lim <- range(c(tm$mean_pred, tm$mean_actual)); lim <- c(floor(lim[1]) - 2, ceiling(lim[2]) + 3)

pA <- ggplot(tm, aes(mean_actual, mean_pred)) +
  geom_abline(slope = 1, intercept = 0, colour = PAL[["neutral"]],
              linetype = "22", linewidth = 0.6) +
  geom_hline(yintercept = grand_pred,   colour = PAL[["warn"]], linewidth = 0.4, alpha = 0.7) +
  geom_vline(xintercept = grand_actual, colour = PAL[["warn"]], linewidth = 0.4, alpha = 0.7) +
  geom_segment(aes(xend = mean_actual, yend = mean_actual),
               colour = "grey75", linewidth = 0.4) +
  geom_point(aes(size = n), colour = PAL[["accent"]], alpha = 0.85) +
  ggrepel::geom_text_repel(aes(label = cancer_type), size = 2.8, colour = "grey25",
                           seed = FIG_SEED, max.overlaps = 30, min.segment.length = 0.25,
                           segment.colour = "grey70") +
  scale_size_area(max_size = 6, name = "n specimens", breaks = c(50, 200, 700)) +
  coord_fixed(xlim = lim, ylim = lim, expand = FALSE) +
  annotate("text", x = grand_actual + 0.6, y = lim[1] + 1.6,
           label = "grand mean", colour = PAL[["warn"]], size = 2.9, hjust = 0) +
  labs(title = "b. Tissue means compressed toward the grand mean",
       subtitle = sprintf("SD of the 30 tissue means:\ntruth %s vs prediction %s HRD units.",
                          fmt(sd_tm_actual, 2), fmt(sd_tm_pred, 2)),
       x = "Observed mean HRDsum in the tissue",
       y = "Predicted mean HRDsum in the tissue") +
  theme_kids26(base_size = 12)

bars <- data.table(
  what = factor(c("Observed HRDsum\n(the biology)", "Model prediction\n(pred_clip)"),
                levels = c("Observed HRDsum\n(the biology)", "Model prediction\n(pred_clip)")),
  r2   = c(R2_ACTUAL, R2_PRED))

pB <- ggplot(bars, aes(what, r2, fill = what)) +
  geom_col(width = 0.56, show.legend = FALSE) +
  geom_text(aes(label = fmt(r2, 4)), vjust = -0.55, size = 5.2, fontface = "bold",
            colour = "grey15") +
  annotate("segment", x = 1, xend = 2, y = R2_PRED + 0.075, yend = R2_PRED + 0.075,
           colour = "grey30", linewidth = 0.5) +
  annotate("text", x = 1.5, y = R2_PRED + 0.155,
           label = sprintf("excess lineage imprinting\n+%s R2 (%s x)",
                           fmt(R2_GAP, 4), fmt(R2_PRED / R2_ACTUAL, 2)),
           size = 3.6, colour = PAL[["warn"]], fontface = "bold", lineheight = 1.1) +
  scale_fill_manual(values = c(PAL[["neutral"]], PAL[["accent"]])) +
  scale_y_continuous(limits = c(0, 0.78), breaks = seq(0, 0.7, 0.1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "a. Tissue identity explains more of the prediction",
       subtitle = paste0("Variance explained by cancer type alone (one-way ANOVA R2),\n",
                         "7,065 specimens, 30 development cancer types."),
       x = NULL, y = expression("Variance explained by cancer type," ~ R^2)) +
  theme_kids26(base_size = 12)

fig05 <- patchwork::wrap_plots(pB, pA, ncol = 2, widths = c(1, 1.08)) +
  patchwork::plot_annotation(
    title = "Lineage imprinting: the model encodes tissue identity more strongly than HRD biology does",
    subtitle = sprintf(
      "R2(cancer type -> prediction) = %s   vs   R2(cancer type -> observed HRDsum) = %s   (gap +%s)",
      fmt(R2_PRED, 4), fmt(R2_ACTUAL, 4), fmt(R2_GAP, 4)),
    caption = paste0(
      "Left: variance in each quantity attributable to cancer type. Right: each tissue's predicted mean against its observed mean;\n",
      "points lie inside the identity line on both tails, i.e. tissue means are shrunk toward the grand mean while tissue identity itself is\n",
      "over-encoded. This is the mechanism behind the C1 per-tissue offset and the reason an unseen tissue cannot yet be scored absolutely."),
    theme = theme_kids26(base_size = 13))
save_fig(fig05, file.path(FIG_DIR, "fig05_lineage_imprinting.png"), width = 14.0, height = 8.0)

## =============================================================================
## fig06 -- few-shot EB-shrinkage calibration (LABELLED path, not zero-shot)
## =============================================================================
sk <- copy(shrink)
stopifnot(all(c("k", "fraction_of_oracle_gain_recovered",
                "fraction_of_oracle_gain_recovered_unshrunk") %in% names(sk)))
skl <- melt(sk, id.vars = "k",
            measure.vars = c("fraction_of_oracle_gain_recovered",
                             "fraction_of_oracle_gain_recovered_unshrunk"),
            variable.name = "method", value.name = "frac")
skl[, method := factor(method,
      levels = c("fraction_of_oracle_gain_recovered",
                 "fraction_of_oracle_gain_recovered_unshrunk"),
      labels = c("Empirical-Bayes shrunk offset", "Unshrunk sample-mean offset"))]

lab3s  <- sk[k == 3,  fraction_of_oracle_gain_recovered]
lab3u  <- sk[k == 3,  fraction_of_oracle_gain_recovered_unshrunk]
lab10s <- sk[k == 10, fraction_of_oracle_gain_recovered]
lab10u <- sk[k == 10, fraction_of_oracle_gain_recovered_unshrunk]

ann <- data.table(
  k     = c(3,   3,    10,    10),
  frac  = c(lab3s + 0.30, lab3u - 0.30, lab10s + 0.22, lab10u - 0.42),
  lab   = c(sprintf("k = 3 shrunk\n%+.0f%% of oracle gain", 100 * lab3s),
            sprintf("k = 3 unshrunk\n%.0f%% - actively harmful", 100 * lab3u),
            sprintf("k = 10 shrunk\n%+.0f%%", 100 * lab10s),
            sprintf("k = 10 unshrunk\n%+.0f%%", 100 * lab10u)),
  col   = c(PAL[["good"]], PAL[["warn"]], PAL[["good"]], PAL[["warn"]]))

fig06 <- ggplot(skl, aes(k, frac, colour = method, shape = method)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
           fill = PAL[["warn"]], alpha = 0.06) +
  geom_hline(yintercept = 0, colour = PAL[["warn"]], linewidth = 0.5) +
  geom_hline(yintercept = 1, colour = PAL[["neutral"]], linetype = "22", linewidth = 0.5) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3.2) +
  geom_label(data = ann, aes(k, frac, label = lab), inherit.aes = FALSE,
             colour = ann$col, size = 3.1, fontface = "bold", lineheight = 1.05,
             fill = "white", linewidth = 0, label.padding = grid::unit(0.12, "lines")) +
  annotate("label", x = 1.05, y = 1.0, label = "oracle (offset known exactly)",
           hjust = 0, size = 3.1, colour = PAL[["neutral"]], fill = "white",
           linewidth = 0, label.padding = grid::unit(0.1, "lines")) +
  annotate("text", x = 19.8, y = -0.17, label = "worse than no correction at all",
           hjust = 1, size = 3.1, colour = PAL[["warn"]]) +
  scale_colour_manual(values = c("Empirical-Bayes shrunk offset" = PAL[["good"]],
                                 "Unshrunk sample-mean offset"   = PAL[["warn"]]),
                      name = NULL) +
  scale_shape_manual(values = c("Empirical-Bayes shrunk offset" = 16,
                                "Unshrunk sample-mean offset"   = 17), name = NULL) +
  scale_x_continuous(breaks = sk$k, transform = "log10",
                     limits = c(0.92, 22)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(-2, 1, 0.5)) +
  labs(
    title = "Few-shot calibration works, but only with labels and shrinkage",
    subtitle = paste0("Fraction of the achievable (oracle) offset-correction gain recovered from k labelled samples of the ",
                      "new cancer type.\nMacro-averaged over 29 eligible development tissues."),
    x = "k = number of LABELLED samples available from the new cancer type (log scale)",
    y = "Fraction of oracle gain recovered",
    caption = paste0(
      "THIS IS THE LABELLED FEW-SHOT PATH, NOT ZERO-SHOT. It does not enable an n-of-1 prediction for a tumour type with no labels.\n",
      "Empirical-Bayes shrinkage toward the cross-tissue offset prior converts the small-k regime from harmful to useful: at k = 3 the\n",
      "unshrunk sample mean is noisier than the offset it estimates. Source: results/c1_shrinkage/shrinkage_macro_by_k.tsv.")
  ) +
  theme_kids26()
save_fig(fig06, file.path(FIG_DIR, "fig06_fewshot_calibration.png"), width = 10.0, height = 6.8)

## =============================================================================
## fig07 -- which genomic scar component does the prediction track?
## =============================================================================
cm <- copy(cmacro)
cm[, c("stat", "component") := tstrsplit(statistic, "_", fixed = TRUE)]
cm[, component := factor(component, levels = c("LOH", "LST", "TAI", "HRDsum"),
                         labels = c("HRD-LOH", "LST", "TAI", "HRDsum\n(composite)"))]
cm[, stat := factor(stat, levels = c("pearson", "spearman"),
                    labels = c("Pearson r", "Spearman rho"))]
setorder(cm, component, stat)

pC1 <- ggplot(cm, aes(component, macro_mean_over_30_tissues, fill = stat)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66) +
  geom_text(aes(label = fmt(macro_mean_over_30_tissues, 3)),
            position = position_dodge(width = 0.72), vjust = -0.5, size = 3.1,
            colour = "grey20") +
  scale_fill_manual(values = c("Pearson r" = PAL[["accent"]],
                               "Spearman rho" = PAL[["accent2"]]), name = NULL) +
  scale_y_continuous(limits = c(0, 0.63), breaks = seq(0, 0.6, 0.1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "a. Macro within-tissue correlation with each component",
       subtitle = "Mean over the 30 development cancer types (correlation computed inside each tissue, then averaged).",
       x = NULL, y = "Macro within-tissue correlation") +
  theme_kids26(base_size = 12) +
  theme(panel.grid.major.x = element_blank())

cp <- copy(cpair)
cp[, comparison := factor(comparison, levels = rev(c("TAI - LOH", "TAI - LST",
                                                     "TAI - HRDsum", "LST - LOH")))]
cp[, sig := boot_lo > 0 | boot_hi < 0]
pC2 <- ggplot(cp, aes(mean_diff, comparison, colour = sig)) +
  geom_vline(xintercept = 0, colour = PAL[["warn"]], linewidth = 0.5) +
  geom_errorbar(aes(xmin = boot_lo, xmax = boot_hi), orientation = "y",
                width = 0.16, linewidth = 0.8) +
  geom_point(size = 3.2) +
  geom_text(aes(label = sprintf("%s [%s, %s]   %d/%d tissues",
                                fmt(mean_diff, 3), fmt(boot_lo, 3), fmt(boot_hi, 3),
                                n_tissues_favouring_first, n_tissues)),
            vjust = -1.25, hjust = 0.5, size = 3.0, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = PAL[["good"]], `FALSE` = PAL[["neutral"]]),
                      guide = "none") +
  scale_x_continuous(limits = c(-0.09, 0.215), breaks = seq(-0.05, 0.15, 0.05)) +
  labs(title = "b. Paired differences over the same 30 tissues (95% bootstrap CI)",
       subtitle = "Positive = the first component is tracked better. CI excluding 0 shown in green.",
       x = "Difference in macro within-tissue Pearson r", y = NULL) +
  theme_kids26(base_size = 12) +
  theme(panel.grid.major.y = element_blank())

fig07 <- patchwork::wrap_plots(pC1, pC2, ncol = 1, heights = c(1, 0.85)) +
  patchwork::plot_annotation(
    title = "The model tracks TAI better than LOH, LST or HRDsum itself",
    subtitle = paste0("Component decomposition of run-01 LOCO predictions.\n",
                      "HRDsum = HRD-LOH + LST + TAI holds exactly for all 7,065 specimens."),
    caption = paste0(
      "Source: results/components_2026-09-18/04_pred_component_macro.tsv and 04b_pred_component_paired.tsv (seed 7317).\n",
      "LST vs LOH is indistinguishable. All statistics are n = 30 tissues; read the CIs, not the point estimates alone."),
    theme = theme_kids26(base_size = 13))
save_fig(fig07, file.path(FIG_DIR, "fig07_components.png"), width = 10.0, height = 10.0)

## =============================================================================
## provenance
## =============================================================================
si <- file.path(FIG_DIR, "sessionInfo.txt")
con <- file(si, open = "wt")
writeLines(c(
  "KIDS26-Team12 -- results/figures provenance",
  paste0("script          : scripts/make_figures.R"),
  paste0("seed (FIG_SEED) : ", FIG_SEED),
  paste0("R               : ", R.version.string),
  "cohort          : 7065 specimens, 30 development cancer types (GBM/LGG excluded and asserted absent)",
  "rule applied    : C2 adopted clipping, pred_clip = pmax(predicted_reference_HRDsum, 0)",
  sprintf("fig05 R2(tissue -> pred_clip)   = %.6f", R2_PRED),
  sprintf("fig05 R2(tissue -> actual)      = %.6f", R2_ACTUAL),
  sprintf("fig05 gap                       = %.6f", R2_GAP),
  "",
  "---- sessionInfo() ----"), con)
capture.output(sessionInfo(), file = con, append = TRUE)
close(con)
message("  wrote ", si)

message("== make_figures.R complete ==")

#!/usr/bin/env Rscript
#
# Visualize expHRD_EBPlusPlus_scores.tsv (output of run_expHRD_EBPlusPlus.R).
#
# Usage:
#   Rscript plot_expHRD_results.R <expHRD_EBPlusPlus_scores.tsv> <output_dir>
#
# Produces:
#   1. expHRD_distribution.png   - histogram + density of expHRD_exploratory
#   2. expHRD_pos_vs_neg.png     - scatter of positive vs negative ssGSEA, colored by expHRD
#   3. expHRD_ranked.png         - samples ranked by expHRD score (waterfall style)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_expHRD_results.R <scores.tsv> <output_dir>")
}
input_file <- args[1]
output_dir <- args[2]
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

d <- fread(input_file)
required_cols <- c("sample_id", "positive_ssGSEA", "negative_ssGSEA", "expHRD_exploratory")
missing_cols <- setdiff(required_cols, names(d))
if (length(missing_cols)) {
  stop("Missing expected columns: ", paste(missing_cols, collapse = ", "))
}

cat("Loaded", nrow(d), "samples.\n")
cat("expHRD_exploratory summary:\n")
print(summary(d$expHRD_exploratory))

theme_set(theme_minimal(base_size = 13))

# ---- 1. Distribution: histogram + density overlay ----
p1 <- ggplot(d, aes(x = expHRD_exploratory)) +
  geom_histogram(aes(y = after_stat(density)), bins = 60,
                  fill = "#4C72B0", color = "white", alpha = 0.85) +
  geom_density(color = "#C44E52", linewidth = 1) +
  geom_vline(xintercept = median(d$expHRD_exploratory, na.rm = TRUE),
             linetype = "dashed", color = "black", linewidth = 0.6) +
  labs(
    title = "Distribution of Exploratory expHRD Scores",
    subtitle = paste0("n = ", nrow(d), " TCGA samples (EBPlusPlus, uncalibrated)"),
    x = "expHRD_exploratory (positive_ssGSEA \u2212 negative_ssGSEA)",
    y = "Density"
  ) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(output_dir, "expHRD_distribution.png"), p1, width = 8, height = 5.5, dpi = 300)

# ---- 2. Scatter: positive vs negative ssGSEA, colored by expHRD score ----
p2 <- ggplot(d, aes(x = negative_ssGSEA, y = positive_ssGSEA, color = expHRD_exploratory)) +
  geom_point(alpha = 0.5, size = 1.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_gradient2(
    low = "#4575B4", mid = "grey90", high = "#D73027",
    midpoint = median(d$expHRD_exploratory, na.rm = TRUE),
    name = "expHRD"
  ) +
  labs(
    title = "HRD-Positive vs HRD-Negative ssGSEA Scores",
    subtitle = "Points above the diagonal have higher HRD-positive signature enrichment",
    x = "negative_ssGSEA (HRD_GENE_NEGATIVE)",
    y = "positive_ssGSEA (HRD_GENE_POSITIVE)"
  ) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(output_dir, "expHRD_pos_vs_neg.png"), p2, width = 7.5, height = 6.5, dpi = 300)

# ---- 3. Ranked waterfall plot of all samples ----
d_ranked <- d[order(d$expHRD_exploratory), ]
d_ranked$rank <- seq_len(nrow(d_ranked))

p3 <- ggplot(d_ranked, aes(x = rank, y = expHRD_exploratory, fill = expHRD_exploratory > 0)) +
  geom_col(width = 1) +
  scale_fill_manual(values = c("TRUE" = "#D73027", "FALSE" = "#4575B4"),
                     labels = c("TRUE" = "HRD-high (>0)", "FALSE" = "HRD-low (\u22640)"),
                     name = NULL) +
  labs(
    title = "Samples Ranked by Exploratory expHRD Score",
    subtitle = paste0(sum(d_ranked$expHRD_exploratory > 0), " / ", nrow(d_ranked),
                       " samples above 0 (not the paper's calibrated threshold)"),
    x = "Sample rank (low \u2192 high expHRD)",
    y = "expHRD_exploratory"
  ) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")

ggsave(file.path(output_dir, "expHRD_ranked.png"), p3, width = 9, height = 5.5, dpi = 300)

cat("\nSaved plots to:", normalizePath(output_dir), "\n")
cat(" - expHRD_distribution.png\n - expHRD_pos_vs_neg.png\n - expHRD_ranked.png\n")

#!/usr/bin/env Rscript
#
# Merge expHRD_EBPlusPlus_scores.tsv with master_samples.tsv (cancer type + HRDsum
# ground truth), then plot:
#   1. expHRD_by_cancertype_boxplot.png - box plot of expHRD_exploratory per TCGA cancer type
#   2. expHRD_vs_HRDsum_scatter.png     - scatter of expHRD_exploratory vs HRDsum (like
#                                          expHRD-vs-scarHRD panel in the paper), with
#                                          Pearson correlation annotated
#
# Merge key: patient_id (12-char TCGA barcode), NOT sample_id, because the RNA-seq
# aliquot barcode differs from the methylation aliquot barcode for the same patient.
#
# Usage:
#   Rscript merge_and_plot_expHRD.R <expHRD_EBPlusPlus_scores.tsv> <master_samples.tsv> <output_dir>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript merge_and_plot_expHRD.R <expHRD_scores.tsv> <master_samples.tsv> <output_dir>")
}
exphrd_file  <- args[1]
master_file  <- args[2]
output_dir   <- args[3]
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

exphrd <- fread(exphrd_file)
master <- fread(master_file)

stopifnot("patient_id" %in% names(exphrd), "patient_id" %in% names(exphrd))
stopifnot(all(c("patient_id", "cancer_type", "HRDsum") %in% names(master)))

cat("expHRD rows:", nrow(exphrd), " | unique patients:", uniqueN(exphrd$patient_id), "\n")
cat("master_samples rows:", nrow(master), " | unique patients:", uniqueN(master$patient_id), "\n")

# master_samples.tsv may have >1 row per patient (e.g. multiple methylation files);
# collapse to one row per patient, preferring the first non-NA cancer_type / HRDsum.
master_dedup <- master[
  , .(cancer_type = cancer_type[which(!is.na(cancer_type))[1]],
      HRDsum      = HRDsum[which(!is.na(HRDsum))[1]],
      HRD_LOHLST  = if ("HRD_LOHLST" %in% names(master)) HRD_LOHLST[which(!is.na(HRD_LOHLST))[1]] else NA,
      TAI         = if ("TAI" %in% names(master)) TAI[which(!is.na(TAI))[1]] else NA),
  by = patient_id
]

merged <- merge(exphrd, master_dedup, by = "patient_id", all.x = FALSE)
cat("Merged rows (patients with both expHRD + master_samples data):", nrow(merged), "\n")

fwrite(merged, file.path(output_dir, "merged_expHRD_master.tsv"), sep = "\t")

theme_set(theme_minimal(base_size = 12))

# ---- 1. Box plot: expHRD_exploratory per cancer type ----
type_order <- merged[, .(med = median(expHRD_exploratory, na.rm = TRUE)), by = cancer_type][order(med)]$cancer_type
merged$cancer_type <- factor(merged$cancer_type, levels = type_order)

n_types <- uniqueN(merged$cancer_type)
p1 <- ggplot(merged, aes(x = cancer_type, y = expHRD_exploratory, fill = cancer_type)) +
  geom_boxplot(outlier.size = 0.6, outlier.alpha = 0.4, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    title = "Exploratory expHRD Score by TCGA Cancer Type",
    subtitle = paste0(nrow(merged), " samples across ", n_types, " cancer types (EBPlusPlus, uncalibrated)"),
    x = "TCGA cancer type", y = "expHRD_exploratory"
  ) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
        plot.title = element_text(face = "bold"))

ggsave(file.path(output_dir, "expHRD_by_cancertype_boxplot.png"), p1,
       width = max(9, n_types * 0.35), height = 6, dpi = 300)

# ---- 2. Scatter: expHRD_exploratory vs HRDsum (ground truth), like paper panel (a) ----
merged_valid <- merged[!is.na(HRDsum) & !is.na(expHRD_exploratory)]
ct <- cor.test(merged_valid$HRDsum, merged_valid$expHRD_exploratory, method = "pearson")
pcc <- round(ct$estimate, 3)
pval <- signif(ct$p.value, 3)

p2 <- ggplot(merged_valid, aes(x = HRDsum, y = expHRD_exploratory)) +
  geom_point(alpha = 0.35, size = 1, color = "#4C72B0") +
  geom_smooth(method = "lm", color = "#C44E52", fill = "grey70") +
  annotate("text", x = min(merged_valid$HRDsum, na.rm = TRUE),
           y = max(merged_valid$expHRD_exploratory, na.rm = TRUE),
           hjust = 0, vjust = 1,
           label = paste0("PCC: ", pcc, "\nP = ", pval, "\nn = ", nrow(merged_valid)),
           size = 4) +
  labs(
    title = "expHRD (exploratory) vs Ground-Truth HRDsum",
    x = "HRDsum (ground truth, master_samples.tsv)",
    y = "expHRD_exploratory (EBPlusPlus, uncalibrated)"
  ) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(output_dir, "expHRD_vs_HRDsum_scatter.png"), p2, width = 7.5, height = 6.5, dpi = 300)

cat("\nPearson correlation (expHRD_exploratory vs HRDsum): r =", pcc, ", P =", pval, ", n =", nrow(merged_valid), "\n")
cat("Saved plots to:", normalizePath(output_dir), "\n")
cat(" - expHRD_by_cancertype_boxplot.png\n - expHRD_vs_HRDsum_scatter.png\n - merged_expHRD_master.tsv\n")

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

stopifnot("patient_id" %in% names(exphrd), "expHRD_exploratory" %in% names(exphrd))
stopifnot(all(c("patient_id", "cancer_type", "HRDsum") %in% names(master)))

cat("expHRD rows:", nrow(exphrd), " | unique patients:", uniqueN(exphrd$patient_id), "\n")
cat("master_samples rows:", nrow(master), " | unique patients:", uniqueN(master$patient_id), "\n")

# Restrict to primary tumor specimens (TCGA sample-type code 01, barcode chars 14-15)
# so that patient-level HRDsum labels are never matched to normal, recurrent or
# metastatic RNA. Multiple primary-tumor rows for the same patient (replicates or
# repeat vials) are then deterministically adjudicated by taking the lowest
# sample_id in sort order, rather than picking by any outcome correlation.
stopifnot("sample_id" %in% names(exphrd))
exphrd[, sample_type := substr(sample_id, 14, 15)]
exphrd_primary <- exphrd[sample_type == "01"]
cat("expHRD primary-tumor (code 01) rows:", nrow(exphrd_primary),
    " | unique patients:", uniqueN(exphrd_primary$patient_id), "\n")
setorder(exphrd_primary, patient_id, sample_id)
exphrd_dedup <- exphrd_primary[, .SD[1], by = patient_id]
exphrd_dedup[, sample_type := NULL]

# master_samples.tsv may have >1 row per patient (e.g. multiple methylation files);
# collapse to one row per patient, preferring the first non-NA cancer_type / HRDsum.
master_dedup <- master[
  , .(cancer_type = cancer_type[which(!is.na(cancer_type))[1]],
      HRDsum      = HRDsum[which(!is.na(HRDsum))[1]],
      HRD_LOH     = if ("HRD_LOH" %in% names(master)) HRD_LOH[which(!is.na(HRD_LOH))[1]] else NA,
      LST         = if ("LST" %in% names(master)) LST[which(!is.na(LST))[1]] else NA,
      TAI         = if ("TAI" %in% names(master)) TAI[which(!is.na(TAI))[1]] else NA),
  by = patient_id
]

merged <- merge(exphrd_dedup, master_dedup, by = "patient_id", all.x = FALSE)
cat("Merged rows (primary-tumor patients with both expHRD + master_samples data):", nrow(merged), "\n")

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

# ---- 2. Scatter: expHRD_exploratory vs HRDsum (reference), like paper panel (a) ----
merged_valid <- merged[!is.na(HRDsum) & !is.na(expHRD_exploratory)]

# Correlate two vectors, returning NA (rather than erroring) when there are too
# few finite pairs for a Pearson test (e.g. only 1-2 cancer types present).
safe_cor_test <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(list(estimate = NA_real_, p.value = NA_real_))
  ct <- cor.test(x[ok], y[ok], method = "pearson")
  list(estimate = ct$estimate, p.value = ct$p.value)
}

# Pooled pan-cancer correlation mixes between-cancer-type differences (E1) with the
# within-cancer-type association (E2); report both separately rather than treating
# the pooled estimate as the sole association (see docs/22_CONFOUNDING_AND_PANCANCER_PLAN.md).
ct_pooled <- safe_cor_test(merged_valid$HRDsum, merged_valid$expHRD_exploratory)
pcc_pooled <- round(ct_pooled$estimate, 3)
pval_pooled <- signif(ct_pooled$p.value, 3)

# E2 (within-cancer-type): center both variables on their cancer-type mean, then
# correlate the residuals. This is the clinically relevant "same tissue" estimand.
merged_valid[, `:=`(
  HRDsum_within  = HRDsum - mean(HRDsum, na.rm = TRUE),
  expHRD_within  = expHRD_exploratory - mean(expHRD_exploratory, na.rm = TRUE)
), by = cancer_type]
ct_within <- safe_cor_test(merged_valid$HRDsum_within, merged_valid$expHRD_within)
pcc_within <- round(ct_within$estimate, 3)
pval_within <- signif(ct_within$p.value, 3)

# E1 (between-cancer-type): correlate the cancer-type-level means.
type_means <- merged_valid[, .(HRDsum_mean = mean(HRDsum, na.rm = TRUE),
                                expHRD_mean = mean(expHRD_exploratory, na.rm = TRUE)),
                            by = cancer_type]
ct_between <- safe_cor_test(type_means$HRDsum_mean, type_means$expHRD_mean)
pcc_between <- round(ct_between$estimate, 3)
pval_between <- signif(ct_between$p.value, 3)

p2 <- ggplot(merged_valid, aes(x = HRDsum, y = expHRD_exploratory)) +
  geom_point(alpha = 0.35, size = 1, color = "#4C72B0") +
  geom_smooth(method = "lm", color = "#C44E52", fill = "grey70") +
  annotate("text", x = min(merged_valid$HRDsum, na.rm = TRUE),
           y = max(merged_valid$expHRD_exploratory, na.rm = TRUE),
           hjust = 0, vjust = 1,
           label = paste0(
             "Pooled PCC: ", pcc_pooled, ", P = ", pval_pooled, "\n",
             "Within-type PCC (E2): ", pcc_within, ", P = ", pval_within, "\n",
             "Between-type PCC (E1): ", pcc_between, ", P = ", pval_between, "\n",
             "n = ", nrow(merged_valid)
           ),
           size = 3.5) +
  labs(
    title = "expHRD (exploratory) vs Reference HRDsum",
    x = "HRDsum (reference, master_samples.tsv)",
    y = "expHRD_exploratory (EBPlusPlus, uncalibrated)"
  ) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(output_dir, "expHRD_vs_HRDsum_scatter.png"), p2, width = 7.5, height = 6.5, dpi = 300)

cat("\nPooled pan-cancer PCC (expHRD_exploratory vs HRDsum): r =", pcc_pooled, ", P =", pval_pooled, ", n =", nrow(merged_valid), "\n")
cat("Within-cancer-type PCC (E2, primary estimand): r =", pcc_within, ", P =", pval_within, "\n")
cat("Between-cancer-type PCC (E1): r =", pcc_between, ", P =", pval_between, ", n types =", nrow(type_means), "\n")
cat("Saved plots to:", normalizePath(output_dir), "\n")
cat(" - expHRD_by_cancertype_boxplot.png\n - expHRD_vs_HRDsum_scatter.png\n - merged_expHRD_master.tsv\n")

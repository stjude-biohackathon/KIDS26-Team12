# adapt_ddr_scores.R
# Reshapes beta_hrd_genes.tsv (wide: probes x samples) into long format,
# joins gene names, and joins sample metadata/HRD scores.
# Output: one tidy data.frame ready for VizModules.

library(data.table)

adapt_ddr_scores <- function(
  beta_path   = "/Users/sdowning/KIDS26-Team12/front_end_data/ddr_scars/beta_hrd_genes.tsv",
  probes_path = "/Users/sdowning/KIDS26-Team12/front_end_data/ddr_scars/hrd_gene_probes.tsv",
  samples_path = "/Users/sdowning/KIDS26-Team12/front_end_data/ddr_scars/master_samples.tsv"
) {
  beta <- fread(beta_path)
  probes <- fread(probes_path, header = FALSE, col.names = c("probe_id", "gene"))
  samples <- fread(samples_path)

  # wide -> long: one row per probe x sample
  beta_long <- melt(
    beta,
    id.vars = "probe_id",
    variable.name = "sample_id",
    value.name = "beta_value"
  )

  # add gene name per probe
  beta_long <- merge(beta_long, probes, by = "probe_id", all.x = TRUE)

  # join sample metadata + HRD scores
  out <- merge(beta_long, samples, by = "sample_id", all.x = TRUE)

  out[]
}

# Example use:
# ddr_data <- adapt_ddr_scores()
# head(ddr_data)

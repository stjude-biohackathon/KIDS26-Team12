library(data.table)

adapt_ddr_scores_repo_root <- function(start = getwd()) {
  candidates <- unique(normalizePath(
    c(start, file.path(start, ".."), file.path(start, "../..")),
    winslash = "/",
    mustWork = FALSE
  ))

  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "front_end_data", "ddr_scars", "master_samples.tsv")) &&
        file.exists(file.path(candidate, "R", "adapters", "adapt_ddr_scores.R"))) {
      return(candidate)
    }
  }

  stop("Could not locate the repository root for DDR score inputs.", call. = FALSE)
}

adapt_ddr_scores <- function(beta_path = NULL,
                             probes_path = NULL,
                             samples_path = NULL) {
  if (is.null(beta_path) || is.null(probes_path) || is.null(samples_path)) {
    repo_root <- adapt_ddr_scores_repo_root()
    base_dir <- file.path(repo_root, "front_end_data", "ddr_scars")
  }

  if (is.null(beta_path)) {
    beta_path <- file.path(base_dir, "beta_hrd_genes.tsv")
  }
  if (is.null(probes_path)) {
    probes_path <- file.path(base_dir, "hrd_gene_probes.tsv")
  }
  if (is.null(samples_path)) {
    samples_path <- file.path(base_dir, "master_samples.tsv")
  }

  paths <- c(beta_path = beta_path, probes_path = probes_path, samples_path = samples_path)
  missing_paths <- paths[!file.exists(paths)]
  if (length(missing_paths)) {
    stop(
      "DDR scores input file(s) missing: ",
      paste(sprintf("%s=%s", names(missing_paths), unname(missing_paths)), collapse = ", "),
      call. = FALSE
    )
  }

  beta <- fread(beta_path)
  probes <- fread(probes_path, header = FALSE, col.names = c("probe_id", "gene"))
  samples <- fread(samples_path)

  if (!nrow(beta) || !nrow(samples)) {
    return(data.table())
  }

  beta_long <- melt(
    beta,
    id.vars = "probe_id",
    variable.name = "sample_id",
    value.name = "beta_value"
  )
  beta_long <- merge(beta_long, probes, by = "probe_id", all.x = TRUE)
  merge(beta_long, samples, by = "sample_id", all.x = TRUE)
}

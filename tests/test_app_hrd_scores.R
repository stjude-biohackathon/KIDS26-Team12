#!/usr/bin/env Rscript

if (!requireNamespace("data.table", quietly = TRUE)) stop("Install dependencies with scripts/setup.R")

source("R/app_hrd_scores_helpers.R")
source("R/adapters/adapt_ddr_scores.R")

pass <- function(n, msg) cat(sprintf("  [%s] PASS  %s\n", n, msg))

cat("tests/test_app_hrd_scores.R - HRD Scores app helpers\n\n")

empty <- prepare_hrd_scores_data(data.table::data.table())
stopifnot(
  identical(names(empty), c("sample_id", "cancer_type", "HRDsum", "HRD_LOH", "LST", "TAI", "purity", "ploidy", "epi_HRD")),
  nrow(empty) == 0
)
pass(1, "empty DDR input returns an empty display table with the expected schema")

seed_input <- data.table::data.table(
  sample_id = sprintf("S%03d", 1:600),
  cancer_type = rep(c("A", "B", "C"), length.out = 600),
  HRDsum = seq_len(600),
  HRD_LOH = seq_len(600) %% 7,
  LST = seq_len(600) %% 11,
  TAI = seq_len(600) %% 13,
  purity = 0.8,
  ploidy = 2.1
)

prepared_a <- prepare_hrd_scores_data(seed_input)
prepared_b <- prepare_hrd_scores_data(seed_input)
stopifnot(
  nrow(prepared_a) == 500,
  identical(prepared_a, prepared_b),
  "epi_HRD" %in% names(prepared_a)
)
pass(2, "display subset and simulated epi_HRD values are reproducible and capped at 500 rows")

selected <- filter_hrd_scores_data(
  prepared_a,
  sample_ids = c(as.character(prepared_a$sample_id[3]), "S999", as.character(prepared_a$sample_id[1])),
  sample_mode = "selected"
)
stopifnot(
  identical(as.character(selected$data$sample_id), c(as.character(prepared_a$sample_id[3]), as.character(prepared_a$sample_id[1]))),
  identical(selected$unavailable, "S999")
)
pass(3, "selected-sample mode preserves requested order and reports unavailable samples")

component_data <- hrd_component_data(selected$data)
stopifnot(
  identical(levels(component_data$component), c("HRD_LOH", "LST", "TAI")),
  identical(levels(component_data$sample_id), as.character(selected$data$sample_id))
)
pass(4, "HRD component data preserves stacked-component order and displayed sample order")

missing_error <- tryCatch({
  adapt_ddr_scores(beta_path = tempfile(fileext = ".tsv"), probes_path = tempfile(fileext = ".tsv"), samples_path = tempfile(fileext = ".tsv"))
  NULL
}, error = conditionMessage)
stopifnot(is.character(missing_error), grepl("missing", missing_error, ignore.case = TRUE))
pass(5, "adapter reports missing DDR input files clearly")

scratch <- tempfile(); dir.create(scratch)
on.exit(unlink(scratch, recursive = TRUE), add = TRUE)

beta_path <- file.path(scratch, "beta.tsv")
probes_path <- file.path(scratch, "probes.tsv")
samples_path <- file.path(scratch, "samples.tsv")
writeLines(c("probe_id\tS1\tS2", "cg1\t0.1\t0.2"), beta_path)
writeLines(c("cg1\tBRCA1"), probes_path)
writeLines(c(
  "sample_id\tcancer_type\tHRDsum\tHRD_LOH\tLST\tTAI\tpurity\tploidy",
  "S1\tBRCA\t10\t2\t3\t5\t0.8\t2.1",
  "S2\tOV\t12\t4\t3\t5\t0.7\t2.3"
), samples_path)

adapted <- adapt_ddr_scores(beta_path = beta_path, probes_path = probes_path, samples_path = samples_path)
stopifnot(
  nrow(adapted) == 2,
  all(c("sample_id", "probe_id", "gene", "beta_value", "cancer_type", "HRDsum") %in% names(adapted))
)
pass(6, "adapter reshapes and joins a minimal DDR scores input set")

cat("\nALL 6 HRD SCORES APP TESTS PASSED\n")

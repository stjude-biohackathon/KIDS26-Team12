#!/usr/bin/env Rscript

stopifnot(file.exists("app/app.R"), file.exists("R/adapters/adapt_ddr_scores.R"), file.exists("scripts/setup.R"))

pass <- function(n, msg) cat(sprintf("  [%s] PASS  %s\n", n, msg))

cat("tests/test_app_governance.R - app integration guards\n\n")

app_src <- paste(readLines("app/app.R", warn = FALSE), collapse = "\n")
adapter_src <- paste(readLines("R/adapters/adapt_ddr_scores.R", warn = FALSE), collapse = "\n")
setup_src <- paste(readLines("scripts/setup.R", warn = FALSE), collapse = "\n")
readme_src <- paste(readLines("README.md", warn = FALSE), collapse = "\n")

stopifnot(grepl("shiny::runApp\\('app'\\)", readme_src))
pass(1, "README keeps the repository launch convention `shiny::runApp('app')`")

stopifnot(!grepl("shiny::runApp\\(", app_src), !grepl("/Users/", app_src), !grepl("/Users/", adapter_src))
pass(2, "app and adapter avoid embedded runApp() calls and hard-coded local clone paths")

stopifnot(
  grepl('source\\(file.path\\(repo_root, "R", "app_hrd_scores_helpers\\.R"\\)', app_src),
  grepl('source\\(file.path\\(repo_root, "R", "adapters", "adapt_ddr_scores\\.R"\\)', app_src),
  grepl('file.path\\(base_dir, "beta_hrd_genes\\.tsv"\\)', adapter_src)
)
pass(3, "the app and adapter resolve HRD inputs through repository-safe paths")

stopifnot(grepl('menuItem\\("HRD Scores"', app_src), !grepl('menuItem\\("Methylation"', app_src))
pass(4, "the sidebar replaces the methylation page with the HRD Scores page")

stopifnot(
  grepl('hrd_scores_defaults <- list\\("x.by" = "HRDsum", "y.by" = "epi_HRD", "color.by" = "cancer_type"\\)', app_src),
  grepl('plotlyOutput\\("hrd_component_bar"', app_src),
  grepl('DTOutput\\("hrd_scores_table"', app_src),
  grepl('radioButtons\\([[:space:]]+"hrd_sample_mode"', app_src)
)
pass(5, "HRD Scores controls and outputs are wired with the requested defaults")

stopifnot(
  grepl('hrd_component_colors <- c\\(HRD_LOH = ', app_src),
  grepl('barmode = "stack"', app_src),
  grepl('categoryarray = levels\\(plot_data\\$sample_id\\)', app_src)
)
pass(6, "the component bar chart stays stacked and preserves displayed sample order")

stopifnot(
  grepl('"shinydashboard"', setup_src),
  grepl('"plotly"', setup_src),
  grepl('"DT"', setup_src),
  grepl('"VizModules"', setup_src)
)
pass(7, "scripts/setup.R declares the app runtime dependencies")

cat("\nALL 7 APP INTEGRATION TESTS PASSED\n")

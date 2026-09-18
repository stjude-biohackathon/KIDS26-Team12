#!/usr/bin/env Rscript

if (!requireNamespace("shiny", quietly = TRUE)) stop("Install dependencies with scripts/setup.R")
if (!requireNamespace("shinydashboard", quietly = TRUE)) stop("Install dependencies with scripts/setup.R")
if (!requireNamespace("data.table", quietly = TRUE)) stop("Install dependencies with scripts/setup.R")
if (!requireNamespace("plotly", quietly = TRUE)) stop("Install dependencies with scripts/setup.R")
if (!requireNamespace("DT", quietly = TRUE)) stop("Install dependencies with scripts/setup.R")
if (!requireNamespace("VizModules", quietly = TRUE)) stop("Install dependencies with scripts/setup.R")

stopifnot(file.exists("app/app.R"), file.exists("R/adapters/adapt_ddr_scores.R"), file.exists("scripts/setup.R"))

pass <- function(n, msg) cat(sprintf("  [%s] PASS  %s\n", n, msg))

cat("tests/test_app_governance.R - app integration guards\n\n")

app_src <- paste(readLines("app/app.R", warn = FALSE), collapse = "\n")
adapter_src <- paste(readLines("R/adapters/adapt_ddr_scores.R", warn = FALSE), collapse = "\n")
setup_src <- paste(readLines("scripts/setup.R", warn = FALSE), collapse = "\n")
readme_src <- paste(readLines("README.md", warn = FALSE), collapse = "\n")

e <- new.env()
suppressWarnings(suppressMessages(source("app/app.R", local = e)))

stopifnot(grepl("shiny::runApp\\('app'\\)", readme_src), !grepl("shiny::runApp\\(", app_src))
pass(1, "README keeps the repository launch convention and app/app.R stays deployable")

stopifnot(!grepl("/Users/", app_src), !grepl("/Users/", adapter_src))
pass(2, "app and adapter avoid hard-coded local clone paths")

root_from_repo <- normalizePath(e$app_repo_root(), winslash = "/")
old_wd <- getwd()
setwd("app")
root_from_app <- normalizePath(e$app_repo_root(), winslash = "/")
setwd(old_wd)
stopifnot(identical(root_from_repo, root_from_app))
pass(3, "repo-root discovery works from both the repository root and the app directory")

stopifnot(
  all(c("sample_id", "cancer_type", "HRDsum", "HRD_LOH", "LST", "TAI", "purity", "ploidy") %in% names(e$sample_data)),
  identical(names(e$display_data), c("sample_id", "cancer_type", "HRDsum", "HRD_LOH", "LST", "TAI", "purity", "ploidy", "epi_HRD")),
  nrow(e$display_data) <= 500
)
pass(4, "sourcing the app prepares the expected runtime HRD sample tables")

stopifnot(grepl('menuItem\\("HRD Scores"', app_src), !grepl('menuItem\\("Methylation"', app_src))
pass(5, "the sidebar replaces the methylation page with the HRD Scores page")

stopifnot(
  identical(e$hrd_scores_defaults, list("x.by" = "HRDsum", "y.by" = "epi_HRD", "color.by" = "cancer_type")),
  identical(names(e$hrd_component_colors), c("HRD_LOH", "LST", "TAI")),
  grepl('barmode = "stack"', app_src),
  grepl('categoryarray = levels\\(plot_data\\$sample_id\\)', app_src)
)
pass(6, "the HRD Scores outputs keep the requested defaults and stacked component ordering")

stopifnot(
  grepl('"shinydashboard"', setup_src),
  grepl('"plotly"', setup_src),
  grepl('"DT"', setup_src),
  grepl('"VizModules"', setup_src)
)
pass(7, "scripts/setup.R declares the app runtime dependencies")

cat("\nALL 7 APP INTEGRATION TESTS PASSED\n")

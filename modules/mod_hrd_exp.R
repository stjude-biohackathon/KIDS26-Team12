# Tab 3 — HRD Transcriptome (exp-HRD)
# Genome HRD vs transcriptome-derived exp-HRD, plus purity QC.
#
# exp_HRD is a PLACEHOLDER until Evan/Kayode deliver the real classifier
# output / expression scores. Replace the placeholder block in
# prepare_hrd_exp_sample_data() when that lands.
#
# Usage (from repo root):
#   source("modules/mod_hrd_exp.R")
#   shiny::runApp(hrd_exp_demo_app())
#
# Standalone VizModules check Suzy suggested:
#   sample_data <- prepare_hrd_exp_sample_data()
#   VizModules::dittoViz_scatterPlotApp(
#     data_list = list("hrd_exp" = as.data.frame(sample_data)),
#     defaults = list(x.by = "HRDsum", y.by = "exp_HRD", color.by = "cancer_type")
#   )

library(shiny)
library(data.table)
library(VizModules)

#' Build per-sample table for the exp-HRD tab.
#'
#' Uses genome HRD from front_end_data via adapt_ddr_scores().
#' Adds placeholder exp_HRD = HRDsum + N(0, 3) until real scores arrive.
prepare_hrd_exp_sample_data <- function(
  data_dir = file.path("front_end_data", "ddr_scars"),
  seed = 2L
) {
  adapter <- file.path("R", "adapters", "adapt_ddr_scores.R")
  if (!file.exists(adapter)) {
    stop("Run from the KIDS26-Team12 repo root so R/adapters/adapt_ddr_scores.R is visible.",
         call. = FALSE)
  }
  source(adapter, local = TRUE)

  ddr_data <- adapt_ddr_scores(data_dir = data_dir)
  sample_data <- unique(ddr_data[, .(sample_id, cancer_type, HRDsum, purity, ploidy)])

  # PLACEHOLDER exp-HRD — REMOVE once real classifier output arrives
  set.seed(seed)
  sample_data[, exp_HRD := HRDsum + rnorm(.N, 0, 3)]

  # Agreement / Bland-Altman helpers for true comparison views
  sample_data[, hrd_mean := (HRDsum + exp_HRD) / 2]
  sample_data[, hrd_diff := exp_HRD - HRDsum]
  sample_data[, exp_HRD_source := "PLACEHOLDER_simulated"]

  as.data.frame(sample_data)
}

#' UI for Tab 3 — transcriptome exp-HRD vs genome HRD.
#'
#' @param sample_data data.frame used to populate VizModules column pickers
hrdExpUI <- function(id, sample_data = NULL) {
  ns <- NS(id)

  if (is.null(sample_data)) {
    sample_data <- data.frame(
      sample_id = character(),
      cancer_type = character(),
      HRDsum = numeric(),
      purity = numeric(),
      ploidy = numeric(),
      exp_HRD = numeric(),
      hrd_mean = numeric(),
      hrd_diff = numeric(),
      exp_HRD_source = character(),
      stringsAsFactors = FALSE
    )
  } else {
    sample_data <- as.data.frame(sample_data)
  }

  tagList(
    h3("HRD Transcriptome (exp-HRD)"),
    p(
      "Compares transcriptome-derived exp-HRD with genome HRDsum.",
      strong("exp_HRD is placeholder noise around HRDsum"),
      "until real classifier / expression data arrive.",
      "Purity is shown for QC only — not used to filter locked partitions."
    ),
    checkboxInput(
      ns("show_genome"),
      "Show genome-HRD panel",
      value = TRUE
    ),
    conditionalPanel(
      condition = sprintf("input['%s']", ns("show_genome")),
      wellPanel(
        h4("Genome HRD (hideable)"),
        p("Reference scar burden from master_samples (HRDsum by cancer type)."),
        fluidRow(
          column(
            4,
            dittoViz_yPlotInputsUI(
              ns("genome"),
              data = sample_data,
              defaults = list(
                var = "HRDsum",
                group.by = "cancer_type",
                color.by = "cancer_type",
                plots = "boxplot"
              ),
              title = "Genome HRD controls",
              columns = 1
            )
          ),
          column(8, dittoViz_yPlotOutputUI(ns("genome")))
        )
      )
    ),
    wellPanel(
      h4("exp-HRD vs genome-HRD (comparison)"),
      p(
        "Scatter colored by cancer type, identity line (agreement reference),",
        "and linear regression (lm). Optional faceting via the split control."
      ),
      fluidRow(
        column(
          4,
          dittoViz_scatterPlotInputsUI(
            ns("compare"),
            data = sample_data,
            defaults = list(
              x.by = "HRDsum",
              y.by = "exp_HRD",
              color.by = "cancer_type",
              add.abline = TRUE,
              abline.slope = 1,
              abline.linetype = "dashed",
              custom.model.enable = TRUE,
              custom.models = list(
                models1 = list(
                  model_type = "lm",
                  formula = "exp_HRD ~ HRDsum",
                  line_colour = "#1F77B4",
                  line_width = 2
                )
              )
            ),
            title = "Comparison controls",
            columns = 1
          )
        ),
        column(8, dittoViz_scatterPlotOutputUI(ns("compare")))
      )
    ),
    wellPanel(
      h4("Agreement plot (mean vs difference)"),
      p("Bland-Altman style: x = mean(HRDsum, exp_HRD), y = exp_HRD − HRDsum."),
      fluidRow(
        column(
          4,
          dittoViz_scatterPlotInputsUI(
            ns("agree"),
            data = sample_data,
            defaults = list(
              x.by = "hrd_mean",
              y.by = "hrd_diff",
              color.by = "cancer_type",
              add.yline = 0,
              yline.linetype = "dashed"
            ),
            title = "Agreement controls",
            columns = 1
          )
        ),
        column(8, dittoViz_scatterPlotOutputUI(ns("agree")))
      )
    ),
    wellPanel(
      h4("QC — purity vs exp-HRD"),
      p("Confound check: does placeholder / future exp-HRD track tumor purity?"),
      fluidRow(
        column(
          4,
          dittoViz_scatterPlotInputsUI(
            ns("qc"),
            data = sample_data,
            defaults = list(
              x.by = "purity",
              y.by = "exp_HRD",
              color.by = "cancer_type",
              custom.model.enable = TRUE,
              custom.models = list(
                models1 = list(
                  model_type = "lm",
                  formula = "exp_HRD ~ purity",
                  line_colour = "#E63946",
                  line_width = 2
                )
              )
            ),
            title = "QC controls",
            columns = 1
          )
        ),
        column(8, dittoViz_scatterPlotOutputUI(ns("qc")))
      )
    )
  )
}

#' Server for Tab 3.
#'
#' @param id module id
#' @param sample_data optional data.frame; if NULL, built via prepare_hrd_exp_sample_data()
hrdExpServer <- function(id, sample_data = NULL) {
  moduleServer(id, function(input, output, session) {
    data_reactive <- reactive({
      if (is.null(sample_data)) {
        prepare_hrd_exp_sample_data()
      } else if (is.reactive(sample_data)) {
        sample_data()
      } else {
        as.data.frame(sample_data)
      }
    })

    genome_defaults <- list(
      var = "HRDsum",
      group.by = "cancer_type",
      color.by = "cancer_type",
      plots = "boxplot"
    )
    compare_defaults <- list(
      x.by = "HRDsum",
      y.by = "exp_HRD",
      color.by = "cancer_type",
      add.abline = TRUE,
      abline.slope = 1,
      abline.linetype = "dashed",
      custom.model.enable = TRUE,
      custom.models = list(
        models1 = list(
          model_type = "lm",
          formula = "exp_HRD ~ HRDsum",
          line_colour = "#1F77B4",
          line_width = 2
        )
      )
    )
    agree_defaults <- list(
      x.by = "hrd_mean",
      y.by = "hrd_diff",
      color.by = "cancer_type",
      add.yline = 0,
      yline.linetype = "dashed"
    )
    qc_defaults <- list(
      x.by = "purity",
      y.by = "exp_HRD",
      color.by = "cancer_type",
      custom.model.enable = TRUE,
      custom.models = list(
        models1 = list(
          model_type = "lm",
          formula = "exp_HRD ~ purity",
          line_colour = "#E63946",
          line_width = 2
        )
      )
    )

    dittoViz_yPlotServer(
      "genome",
      data = data_reactive,
      defaults = genome_defaults
    )
    dittoViz_scatterPlotServer(
      "compare",
      data = data_reactive,
      defaults = compare_defaults
    )
    dittoViz_scatterPlotServer(
      "agree",
      data = data_reactive,
      defaults = agree_defaults
    )
    dittoViz_scatterPlotServer(
      "qc",
      data = data_reactive,
      defaults = qc_defaults
    )
  })
}

#' Demo app wrapping the Tab 3 module.
hrd_exp_demo_app <- function(sample_data = NULL) {
  if (is.null(sample_data)) {
    sample_data <- prepare_hrd_exp_sample_data()
  }

  shinyApp(
    ui = fluidPage(
      titlePanel("KIDS26 Tab 3 — HRD Transcriptome (exp-HRD)"),
      hrdExpUI("tab3", sample_data = sample_data)
    ),
    server = function(input, output, session) {
      hrdExpServer("tab3", sample_data = sample_data)
    }
  )
}

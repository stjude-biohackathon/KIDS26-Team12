# modules/mod_hrd_exp.R
# Tab 3 — HRD Transcriptome (exp-HRD)
#
# Data recipe from the task guide. exp_HRD is a PLACEHOLDER until real
# classifier output arrives — replace that one assignment when it lands.
#
# From repo root:
#   source("modules/mod_hrd_exp.R")
#   shiny::runApp(hrd_exp_demo_app())

library(shiny)
library(data.table)
library(VizModules)

prepare_hrd_exp_sample_data <- function() {
  source("R/adapters/adapt_ddr_scores.R", local = TRUE)

  ddr_data <- adapt_ddr_scores()
  sample_data <- unique(ddr_data[, .(sample_id, cancer_type, HRDsum, purity, ploidy)])

  # placeholder exp-HRD score — REMOVE once real classifier output arrives
  set.seed(2)
  sample_data$exp_HRD <- sample_data$HRDsum + rnorm(nrow(sample_data), 0, 3)

  as.data.frame(sample_data)
}

hrdExpUI <- function(id, sample_data) {
  ns <- NS(id)

  tagList(
    checkboxInput(ns("show_genome"), "Show genome-HRD panel", value = TRUE),
    conditionalPanel(
      condition = "input.show_genome",
      ns = ns,
      h4("Genome HRD"),
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
            )
          )
        ),
        column(8, dittoViz_yPlotOutputUI(ns("genome")))
      )
    ),
    h4("exp-HRD vs genome-HRD"),
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
            # dittoViz: add.abline is the intercept (use 0 for y = x, not TRUE)
            add.abline = 0,
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
        )
      ),
      column(8, dittoViz_scatterPlotOutputUI(ns("compare")))
    ),
    h4("QC — purity vs exp-HRD"),
    fluidRow(
      column(
        4,
        dittoViz_scatterPlotInputsUI(
          ns("qc"),
          data = sample_data,
          defaults = list(
            x.by = "purity",
            y.by = "exp_HRD",
            color.by = "cancer_type"
          )
        )
      ),
      column(8, dittoViz_scatterPlotOutputUI(ns("qc")))
    )
  )
}

hrdExpServer <- function(id, sample_data) {
  moduleServer(id, function(input, output, session) {
    data_reactive <- reactive(sample_data)

    dittoViz_yPlotServer(
      "genome",
      data = data_reactive,
      defaults = list(
        var = "HRDsum",
        group.by = "cancer_type",
        color.by = "cancer_type",
        plots = "boxplot"
      )
    )

    dittoViz_scatterPlotServer(
      "compare",
      data = data_reactive,
      defaults = list(
        x.by = "HRDsum",
        y.by = "exp_HRD",
        color.by = "cancer_type",
        add.abline = 0,
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
    )

    dittoViz_scatterPlotServer(
      "qc",
      data = data_reactive,
      defaults = list(
        x.by = "purity",
        y.by = "exp_HRD",
        color.by = "cancer_type"
      )
    )
  })
}

hrd_exp_demo_app <- function() {
  sample_data <- prepare_hrd_exp_sample_data()
  shinyApp(
    ui = fluidPage(
      titlePanel("HRD Transcriptome (exp-HRD)"),
      hrdExpUI("tab3", sample_data)
    ),
    server = function(input, output, session) {
      hrdExpServer("tab3", sample_data)
    }
  )
}

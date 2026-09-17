# modules/mod_hrd_exp.R
# Tab 3 — HRD Transcriptome (exp-HRD)
#
# Data recipe from the task guide. exp_HRD is a PLACEHOLDER until real
# classifier output arrives — replace that one assignment when it lands.
#
# From repo root:
#   source("modules/mod_hrd_exp.R")
#   shiny::runApp(hrd_exp_demo_app())
#
# Standalone VizModules check from the task guide:
#   sample_data <- prepare_hrd_exp_sample_data()
#   VizModules::dittoViz_scatterPlotApp(
#     data_list = list(hrd_exp = sample_data),
#     defaults = list(x.by = "HRDsum", y.by = "exp_HRD", color.by = "cancer_type")
#   )

library(shiny)
library(data.table)
library(dittoViz)
library(plotly)

prepare_hrd_exp_sample_data <- function() {
  source("R/adapters/adapt_ddr_scores.R", local = TRUE)

  ddr_data <- adapt_ddr_scores()
  sample_data <- unique(ddr_data[, .(sample_id, cancer_type, HRDsum, purity, ploidy)])

  # placeholder exp-HRD score — REMOVE once real classifier output arrives
  set.seed(2)
  sample_data$exp_HRD <- sample_data$HRDsum + rnorm(nrow(sample_data), 0, 3)

  as.data.frame(sample_data)
}

# Fast plotly wrapper: dittoViz ggplot -> plotly, no hover lag
.to_plotly <- function(p) {
  ggplotly(p, tooltip = "none") |>
    layout(margin = list(l = 40, r = 10, t = 30, b = 40)) |>
    config(displayModeBar = FALSE)
}

hrdExpUI <- function(id) {
  ns <- NS(id)

  tagList(
    tags$style(HTML(sprintf("
      #%s {
        height: 100vh;
        overflow: hidden;
        padding: 8px 12px;
        box-sizing: border-box;
      }
      #%s .hrd-plots {
        height: calc(100vh - 48px);
      }
      #%s .hrd-plots > .row, #%s .hrd-plots .col {
        height: 100%%;
      }
      #%s h5 { margin: 0 0 4px 0; font-size: 13px; }
    ", ns("wrap"), ns("wrap"), ns("wrap"), ns("wrap"), ns("wrap")))),
    div(
      id = ns("wrap"),
      checkboxInput(ns("show_genome"), "Show genome-HRD panel", value = TRUE),
      uiOutput(ns("plots"))
    )
  )
}

hrdExpServer <- function(id, sample_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$plots <- renderUI({
      show <- isTRUE(input$show_genome)
      w <- if (show) 4 else 6

      cols <- list()
      if (show) {
        cols <- c(cols, list(column(
          w,
          h5("Genome HRD"),
          plotlyOutput(ns("genome"), height = "calc(100vh - 72px)")
        )))
      }
      cols <- c(cols, list(
        column(
          w,
          h5("exp-HRD vs genome-HRD"),
          plotlyOutput(ns("compare"), height = "calc(100vh - 72px)")
        ),
        column(
          w,
          h5("QC — purity vs exp-HRD"),
          plotlyOutput(ns("qc"), height = "calc(100vh - 72px)")
        )
      ))

      div(class = "hrd-plots", do.call(fluidRow, cols))
    })

    output$genome <- renderPlotly({
      req(isTRUE(input$show_genome))
      p <- dittoViz::yPlot(
        sample_data,
        var = "HRDsum",
        group.by = "cancer_type",
        color.by = "cancer_type",
        plots = "boxplot",
        do.hover = FALSE,
        legend.show = FALSE,
        main = NULL,
        xlab = NULL
      )
      .to_plotly(p)
    })

    output$compare <- renderPlotly({
      p <- dittoViz::scatterPlot(
        sample_data,
        x.by = "HRDsum",
        y.by = "exp_HRD",
        color.by = "cancer_type",
        # dittoViz: add.abline is the intercept (0 => y = x)
        add.abline = 0,
        abline.slope = 1,
        abline.linetype = "dashed",
        do.hover = FALSE,
        do.raster = TRUE,
        legend.show = FALSE,
        main = NULL
      )
      fit <- stats::lm(exp_HRD ~ HRDsum, data = sample_data)
      p <- p + ggplot2::geom_abline(
        intercept = stats::coef(fit)[[1]],
        slope = stats::coef(fit)[[2]],
        colour = "#1F77B4",
        linewidth = 1
      )
      .to_plotly(p)
    })

    output$qc <- renderPlotly({
      p <- dittoViz::scatterPlot(
        sample_data,
        x.by = "purity",
        y.by = "exp_HRD",
        do.hover = FALSE,
        do.raster = TRUE,
        legend.show = FALSE,
        main = NULL
      )
      .to_plotly(p)
    })
  })
}

hrd_exp_demo_app <- function() {
  sample_data <- prepare_hrd_exp_sample_data()
  shinyApp(
    ui = fluidPage(
      theme = NULL,
      hrdExpUI("tab3")
    ),
    server = function(input, output, session) {
      hrdExpServer("tab3", sample_data)
    },
    options = list(launch.browser = TRUE)
  )
}

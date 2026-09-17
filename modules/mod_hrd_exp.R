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
# Standalone VizModules check:
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
#' Matches Suzy's recipe: adapt_ddr_scores() -> unique sample rows ->
#' placeholder exp_HRD = HRDsum + N(0, 3).
prepare_hrd_exp_sample_data <- function(
  data_dir = file.path("front_end_data", "ddr_scars"),
  seed = 2L,
  facet_top_n = 6L
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

  # Agreement helpers (Bland-Altman)
  sample_data[, hrd_mean := (HRDsum + exp_HRD) / 2]
  sample_data[, hrd_diff := exp_HRD - HRDsum]
  sample_data[, exp_HRD_source := "PLACEHOLDER_simulated"]

  # Faceted comparison uses the largest tissues only (32 full facets is unreadable)
  top <- names(sort(table(sample_data$cancer_type), decreasing = TRUE))
  top <- top[seq_len(min(facet_top_n, length(top)))]
  sample_data[, cancer_facet := fifelse(cancer_type %in% top, cancer_type, NA_character_)]

  as.data.frame(sample_data)
}

.compare_defaults <- function() {
  list(
    x.by = "HRDsum",
    y.by = "exp_HRD",
    color.by = "cancer_type",
    # dittoViz: add.abline is the INTERCEPT (numeric). TRUE becomes 1 -> wrong line.
    add.abline = 0,
    abline.slope = 1,
    abline.linetype = "dashed",
    do.hover = FALSE,
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
}

.agree_defaults <- function() {
  list(
    x.by = "hrd_mean",
    y.by = "hrd_diff",
    color.by = "cancer_type",
    add.yline = 0,
    yline.linetype = "dashed",
    do.hover = FALSE
  )
}

.qc_defaults <- function() {
  list(
    x.by = "purity",
    y.by = "exp_HRD",
    color.by = "cancer_type",
    do.hover = FALSE,
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
}

.genome_defaults <- function() {
  list(
    var = "HRDsum",
    group.by = "cancer_type",
    color.by = "cancer_type",
    plots = "boxplot",
    do.hover = FALSE
  )
}

.facet_defaults <- function() {
  list(
    x.by = "HRDsum",
    y.by = "exp_HRD",
    color.by = "cancer_type",
    split.by = "cancer_facet",
    split.ncol = 3,
    add.abline = 0,
    abline.slope = 1,
    abline.linetype = "dashed",
    do.hover = FALSE,
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
      cancer_facet = character(),
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
      condition = "input.show_genome",
      ns = ns,
      wellPanel(
        h4("Genome HRD (hideable)"),
        p("Reference scar burden from master_samples (HRDsum by cancer type)."),
        fluidRow(
          column(
            4,
            dittoViz_yPlotInputsUI(
              ns("genome"),
              data = sample_data,
              defaults = .genome_defaults(),
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
        "Scatter colored by cancer type, identity line (y = x),",
        "and linear regression (lm). Use Facet tab / split to facet further."
      ),
      fluidRow(
        column(
          4,
          dittoViz_scatterPlotInputsUI(
            ns("compare"),
            data = sample_data,
            defaults = .compare_defaults(),
            title = "Comparison controls",
            columns = 1
          )
        ),
        column(8, dittoViz_scatterPlotOutputUI(ns("compare")))
      )
    ),
    wellPanel(
      h4("Faceted comparison (largest cancer types)"),
      p(
        "Same axes as above, split by cancer_facet (top tissues by N).",
        "Smaller types are omitted from this panel so facets stay readable."
      ),
      fluidRow(
        column(
          4,
          dittoViz_scatterPlotInputsUI(
            ns("facet"),
            data = sample_data[!is.na(sample_data$cancer_facet), , drop = FALSE],
            defaults = .facet_defaults(),
            title = "Facet controls",
            columns = 1
          )
        ),
        column(8, dittoViz_scatterPlotOutputUI(ns("facet")))
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
            defaults = .agree_defaults(),
            title = "Agreement controls",
            columns = 1
          )
        ),
        column(8, dittoViz_scatterPlotOutputUI(ns("agree")))
      )
    ),
    wellPanel(
      h4("QC — purity vs exp-HRD"),
      p(
        "Confound check: does placeholder / future exp-HRD track tumor purity?",
        "Rows with missing purity are dropped by the plotter."
      ),
      fluidRow(
        column(
          4,
          dittoViz_scatterPlotInputsUI(
            ns("qc"),
            data = sample_data,
            defaults = .qc_defaults(),
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
#' @param sample_data optional data.frame; if NULL, built once via prepare_hrd_exp_sample_data()
hrdExpServer <- function(id, sample_data = NULL) {
  moduleServer(id, function(input, output, session) {
    cached <- if (is.null(sample_data)) {
      prepare_hrd_exp_sample_data()
    } else if (is.reactive(sample_data)) {
      NULL
    } else {
      as.data.frame(sample_data)
    }

    data_reactive <- reactive({
      if (is.reactive(sample_data)) {
        as.data.frame(sample_data())
      } else {
        cached
      }
    })

    facet_reactive <- reactive({
      d <- data_reactive()
      d[!is.na(d$cancer_facet), , drop = FALSE]
    })

    dittoViz_yPlotServer(
      "genome",
      data = data_reactive,
      defaults = .genome_defaults()
    )
    dittoViz_scatterPlotServer(
      "compare",
      data = data_reactive,
      defaults = .compare_defaults()
    )
    dittoViz_scatterPlotServer(
      "facet",
      data = facet_reactive,
      defaults = .facet_defaults()
    )
    dittoViz_scatterPlotServer(
      "agree",
      data = data_reactive,
      defaults = .agree_defaults()
    )
    dittoViz_scatterPlotServer(
      "qc",
      data = data_reactive,
      defaults = .qc_defaults()
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

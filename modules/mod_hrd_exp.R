# modules/mod_hrd_exp.R
# Tab 3: HRD Transcriptome (exp-HRD)
#
# Data recipe from the task guide. exp_HRD is a PLACEHOLDER until real
# classifier output arrives - replace that one assignment when it lands.
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
library(ggplot2)

prepare_hrd_exp_sample_data <- function() {
  # Sample-level fields only (Suzy's unique(...)). Prefer the samples table
  # directly so we do not melt the probe x sample beta matrix on every launch;
  # those columns are what adapt_ddr_scores() joins from master_samples.tsv.
  samples_path <- file.path("front_end_data", "ddr_scars", "master_samples.tsv")
  if (!file.exists(samples_path)) {
    stop(
      "Missing required sample table: ", samples_path,
      "\nRun from the repository root with front_end_data/ddr_scars available.",
      call. = FALSE
    )
  }
  samples <- fread(samples_path)
  sample_data <- unique(samples[, .(sample_id, cancer_type, HRDsum, purity, ploidy)])

  # placeholder exp-HRD score - REMOVE once real classifier output arrives
  set.seed(2)
  sample_data$exp_HRD <- sample_data$HRDsum + rnorm(nrow(sample_data), 0, 3)

  as.data.frame(sample_data)
}

# --- plot builders (built once; renderPlot just draws) -----------------------

.theme_panel <- function(base_size = 13) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_blank(),
      plot.margin = margin(8, 10, 8, 8),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 2)
    )
}

build_genome_plot <- function(sample_data) {
  med <- tapply(sample_data$HRDsum, sample_data$cancer_type, median, na.rm = TRUE)
  ord <- names(sort(med, decreasing = FALSE))
  d <- sample_data
  d$cancer_type <- factor(d$cancer_type, levels = ord)

  dittoViz::yPlot(
    d,
    var = "HRDsum",
    group.by = "cancer_type",
    color.by = "cancer_type",
    plots = "boxplot",
    do.hover = FALSE,
    legend.show = FALSE,
    main = NULL,
    xlab = NULL,
    ylab = "HRDsum"
  ) +
    coord_flip() +
    .theme_panel(14) +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 11)
    )
}

build_compare_plot <- function(sample_data, cancer_type) {
  cancer_type <- as.character(cancer_type)[[1]]
  if (!nzchar(cancer_type) || !cancer_type %in% sample_data$cancer_type) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "Select a cancer type", size = 5, colour = "#666") +
        theme_void()
    )
  }

  d <- sample_data[sample_data$cancer_type == cancer_type, , drop = FALSE]
  d$cancer_type <- factor(d$cancer_type, levels = cancer_type)
  lim <- range(c(sample_data$HRDsum, sample_data$exp_HRD), finite = TRUE)

  # One cancer facet at a time: identity line + lm
  ggplot(d, aes(.data$HRDsum, .data$exp_HRD)) +
    geom_point(size = 1.6, alpha = 0.55, stroke = 0, shape = 16, colour = "#2a6f7c") +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      colour = "grey40",
      linewidth = 0.6
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      colour = "#222222",
      linewidth = 0.9,
      na.rm = TRUE
    ) +
    facet_wrap(~cancer_type) +
    coord_cartesian(xlim = lim, ylim = lim) +
    labs(x = "HRDsum", y = "exp_HRD") +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_blank(),
      plot.background = element_rect(fill = "#ffffff", colour = NA),
      panel.background = element_rect(fill = "#ffffff", colour = NA),
      panel.border = element_blank(),
      axis.line = element_line(colour = "#bbbbbb", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#f0f0f0"),
      legend.position = "none",
      strip.text = element_text(size = 13, face = "bold", colour = "#111111", hjust = 0.5),
      strip.background = element_blank(),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      plot.margin = margin(12, 20, 18, 14)
    )
}

default_cancer_type <- function(sample_data) {
  counts <- sort(table(sample_data$cancer_type), decreasing = TRUE)
  names(counts)[[1]]
}

build_qc_plot <- function(sample_data) {
  d <- sample_data[is.finite(sample_data$purity) & is.finite(sample_data$exp_HRD), , drop = FALSE]

  dittoViz::scatterPlot(
    d,
    x.by = "purity",
    y.by = "exp_HRD",
    do.hover = FALSE,
    do.raster = FALSE,
    legend.show = FALSE,
    main = NULL,
    xlab = "purity",
    ylab = "exp_HRD",
    opacity = 0.55,
    size = 0.7
  ) +
    .theme_panel(14)
}

# --- shiny module ------------------------------------------------------------

hrdExpUI <- function(id) {
  ns <- NS(id)

  tagList(
    tags$style(HTML(sprintf("
      #%s {
        min-height: 100vh;
        overflow: hidden;
        padding: 10px 20px 12px;
        box-sizing: border-box;
        background: #ffffff;
        font-family: 'Segoe UI', 'Helvetica Neue', Helvetica, Arial, sans-serif;
        accent-color: #222;
        display: flex;
        flex-direction: column;
      }
      #%s .hrd-placeholder-warning {
        margin: 0 0 10px;
        padding: 8px 12px;
        border: 1px solid #f2c879;
        background: #fff6df;
        color: #6c4a00;
        font-size: 12px;
        line-height: 1.4;
        text-align: center;
      }
      #%s .hrd-bar {
        display: flex;
        align-items: center;
        justify-content: center;
        gap: 12px;
        min-height: 28px;
        margin-bottom: 8px;
      }
      #%s .hrd-bar .shiny-input-container,
      #%s .hrd-bar .form-group {
        margin: 0 !important;
        padding: 0 !important;
        width: auto !important;
        max-width: none !important;
      }
      #%s .hrd-bar label {
        font-size: 12px;
        margin: 0;
        font-weight: 400;
        color: #555;
      }
      #%s .tabbable {
        display: flex;
        flex-direction: column;
        flex: 1 1 auto;
        min-height: 0;
        background: #ffffff;
        border: none;
      }
      #%s .nav-tabs {
        margin: 0 auto;
        border: none;
        border-bottom: 1px solid #ebebeb;
        padding: 0;
        flex: 0 0 auto;
        background: #ffffff;
        display: flex !important;
        justify-content: center;
        width: 100%%;
      }
      #%s .nav-tabs:before,
      #%s .nav-tabs:after { display: none !important; content: none !important; }
      #%s .nav-tabs > li {
        margin-bottom: -1px;
        float: none !important;
      }
      #%s .nav-tabs > li > a {
        font-size: 13px;
        color: #757575;
        border: none !important;
        border-radius: 0;
        padding: 10px 16px 11px;
        background: transparent !important;
        margin-right: 4px;
      }
      #%s .nav-tabs > li > a:hover {
        color: #222;
        background: transparent !important;
        border: none !important;
      }
      #%s .nav-tabs > li.active > a,
      #%s .nav-tabs > li.active > a:hover,
      #%s .nav-tabs > li.active > a:focus {
        color: #111;
        font-weight: 600;
        background: #ffffff !important;
        border: none !important;
        border-bottom: 2px solid #111 !important;
        padding-bottom: 9px;
      }
      #%s .tab-content {
        flex: 1 1 auto;
        overflow: hidden;
        background: #ffffff;
        border: none;
        padding: 10px 16px 10px;
        box-sizing: border-box;
      }
      #%s .tab-pane.active {
        height: 100%%;
        display: flex;
        flex-direction: column;
        align-items: center;
        min-height: 0;
      }
      #%s .hrd-page {
        height: 100%%;
        width: 100%%;
        max-width: 880px;
        margin: 0 auto;
        display: flex;
        flex-direction: column;
        align-items: center;
        min-height: 0;
        background: #ffffff;
      }
      #%s .hrd-page-with-select {
        display: flex;
        flex-direction: column;
        justify-content: flex-start;
        gap: 8px;
      }
      #%s .hrd-page-plot {
        width: 100%%;
        max-width: 880px;
        flex: 1 1 auto;
        min-height: 0;
        background: #ffffff;
        overflow: hidden;
      }
      #%s .hrd-page-plot .shiny-plot-output {
        width: 100%% !important;
        height: 100%% !important;
        max-height: none !important;
        min-height: 0 !important;
        background: #ffffff !important;
        border: none !important;
        overflow: hidden;
      }
      #%s .hrd-page-plot .shiny-plot-output img {
        display: block;
        margin: 0 auto;
        width: 100%%;
        height: 100%% !important;
        object-fit: contain;
        object-position: center bottom;
        background: #ffffff;
        border: none;
      }
      #%s .hrd-compare-select {
        display: flex;
        justify-content: center;
        align-items: center;
        width: 100%%;
        padding: 0 0 4px;
        background: #ffffff;
        position: relative;
        z-index: 2;
      }
      #%s .hrd-compare-select .shiny-input-container {
        margin: 0;
        width: auto;
        text-align: center;
      }
      #%s .hrd-compare-select label {
        font-size: 12px;
        font-weight: 500;
        color: #555;
        margin-bottom: 4px;
        display: block;
        text-align: center;
        width: 100%%;
      }
      #%s .hrd-compare-select select {
        font-size: 13px;
        height: 34px;
        min-width: 160px;
        border: 1px solid #d0d0d0;
        border-radius: 2px;
        padding: 2px 10px;
        background: #ffffff;
        color: #111;
      }
    ", ns("wrap"), ns("wrap"), ns("wrap"), ns("wrap"), ns("wrap"),
       ns("wrap"), ns("wrap"), ns("wrap"), ns("wrap"), ns("wrap"),
       ns("wrap"), ns("wrap"), ns("wrap"), ns("wrap"), ns("wrap"),
       ns("wrap"), ns("wrap"), ns("wrap"), ns("wrap"), ns("wrap"),
       ns("wrap"), ns("wrap"), ns("wrap"), ns("wrap"), ns("wrap"),
       ns("wrap"), ns("wrap")))),
    div(
      id = ns("wrap"),
      div(
        class = "hrd-bar",
        checkboxInput(ns("show_genome"), "Show genome-HRD", value = TRUE)
      ),
      div(
        class = "hrd-placeholder-warning",
        "Synthetic placeholder only: exp_HRD is generated from HRDsum + noise and is not classifier output."
      ),
      uiOutput(ns("pages"))
    )
  )
}

hrdExpServer <- function(id, sample_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    cancer_choices <- sort(unique(as.character(sample_data$cancer_type)))
    default_cancer <- default_cancer_type(sample_data)

    genome_gg <- build_genome_plot(sample_data)
    qc_gg <- build_qc_plot(sample_data)

    output$pages <- renderUI({
      selected_cancer <- isolate(input$cancer)
      if (is.null(selected_cancer) || !selected_cancer %in% cancer_choices) {
        selected_cancer <- default_cancer
      }
      # Keep current tab when rebuilding for genome toggle
      selected_page <- isolate(input$page)
      if (is.null(selected_page)) selected_page <- "comparison"
      if (!isTRUE(input$show_genome) && identical(selected_page, "genome")) {
        selected_page <- "comparison"
      }

      tabs <- list(
        tabPanel(
          title = "Exp-HRD vs genome-HRD",
          value = "comparison",
          div(
            class = "hrd-page hrd-page-with-select",
            div(
              class = "hrd-page-plot",
              plotOutput(ns("compare"), height = "100%")
            ),
            div(
              class = "hrd-compare-select",
              selectInput(
                ns("cancer"),
                label = "Cancer type",
                choices = cancer_choices,
                selected = selected_cancer,
                width = "200px",
                selectize = FALSE
              )
            )
          )
        ),
        tabPanel(
          title = "QC",
          value = "qc",
          div(
            class = "hrd-page",
            div(
              class = "hrd-page-plot",
              plotOutput(ns("qc"), height = "100%")
            )
          )
        )
      )
      if (isTRUE(input$show_genome)) {
        tabs <- c(tabs, list(tabPanel(
          title = "Genome-HRD",
          value = "genome",
          div(
            class = "hrd-page",
            div(
              class = "hrd-page-plot",
              plotOutput(ns("genome"), height = "100%")
            )
          )
        )))
      }
      do.call(
        tabsetPanel,
        c(tabs, list(id = ns("page"), type = "tabs", selected = selected_page))
      )
    })

    output$genome <- renderPlot({
      req(isTRUE(input$show_genome), identical(input$page, "genome"))
      genome_gg
    }, res = 110, width = 900, height = 700)

    output$compare <- renderPlot({
      req(input$cancer)
      build_compare_plot(sample_data, input$cancer)
    }, res = 110, width = 900, height = 700)

    output$qc <- renderPlot({
      req(identical(input$page, "qc"))
      qc_gg
    }, res = 110, width = 900, height = 700)
  })
}

hrd_exp_demo_app <- function() {
  sample_data <- prepare_hrd_exp_sample_data()
  shinyApp(
    ui = fluidPage(
      tags$style(HTML("
        html, body {
          background: #ffffff !important;
          margin: 0;
          padding: 0;
        }
        .container-fluid {
          padding-left: 0 !important;
          padding-right: 0 !important;
          background: #ffffff !important;
        }
      ")),
      hrdExpUI("tab3")
    ),
    server = function(input, output, session) {
      hrdExpServer("tab3", sample_data)
    }
  )
}

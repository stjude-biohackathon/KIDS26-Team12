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
library(ggplot2)

prepare_hrd_exp_sample_data <- function() {
  source("R/adapters/adapt_ddr_scores.R", local = TRUE)

  # Sample-level fields only (Suzy's unique(...)). Prefer the samples table
  # directly so we do not melt the probe x sample beta matrix on every launch;
  # those columns are what adapt_ddr_scores() joins from master_samples.tsv.
  samples_path <- file.path("front_end_data", "ddr_scars", "master_samples.tsv")
  if (!file.exists(samples_path)) {
    ddr_data <- adapt_ddr_scores()
    sample_data <- unique(ddr_data[, .(sample_id, cancer_type, HRDsum, purity, ploidy)])
  } else {
    samples <- fread(samples_path)
    sample_data <- unique(samples[, .(sample_id, cancer_type, HRDsum, purity, ploidy)])
  }

  # placeholder exp-HRD score — REMOVE once real classifier output arrives
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

build_compare_plot <- function(sample_data) {
  lim <- range(c(sample_data$HRDsum, sample_data$exp_HRD), finite = TRUE)
  fit <- stats::lm(exp_HRD ~ HRDsum, data = sample_data)

  dittoViz::scatterPlot(
    sample_data,
    x.by = "HRDsum",
    y.by = "exp_HRD",
    color.by = "cancer_type",
    # dittoViz: add.abline is the intercept (0 => y = x)
    add.abline = 0,
    abline.slope = 1,
    abline.linetype = "dashed",
    abline.color = "grey40",
    do.hover = FALSE,
    do.raster = TRUE,
    legend.show = TRUE,
    main = NULL,
    xlab = "HRDsum (genome)",
    ylab = "exp_HRD"
  ) +
    geom_abline(
      intercept = stats::coef(fit)[[1]],
      slope = stats::coef(fit)[[2]],
      colour = "#1F77B4",
      linewidth = 1.2
    ) +
    coord_cartesian(xlim = lim, ylim = lim) +
    .theme_panel(14) +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 9),
      legend.key.size = grid::unit(0.4, "cm"),
      legend.spacing.y = grid::unit(0.08, "cm")
    ) +
    guides(colour = guide_legend(ncol = 2, override.aes = list(size = 2.2, alpha = 1)))
}

build_qc_plot <- function(sample_data) {
  d <- sample_data[is.finite(sample_data$purity) & is.finite(sample_data$exp_HRD), , drop = FALSE]

  dittoViz::scatterPlot(
    d,
    x.by = "purity",
    y.by = "exp_HRD",
    do.hover = FALSE,
    do.raster = TRUE,
    legend.show = FALSE,
    main = NULL,
    xlab = "purity",
    ylab = "exp_HRD",
    opacity = 0.55
  ) +
    .theme_panel(14)
}

# --- shiny module ------------------------------------------------------------

hrdExpUI <- function(id) {
  ns <- NS(id)

  tagList(
    tags$style(HTML(sprintf("
      #%s {
        height: 100vh;
        overflow: hidden;
        padding: 8px 12px;
        box-sizing: border-box;
        background: #f7f7f7;
      }
      #%s .hrd-bar {
        display: flex;
        align-items: center;
        gap: 18px;
        height: 32px;
        margin-bottom: 6px;
      }
      #%s .hrd-bar .shiny-input-container { margin: 0; padding: 0; }
      #%s .hrd-bar label { font-size: 13px; margin: 0; font-weight: 500; }
      #%s .hrd-note { font-size: 12px; color: #666; margin: 0; }
      #%s .tabbable, #%s .tab-content, #%s .tab-pane {
        height: calc(100vh - 52px);
      }
      #%s .tab-content {
        overflow: hidden;
        background: #fff;
        border: 1px solid #ddd;
        border-top: none;
        border-radius: 0 0 4px 4px;
        padding: 8px 10px;
        box-sizing: border-box;
      }
      #%s .nav-tabs { margin-bottom: 0; }
      #%s .shiny-plot-output {
        height: calc(100vh - 110px) !important;
      }
    ", ns("wrap"), ns("wrap"), ns("wrap"), ns("wrap"), ns("wrap"),
       ns("wrap"), ns("wrap"), ns("wrap"), ns("wrap"), ns("wrap"), ns("wrap")))),
    div(
      id = ns("wrap"),
      div(
        class = "hrd-bar",
        checkboxInput(ns("show_genome"), "Show genome-HRD page", value = TRUE),
        span(class = "hrd-note", "exp_HRD is placeholder until classifier output arrives")
      ),
      uiOutput(ns("pages"))
    )
  )
}

hrdExpServer <- function(id, sample_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    genome_gg <- build_genome_plot(sample_data)
    compare_gg <- build_compare_plot(sample_data)
    qc_gg <- build_qc_plot(sample_data)

    output$pages <- renderUI({
      tabs <- list(
        tabPanel(
          "Comparison",
          plotOutput(ns("compare"), height = "100%")
        ),
        tabPanel(
          "QC",
          plotOutput(ns("qc"), height = "100%")
        )
      )
      if (isTRUE(input$show_genome)) {
        tabs <- c(
          list(tabPanel(
            "Genome HRD",
            plotOutput(ns("genome"), height = "100%")
          )),
          tabs
        )
      }
      do.call(tabsetPanel, c(tabs, list(id = ns("page"), type = "tabs")))
    })

    output$genome <- renderPlot({
      req(isTRUE(input$show_genome))
      genome_gg
    }, res = 120)

    output$compare <- renderPlot({
      compare_gg
    }, res = 120)

    output$qc <- renderPlot({
      qc_gg
    }, res = 120)
  })
}

hrd_exp_demo_app <- function() {
  sample_data <- prepare_hrd_exp_sample_data()
  shinyApp(
    ui = fluidPage(hrdExpUI("tab3")),
    server = function(input, output, session) {
      hrdExpServer("tab3", sample_data)
    }
  )
}

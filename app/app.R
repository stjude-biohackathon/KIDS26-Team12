library(shiny)
library(shinydashboard)
library(VizModules)
library(data.table)
library(plotly)
library(DT)

app_repo_root <- function(start = getwd()) {
  candidates <- unique(normalizePath(
    c(start, file.path(start, ".."), file.path(start, "../..")),
    winslash = "/",
    mustWork = FALSE
  ))

  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "app", "app.R")) &&
        file.exists(file.path(candidate, "R", "adapters", "adapt_ddr_scores.R"))) {
      return(candidate)
    }
  }

  stop("Could not locate the repository root for the Shiny app.", call. = FALSE)
}

repo_root <- app_repo_root()
source(file.path(repo_root, "R", "app_hrd_scores_helpers.R"), local = TRUE)
source(file.path(repo_root, "R", "adapters", "adapt_ddr_scores.R"), local = TRUE)

required_cols <- c("sample_id", "cancer_type", "HRDsum", "HRD_LOH", "LST", "TAI", "purity", "ploidy")
empty_base_sample_data <- empty_hrd_scores_data()[, .(sample_id, cancer_type, HRDsum, HRD_LOH, LST, TAI, purity, ploidy)]

load_ddr_scores <- function() {
  tryCatch(adapt_ddr_scores(), error = function(err) {
    empty <- data.table()
    attr(empty, "load_error") <- conditionMessage(err)
    empty
  })
}

base_sample_data <- function(ddr_data) {
  if (is.null(ddr_data) || !nrow(ddr_data)) {
    return(copy(empty_base_sample_data))
  }

  missing_cols <- setdiff(required_cols, names(ddr_data))
  if (length(missing_cols)) {
    stop(
      "DDR input is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  sample_data <- unique(as.data.table(ddr_data)[, ..required_cols])
  setorder(sample_data, sample_id)
  sample_data
}

no_data_ui <- function(title, detail) {
  div(
    style = "padding: 30px; background: #fff; border-radius: 4px; margin: 15px;",
    h4(title, style = "color: #1b2a4a;"),
    p(detail, style = "margin-bottom: 0;")
  )
}

ddr_data <- load_ddr_scores()
ddr_load_error <- attr(ddr_data, "load_error")
sample_data <- tryCatch(base_sample_data(ddr_data), error = function(err) {
  ddr_load_error <<- conditionMessage(err)
  copy(empty_base_sample_data)
})
display_data <- tryCatch(prepare_hrd_scores_data(ddr_data), error = function(err) {
  ddr_load_error <<- conditionMessage(err)
  empty_hrd_scores_data()
})
cancer_counts <- if (nrow(sample_data)) sample_data[, .N, by = cancer_type] else data.table(cancer_type = character(), N = integer())

has_sample_data <- nrow(sample_data) > 0
has_display_data <- nrow(display_data) > 0
hrd_scores_defaults <- list("x.by" = "HRDsum", "y.by" = "epi_HRD", "color.by" = "cancer_type")
hrd_component_colors <- c(HRD_LOH = "#1B9E77", LST = "#D95F02", TAI = "#7570B3")

logo_svg <- HTML('
<svg width="40" height="40" viewBox="0 0 100 100" xmlns="http://www.w3.org/2000/svg">
  <path d="M20,10 Q35,25 20,40" stroke="#c41230" stroke-width="4" fill="none" stroke-linecap="round"/>
  <path d="M80,10 Q65,25 80,40" stroke="#c41230" stroke-width="4" fill="none" stroke-linecap="round"/>
  <line x1="24" y1="14" x2="76" y2="14" stroke="#8c8c8c" stroke-width="3"/>
  <line x1="27" y1="24" x2="73" y2="24" stroke="#8c8c8c" stroke-width="3"/>
  <line x1="22" y1="34" x2="78" y2="34" stroke="#8c8c8c" stroke-width="3"/>
  <path d="M20,60 Q35,75 20,90" stroke="#c41230" stroke-width="4" fill="none" stroke-linecap="round"/>
  <path d="M80,60 Q65,75 80,90" stroke="#c41230" stroke-width="4" fill="none" stroke-linecap="round"/>
  <line x1="24" y1="64" x2="76" y2="64" stroke="#8c8c8c" stroke-width="3"/>
  <line x1="27" y1="74" x2="73" y2="74" stroke="#8c8c8c" stroke-width="3"/>
  <line x1="22" y1="84" x2="78" y2="84" stroke="#8c8c8c" stroke-width="3"/>
  <path d="M15,48 L30,42 L20,52 L35,50" stroke="#c41230" stroke-width="3" fill="none" stroke-linecap="round" stroke-linejoin="round"/>
  <path d="M85,48 L70,42 L80,52 L65,50" stroke="#c41230" stroke-width="3" fill="none" stroke-linecap="round" stroke-linejoin="round"/>
</svg>')

ui <- dashboardPage(
  skin = "red",
  dashboardHeader(title = span(logo_svg, " Genomic SCARS")),
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem("Welcome", tabName = "welcome", icon = icon("home")),
      menuItem("HRD (Genome)", tabName = "hrd", icon = icon("dna")),
      menuItem("HRD Scores", tabName = "hrd_scores", icon = icon("table")),
      menuItem("Transcriptomics", tabName = "trans", icon = icon("chart-bar"))
    )
  ),
  dashboardBody(
    tags$head(tags$style(HTML('
      .skin-red .main-header .navbar { background-color: #1b2a4a; }
      .skin-red .main-header .logo { background-color: #1b2a4a; color: #fff; }
      .skin-red .main-header .logo:hover { background-color: #243560; }
      .skin-red .main-sidebar { background-color: #1b2a4a; }
      .skin-red .sidebar-menu > li.active > a { background-color: #c41230; color: #fff; }
      .skin-red .sidebar-menu > li > a:hover { background-color: #2a3f6a; color: #fff; }
      .content-wrapper { background-color: #f5f5f5; }
      .nav-tabs > li > a { color: #1b2a4a; }
      .nav-tabs > li.active > a,
      .nav-tabs > li.active > a:hover,
      .nav-tabs > li.active > a:focus {
        color: #c41230;
        border-bottom: 2px solid #c41230;
      }
    '))),
    tabItems(
      tabItem(
        tabName = "welcome",
        div(style = "padding: 40px; max-width: 800px; margin: auto; text-align: center;",
          div(style = "display: flex; justify-content: center;", logo_svg),
          h1("Genomic SCARS", style = "font-weight: bold; letter-spacing: 2px; color: #1b2a4a;"),
          h4("Assessing HRD Across Genome, Scores & Transcriptome",
             style = "color: #666; font-weight: normal;"),
          hr(),
          p(style = "text-align: left; font-size: 16px; line-height: 1.6;",
            "This app explores HRD (Homologous Recombination Deficiency) signals in pediatric and pan-cancer WGS data. HRDsum combines three genomic scar signatures - LOH, LST, and TAI - that together indicate a tumor's loss of homologous recombination repair capacity."
          ),
          p(style = "text-align: left; font-size: 16px; line-height: 1.6;",
            "Use the sidebar to move between genome, HRD score, and transcriptome views. The HRD Scores page compares genomic HRDsum with reproducibly simulated epi_HRD values for a display subset of samples."
          )
        )
      ),
      tabItem(
        tabName = "hrd",
        if (has_sample_data) {
          tabsetPanel(
            tabPanel("HRDsum by Cancer Type",
              plotthis_BoxPlotInputsUI("hrdsum_box", data = sample_data),
              plotthis_BoxPlotOutputUI("hrdsum_box")
            ),
            tabPanel("HRD Sub-Scores",
              plotthis_BoxPlotInputsUI("subscore_box", data = sample_data),
              plotthis_BoxPlotOutputUI("subscore_box")
            ),
            tabPanel("Purity vs HRDsum",
              dittoViz_scatterPlotInputsUI("purity_scatter", data = sample_data),
              dittoViz_scatterPlotOutputUI("purity_scatter")
            ),
            tabPanel("Ploidy vs HRDsum",
              dittoViz_scatterPlotInputsUI("ploidy_scatter", data = sample_data),
              dittoViz_scatterPlotOutputUI("ploidy_scatter")
            ),
            tabPanel("Sample Counts",
              plotthis_BarPlotInputsUI("cancer_bar", data = cancer_counts),
              plotthis_BarPlotOutputUI("cancer_bar")
            )
          )
        } else {
          no_data_ui(
            "Genome HRD data unavailable",
            if (!is.null(ddr_load_error) && nzchar(ddr_load_error)) ddr_load_error else "No DDR score data were loaded."
          )
        }
      ),
      tabItem(
        tabName = "hrd_scores",
        if (has_display_data) {
          fluidRow(
            box(
              width = 12,
              title = "HRD Scores controls",
              status = "primary",
              solidHeader = TRUE,
              radioButtons(
                "hrd_sample_mode",
                "Samples to display",
                choices = c(
                  "All display-subset samples" = "all",
                  "User-selected samples" = "selected"
                ),
                selected = "all",
                inline = TRUE
              ),
              selectizeInput(
                "hrd_selected_samples",
                "User-selected samples",
                choices = sample_data$sample_id,
                selected = NULL,
                multiple = TRUE,
                options = list(placeholder = "Select samples to focus the HRD Scores outputs")
              ),
              helpText(sprintf(
                "A reproducible random display subset of %d sample%s is shown here (maximum 500 rows).",
                nrow(display_data),
                if (nrow(display_data) == 1) "" else "s"
              )),
              uiOutput("hrd_scores_status")
            ),
            box(
              width = 12,
              tabsetPanel(
                tabPanel(
                  "Scatter Plot",
                  dittoViz_scatterPlotInputsUI(
                    "hrd_scores_scatter",
                    data = display_data,
                    defaults = hrd_scores_defaults
                  ),
                  dittoViz_scatterPlotOutputUI("hrd_scores_scatter")
                ),
                tabPanel(
                  "HRD Components",
                  plotlyOutput("hrd_component_bar", height = "600px")
                ),
                tabPanel(
                  "Sample Table",
                  DTOutput("hrd_scores_table")
                )
              )
            )
          )
        } else {
          no_data_ui(
            "HRD Scores data unavailable",
            if (!is.null(ddr_load_error) && nzchar(ddr_load_error)) ddr_load_error else "No DDR score rows were available after loading."
          )
        }
      ),
      tabItem(
        tabName = "trans",
        div(style = "padding: 40px; color: #888;",
          h3("Transcriptomics x HRD", style = "color: #1b2a4a;"),
          p("Pending exp-HRD data. Tabs will be added here once available.")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  if (has_sample_data) {
    plotthis_BoxPlotServer("hrdsum_box", data = reactive(sample_data))
    plotthis_BoxPlotServer("subscore_box", data = reactive(sample_data))
    dittoViz_scatterPlotServer("purity_scatter", data = reactive(sample_data))
    dittoViz_scatterPlotServer("ploidy_scatter", data = reactive(sample_data))
    plotthis_BarPlotServer("cancer_bar", data = reactive(cancer_counts))
  }

  if (has_display_data) {
    active_hrd_scores <- reactive({
      filter_hrd_scores_data(
        display_data,
        sample_ids = input$hrd_selected_samples,
        sample_mode = if (is.null(input$hrd_sample_mode)) "all" else input$hrd_sample_mode
      )
    })

    output$hrd_scores_status <- renderUI({
      active <- active_hrd_scores()

      if (identical(input$hrd_sample_mode, "selected")) {
        if (length(active$unavailable)) {
          return(tags$p(
            style = "color: #c41230; margin-top: 10px;",
            paste(
              "These selected samples are not in the display subset and were skipped:",
              paste(active$unavailable, collapse = ", ")
            )
          ))
        }

        if (!nrow(active$data)) {
          return(tags$p(
            style = "color: #666; margin-top: 10px;",
            "Select one or more samples from the display subset to populate the plots and table."
          ))
        }
      }

      tags$p(
        style = "color: #666; margin-top: 10px;",
        sprintf("Displaying %d sample%s.", nrow(active$data), if (nrow(active$data) == 1) "" else "s")
      )
    })

    dittoViz_scatterPlotServer(
      "hrd_scores_scatter",
      data = reactive(active_hrd_scores()$data),
      defaults = hrd_scores_defaults
    )

    output$hrd_component_bar <- renderPlotly({
      active <- active_hrd_scores()$data
      validate(need(nrow(active) > 0, "No samples available for the HRD component chart."))

      plot_data <- hrd_component_data(active)
      chart <- plot_ly(
        data = plot_data,
        x = ~sample_id,
        y = ~score,
        color = ~component,
        colors = hrd_component_colors,
        type = "bar",
        hovertemplate = paste(
          "Sample: %{x}",
          "<br>Component: %{fullData.name}",
          "<br>Score: %{y}",
          "<extra></extra>"
        )
      )
      layout(
        chart,
        barmode = "stack",
        xaxis = list(
          title = "Sample",
          categoryorder = "array",
          categoryarray = levels(plot_data$sample_id),
          tickangle = -45
        ),
        yaxis = list(title = "HRD component score"),
        legend = list(title = list(text = "Component"))
      )
    })

    output$hrd_scores_table <- renderDT({
      active <- active_hrd_scores()$data
      datatable(
        active,
        filter = "top",
        rownames = FALSE,
        options = list(pageLength = 25, scrollX = TRUE)
      )
    })
  }
}

shinyApp(ui, server)

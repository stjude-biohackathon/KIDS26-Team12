library(shiny)
library(shinydashboard)
library(VizModules)
library(data.table)

source("/Users/sdowning/KIDS26-Team12/R/adapters/adapt_ddr_scores.R")

ddr_data <- adapt_ddr_scores()

sample_data   <- unique(ddr_data[, .(sample_id, cancer_type, HRD_LOH, LST, TAI, HRDsum, purity, ploidy)])
meth_data     <- unique(ddr_data[, .(sample_id, probe_id, gene, beta_value, cancer_type, HRDsum)])
cancer_counts <- sample_data[, .N, by = cancer_type]

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
      menuItem("Methylation", tabName = "meth", icon = icon("flask")),
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

      tabItem(tabName = "welcome",
        div(style = "padding: 40px; max-width: 800px; margin: auto; text-align: center;",
          div(style = "display: flex; justify-content: center;", logo_svg),
          h1("Genomic SCARS", style = "font-weight: bold; letter-spacing: 2px; color: #1b2a4a;"),
          h4("Assessing HRD Across Genome, Methylation & Transcriptome",
             style = "color: #666; font-weight: normal;"),
          hr(),
          p(style = "text-align: left; font-size: 16px; line-height: 1.6;",
            "This app explores HRD (Homologous Recombination Deficiency) signals in
             pediatric and pan-cancer WGS data. HRDsum combines three genomic scar
             signatures - LOH, LST, and TAI - that together indicate a tumor's loss
             of homologous recombination repair capacity."
          ),
          p(style = "text-align: left; font-size: 16px; line-height: 1.6;",
            "Use the sidebar to move between genome, methylation, and transcriptome
             HRD views. Each page has its own set of visualization tabs."
          )
        )
      ),

      tabItem(tabName = "hrd",
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
      ),

      tabItem(tabName = "meth",
        tabsetPanel(
          tabPanel("Methylation by Gene",
            plotthis_BoxPlotInputsUI("meth_box", data = meth_data),
            plotthis_BoxPlotOutputUI("meth_box")
          ),
          tabPanel("Methylation vs HRDsum",
            dittoViz_scatterPlotInputsUI("meth_hrd_scatter", data = meth_data),
            dittoViz_scatterPlotOutputUI("meth_hrd_scatter")
          )
        )
      ),

      tabItem(tabName = "trans",
        div(style = "padding: 40px; color: #888;",
          h3("Transcriptomics x HRD", style = "color: #1b2a4a;"),
          p("Pending exp-HRD data. Tabs will be added here once available.")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  plotthis_BoxPlotServer("hrdsum_box", data = reactive(sample_data))
  plotthis_BoxPlotServer("subscore_box", data = reactive(sample_data))
  dittoViz_scatterPlotServer("purity_scatter", data = reactive(sample_data))
  dittoViz_scatterPlotServer("ploidy_scatter", data = reactive(sample_data))
  plotthis_BoxPlotServer("meth_box", data = reactive(meth_data))
  dittoViz_scatterPlotServer("meth_hrd_scatter", data = reactive(meth_data))
  plotthis_BarPlotServer("cancer_bar", data = reactive(cancer_counts))
}

shinyApp(ui, server)

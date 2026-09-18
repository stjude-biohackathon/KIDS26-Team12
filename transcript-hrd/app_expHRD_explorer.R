#!/usr/bin/env Rscript
#
# Shiny app: interactive explorer for merged expHRD + master_samples data.
# Lets you pick one or more patients (optionally filtered by cancer type) and see:
#   - a searchable/sortable table of their scores
#   - their position highlighted on the expHRD-vs-HRDsum scatter plot
#   - their position highlighted on the cancer-type box plot
#
# Data: expects merged_expHRD_master.tsv produced by merge_and_plot_expHRD.R, with
# columns: patient_id, sample_id, cancer_type, positive_ssGSEA, negative_ssGSEA,
#          expHRD_exploratory, HRDsum, HRD_LOH, LST, TAI, interpretation
#
# Run:
#   Rscript -e 'shiny::runApp("app_expHRD_explorer.R", host="0.0.0.0", port=8787)'
# or, if run interactively / via Rscript directly:
#   Rscript app_expHRD_explorer.R
#
# Data path: set the DATA_FILE environment variable, or edit DATA_FILE below.
# On a remote cluster, forward the port over SSH to view it locally, e.g.:
#   ssh -L 8787:localhost:8787 jvu@<cluster_host>
# then open http://localhost:8787 in your browser.

suppressPackageStartupMessages({
  library(shiny)
  library(data.table)
  library(ggplot2)
  library(DT)
})
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
if (!has_ggrepel) {
  message("Package 'ggrepel' not found - point labels on the scatter plot will be disabled. ",
          "Install with: install.packages('ggrepel') for labeled points.")
}

DATA_FILE <- Sys.getenv("DATA_FILE", unset = "merged_expHRD_master.tsv")
if (!file.exists(DATA_FILE)) {
  stop("Data file not found: ", DATA_FILE,
       "\nSet DATA_FILE env var or place merged_expHRD_master.tsv next to this app.")
}

full_data <- fread(DATA_FILE)
req_cols <- c("patient_id", "cancer_type", "expHRD_exploratory", "HRDsum")
missing <- setdiff(req_cols, names(full_data))
if (length(missing)) stop("Missing expected columns: ", paste(missing, collapse = ", "))

full_data <- full_data[!is.na(cancer_type)]
cancer_types <- sort(unique(full_data$cancer_type))

ui <- fluidPage(
  titlePanel("expHRD Explorer Pan-Cancer TCGA (EBPlusPlus, exploratory)"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("cancer_filter", "Cancer type (patient choices & highlight)",
                  choices = c("All", cancer_types), selected = "All"),
      selectizeInput("patient_select", "Select patient(s)",
                     choices = NULL, multiple = TRUE,
                     options = list(placeholder = "Type a TCGA patient ID...")),
      hr(),
      sliderInput("hrdsum_range", "Filter by HRDsum range",
                  min = floor(min(full_data$HRDsum, na.rm = TRUE)),
                  max = ceiling(max(full_data$HRDsum, na.rm = TRUE)),
                  value = c(floor(min(full_data$HRDsum, na.rm = TRUE)),
                            ceiling(max(full_data$HRDsum, na.rm = TRUE)))),
      helpText("Selected patients are highlighted in red on both plots."),
      helpText(paste0("Loaded ", nrow(full_data), " samples across ",
                       length(cancer_types), " cancer types."))
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Selected patient table", DTOutput("patient_table")),
        tabPanel("expHRD vs HRDsum", plotOutput("scatter_plot", height = "550px")),
        tabPanel("expHRD by cancer type", plotOutput("box_plot", height = "600px"))
      )
    )
  )
)

server <- function(input, output, session) {

  filtered_by_type <- reactive({
    if (input$cancer_filter == "All") full_data else full_data[cancer_type == input$cancer_filter]
  })

  observeEvent(filtered_by_type(), {
    updateSelectizeInput(session, "patient_select",
                          choices = sort(unique(filtered_by_type()$patient_id)),
                          server = TRUE)
  })

  range_filtered <- reactive({
    full_data[HRDsum >= input$hrdsum_range[1] & HRDsum <= input$hrdsum_range[2]]
  })

  selected_data <- reactive({
    d <- range_filtered()
    if (length(input$patient_select)) d[patient_id %in% input$patient_select] else d[0]
  })

  output$patient_table <- renderDT({
    d <- selected_data()
    cols <- intersect(c("patient_id", "sample_id", "cancer_type", "positive_ssGSEA",
                         "negative_ssGSEA", "expHRD_exploratory", "HRDsum",
                         "HRD_LOH", "LST", "TAI", "interpretation"), names(d))
    datatable(d[, ..cols], options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE)
  })

  output$scatter_plot <- renderPlot({
    base <- range_filtered()
    sel <- selected_data()
    ggplot(base, aes(x = HRDsum, y = expHRD_exploratory)) +
      geom_point(alpha = 0.25, color = "#4C72B0", size = 1.2) +
      { if (nrow(sel)) geom_point(data = sel, aes(x = HRDsum, y = expHRD_exploratory),
                                   color = "red", size = 3.2) } +
      { if (nrow(sel) && has_ggrepel) ggrepel::geom_text_repel(data = sel,
                                                  aes(x = HRDsum, y = expHRD_exploratory, label = patient_id),
                                                  size = 3.5, color = "black", max.overlaps = 20) } +
      { if (nrow(sel) && !has_ggrepel) geom_text(data = sel,
                                                  aes(x = HRDsum, y = expHRD_exploratory, label = patient_id),
                                                  size = 3.5, color = "black", vjust = -1) } +
      labs(title = "expHRD (exploratory) vs Reference HRDsum",
           subtitle = paste0(nrow(sel), " patient(s) highlighted"),
           x = "HRDsum (reference)", y = "expHRD_exploratory") +
      theme_minimal(base_size = 13)
  })

  output$box_plot <- renderPlot({
    base <- range_filtered()
    sel <- selected_data()
    type_order <- base[, .(med = median(expHRD_exploratory, na.rm = TRUE)), by = cancer_type][order(med)]$cancer_type
    base$cancer_type <- factor(base$cancer_type, levels = type_order)

    # Highlight the box for the cancer type currently chosen in the "Filter by cancer
    # type" dropdown: distinct fill + thicker black outline vs. the default gray boxes.
    base$is_selected_type <- if (input$cancer_filter != "All") {
      base$cancer_type == input$cancer_filter
    } else {
      FALSE
    }

    p <- ggplot(base, aes(x = cancer_type, y = expHRD_exploratory,
                           fill = is_selected_type, color = is_selected_type,
                           linewidth = is_selected_type)) +
      geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
      scale_fill_manual(values = c("FALSE" = "#B0C4DE", "TRUE" = "#F2C464"), guide = "none") +
      scale_color_manual(values = c("FALSE" = "grey40", "TRUE" = "black"), guide = "none") +
      scale_linewidth_manual(values = c("FALSE" = 0.4, "TRUE" = 1.3), guide = "none") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
      labs(title = "expHRD by Cancer Type", x = "TCGA cancer type", y = "expHRD_exploratory")

    # Bold the selected cancer type's axis label too.
    if (any(base$is_selected_type)) {
      face_vec <- ifelse(levels(base$cancer_type) == input$cancer_filter, "bold", "plain")
      size_vec <- ifelse(levels(base$cancer_type) == input$cancer_filter, 11, 8)
      p <- p + theme(axis.text.x = element_text(angle = 60, hjust = 1,
                                                  face = face_vec, size = size_vec))
    }

    # Selected individual patients: shown as very large, high-contrast points so they
    # stand out clearly against the boxes, with a black outline and a text label.
    if (nrow(sel)) {
      sel$cancer_type <- factor(sel$cancer_type, levels = type_order)
      p <- p + geom_point(data = sel, aes(x = cancer_type, y = expHRD_exploratory),
                           inherit.aes = FALSE, color = "red", fill = "yellow",
                           shape = 21, size = 9, stroke = 1.6) +
        geom_text(data = sel, aes(x = cancer_type, y = expHRD_exploratory, label = patient_id),
                  inherit.aes = FALSE, size = 3.5, fontface = "bold",
                  vjust = -1.4, color = "black")
    }
    p
  })
}

shinyApp(ui, server)

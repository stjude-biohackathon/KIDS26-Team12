library(shiny)
library(data.table)

# ---------------------------------------------------------------------
# Load real committed front-end data
# ---------------------------------------------------------------------

data_dir <- file.path("front_end_data", "ddr_scars")

beta <- fread(
    file.path(data_dir, "beta_hrd_genes.tsv")
)

probes <- fread(
    file.path(data_dir, "hrd_gene_probes.tsv"),
    header = FALSE,
    col.names = c("probe_id", "gene_name")
)

sample_metadata <- fread(
    file.path(data_dir, "master_samples.tsv")
)

# Sample columns are every beta-matrix column except probe_id.
sample_columns <- setdiff(names(beta), "probe_id")

# Extract individual gene symbols from entries such as "ATM;NPAT".
gene_choices <- sort(unique(unlist(
    strsplit(probes$gene_name, ";", fixed = TRUE)
)))

gene_choices <- trimws(gene_choices)
gene_choices <- gene_choices[nzchar(gene_choices)]

# ---------------------------------------------------------------------
# User interface
# ---------------------------------------------------------------------

ui <- fluidPage(
    titlePanel("KIDS26 DDR/HRD methylation viewer"),
    
    p(
        "Real committed beta values from front_end_data/ddr_scars.",
        "The boxplot shows the distribution of mean beta values across samples",
        "for probes assigned to the selected gene."
    ),
    
    selectInput(
        inputId = "gene",
        label = "Select gene",
        choices = gene_choices,
        selected = "BRCA1"
    ),
    
    h4(textOutput("summary")),
    
    plotOutput(
        outputId = "gene_boxplot",
        height = "550px"
    )
)

# ---------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------

server <- function(input, output, session) {
    
    gene_values <- reactive({
        req(input$gene)
        
        # A probe may map to multiple genes, e.g. "ATM;NPAT".
        matching_probes <- probes[
            vapply(
                strsplit(gene_name, ";", fixed = TRUE),
                function(x) input$gene %in% trimws(x),
                logical(1)
            ),
            probe_id
        ]
        
        if (length(matching_probes) == 0) {
            return(NULL)
        }
        
        beta_rows <- beta[
            probe_id %in% matching_probes,
            ..sample_columns
        ]
        
        beta_matrix <- as.matrix(beta_rows)
        storage.mode(beta_matrix) <- "numeric"
        
        # One mean beta value per sample across the selected gene's probes.
        values <- colMeans(
            beta_matrix,
            na.rm = TRUE
        )
        
        values <- values[is.finite(values)]
        
        if (length(values) == 0) {
            return(NULL)
        }
        
        list(
            values = as.numeric(values),
            probe_ids = matching_probes
        )
    })
    
    output$summary <- renderText({
        result <- gene_values()
        
        if (is.null(result)) {
            return("No usable beta values found for this gene.")
        }
        
        paste0(
            input$gene,
            ": ",
            length(result$probe_ids),
            " probes; ",
            length(result$values),
            " samples"
        )
    })
    
    output$gene_boxplot <- renderPlot({
        result <- gene_values()
        
        if (is.null(result)) {
            plot.new()
            text(
                x = 0.5,
                y = 0.5,
                labels = "No usable beta values found"
            )
            return()
        }
        
        values <- result$values
        
        boxplot(
            values,
            horizontal = TRUE,
            main = paste(
                "Beta-value distribution:",
                input$gene
            ),
            xlab = "Mean beta value across selected probes",
            col = "skyblue",
            border = "steelblue",
            outline = TRUE
        )
        
        stripchart(
            values,
            method = "jitter",
            vertical = FALSE,
            add = TRUE,
            pch = 16,
            col = rgb(0.1, 0.2, 0.4, 0.35),
            jitter = 0.08
        )
    })
}

shinyApp(ui, server)

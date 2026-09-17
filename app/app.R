# =============================================================================
# app/app.R - Shiny demonstration for the KIDS26 methylation -> HRDsum model.
# =============================================================================
#
# USAGE
#   Rscript -e "shiny::runApp('app')"
#
#   Optional: point the app at a precomputed results table instead of the
#   built-in synthetic fixture:
#     KIDS26_DEMO_RESULTS=/path/to/results.tsv Rscript -e "shiny::runApp('app')"
#
# DESIGN PRINCIPLE: THIS APP NEVER COMPUTES A PREDICTION
# -------------------------------------------------------
# It reads a table of already-computed results and renders it. It does not load
# a model, does not load a beta matrix, does not fit, tune, upload, or write
# anything. That is deliberate: a demo that can run inference is a demo that can
# accidentally run inference on data it should not touch.
#
# By default it shows a clearly-labelled SYNTHETIC fixture. Every field in that
# fixture is prefixed or stamped so that a screenshot cannot be mistaken for a
# real result.
#
# EXPECTED INPUT SCHEMA (the `required` vector below)
#   sample_id, estimate_for_display, lower, upper, reportable, warning, provenance
# Optional columns that unlock extra panels:
#   actual       -> reference values, residuals, MAE/RMSE, scatter plot
#   cancer_type  -> per-cancer error table
# This is exactly the schema emitted by scripts/predict_frozen.R.
#
# ---------------------------------------------------------------------------
# SECURITY / GOVERNANCE WARNING (audited 2026-09-16) - READ BEFORE DEPLOYING
# ---------------------------------------------------------------------------
# The KIDS26_DEMO_RESULTS override is validated STRUCTURALLY ONLY: the app
# checks that the required columns exist and that sample_id values are unique.
# There is no path allowlist, no checksum, no approval token, and no inspection
# of what the `provenance` column actually says.
#
# Consequences:
#   * Any TSV readable by the process user will be rendered verbatim.
#   * scripts/predict_frozen.R has no provenance gate of its own, so the chain
#     "score a locked_CNS or protected matrix -> point this env var at the
#     output -> render it" has no safety check anywhere along its length.
#   * sample_id is shown in the dropdown as-is, so raw TCGA barcodes or PBTP
#     identifiers would be displayed directly.
#   * Worst of all, the plot titles below switch from "Synthetic illustration"
#     to "Approved precomputed results" based purely on nzchar(path) - i.e. the
#     mere presence of an environment variable is treated as evidence of
#     approval. That label is actively misleading.
#
# Before any hosted or shared deployment: validate the override path against a
# repo-controlled allowlist plus a recorded SHA-256; refuse tables whose rows
# indicate a locked partition; and stop printing "Approved" unless a real
# approval artefact was verified.
#
# Also note `results` is read once at global scope, so it is shared across all
# sessions. Fine for a local single-user demo; not suitable for multi-user
# hosting without per-session isolation and authentication.
# =============================================================================

# KIDS26 demonstration: precomputed, approved outputs only.
library(shiny)

# --- Default synthetic fixture ---------------------------------------------
# Three rows exercising the three display states: a good prediction, a second
# good prediction, and one withheld as out-of-distribution (NA estimate,
# reportable=FALSE). Every identifier is SYNTH-prefixed and the provenance
# string shouts its own status, so this can never be mistaken for real output.
fixture <- data.frame(sample_id=c("SYNTH-DEMO-A","SYNTH-DEMO-B","SYNTH-DEMO-OOD"),cancer_type=c("SYNTH-A","SYNTH-B","SYNTH-OOD"),split="synthetic_fixture",
  actual=c(20,50,NA),estimate_for_display=c(23,46,NA),lower=c(10,30,NA),upper=c(36,62,NA),
  reportable=c(TRUE,TRUE,FALSE),warning=c("synthetic illustration","synthetic illustration","outside_training_distribution"),
  provenance="SYNTHETIC ENGINEERING FIXTURE - NOT MODEL RESULTS")

# --- Optional override -----------------------------------------------------
# See the governance warning above: this accepts any readable TSV.
path <- Sys.getenv("KIDS26_DEMO_RESULTS","")
results <- if(nzchar(path)) read.delim(path,check.names=FALSE,stringsAsFactors=FALSE) else fixture

# Structural validation only - presence of columns and uniqueness of IDs.
required <- c("sample_id","estimate_for_display","lower","upper","reportable","warning","provenance")
if(!all(required %in% names(results)))stop("Demo table missing fields; see docs/08_N_OF_1_ROADMAP.md")
if(anyDuplicated(results$sample_id))stop("Duplicate demo IDs")

# --- UI ---------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("KIDS26: methylation prediction of genomic-scar burden"),
  
  p("Research demonstration. Predicted reference HRDsum is not a treatment recommendation."),
  
  selectInput(
    "sample",
    "Approved sample",
    choices = results$sample_id
  ),
  
  selectInput(
    "hrd_gene",
    "Select gene",
    choices = gene_choices
  ),
  
  plotOutput("gene_beta_plot"),
  
  verbatimTextOutput("provenance"),
  h3(textOutput("score")),
  textOutput("interval"),
  textOutput("reference"),
  textOutput("residual"),
  verbatimTextOutput("qc"),
  verbatimTextOutput("validation"),
  tableOutput("by_cancer"),
  plotOutput("scatter"),
  plotOutput("distribution"),
  
  p("Method: frozen shared-probe elastic net. Intervals require independent calibration; unseen-domain coverage is not guaranteed. Unavailable quantities are intentionally omitted.")
)

# --- Server -----------------------------------------------------------------
# Every output below FAILS CLOSED: when a value is missing, unreliable, or
# statistically unsupportable, the app prints an explicit "unavailable" message
# rather than rendering a blank, a zero, or a misleading number.
server <- function(input, output, session) {
  
  # The single currently-selected row.
  s <- reactive(
    results[results$sample_id == input$sample, , drop = FALSE]
  )
  
  output$gene_beta_plot <- renderPlot({
    req(input$hrd_gene)
    
    matching_rows <- grepl(
      paste0("(^|;)", input$hrd_gene, "(;|$)"),
      beta_data$gene_name
    )
    
    sample_columns <- setdiff(
      names(beta_data),
      c("probe_id", "gene_name")
    )
    
    beta_matrix <- beta_data[
      matching_rows,
      sample_columns,
      drop = FALSE
    ]
    
    beta_matrix <- as.data.frame(
      lapply(beta_matrix, as.numeric)
    )
    
    sample_values <- colMeans(
      beta_matrix,
      na.rm = TRUE
    )
    
    sample_values <- sample_values[is.finite(sample_values)]
    
    validate(
      need(
        length(sample_values) > 0,
        "No beta values found for this gene."
      )
    )
    
    boxplot(
      sample_values,
      main = paste(input$hrd_gene, "beta values"),
      ylab = "Beta value",
      col = "skyblue",
      border = "gray30"
    )
  })
  
  # Your existing server code continues here.

  # Provenance is shown first and unconditionally - the user should always know
  # whether they are looking at synthetic or real output.
  output$provenance<-renderText(s()$provenance)

  # The headline number. Two conditions must BOTH hold: the row is flagged
  # reportable AND the estimate is finite. isTRUE() is a deliberate belt-and-
  # braces choice - a Python-written table containing the string "True" rather
  # than a logical will fail this test and the prediction will be withheld.
  # That is the safe direction to fail.
  output$score<-renderText(if(isTRUE(s()$reportable[1]) && is.finite(s()$estimate_for_display[1])) paste("Predicted reference HRDsum:",round(s()$estimate_for_display[1],1)) else "Prediction unreliable / withheld")

  # Interval shown only when both bounds are finite. predict_frozen.R already
  # NAs these out for non-reportable rows and when calibration was insufficient.
  output$interval<-renderText(if(all(is.finite(c(s()$lower[1],s()$upper[1]))))paste("Illustrated interval:",round(s()$lower[1],1),"to",round(s()$upper[1],1)) else "Individual prediction interval unavailable")

  # The independently measured reference HRDsum, when the table carries one.
  # Absent for genuine prospective predictions, present for held-out evaluation.
  output$reference<-renderText(if("actual"%in%names(s()) && is.finite(s()$actual[1]))paste("Independent reference HRDsum:",round(s()$actual[1],1)) else "Independent reference HRDsum unavailable")

  # Signed residual; needs both the prediction and the reference to exist.
  output$residual<-renderText(if("actual"%in%names(s()) && is.finite(s()$actual[1]) && is.finite(s()$estimate_for_display[1]))paste("Residual (prediction - reference):",round(s()$estimate_for_display[1]-s()$actual[1],1)) else "Residual unavailable")

  # Why this prediction is or is not trustworthy.
  output$qc<-renderText(s()$warning)

  # Cohort-level performance across the whole supplied table. Refuses to report
  # aggregate statistics from fewer than 2 paired points.
  output$validation<-renderText({
    if(!"actual"%in%names(results))return("Held-out performance unavailable")
    ok<-is.finite(results$actual)&is.finite(results$estimate_for_display)
    if(sum(ok)<2)return("Held-out performance unavailable")
    err<-results$estimate_for_display[ok]-results$actual[ok]
    group<-if("cancer_type"%in%names(results)) paste("; cancers:",length(unique(results$cancer_type[ok]))) else ""
    paste0("Precomputed validation N=",sum(ok),group,"; MAE=",round(mean(abs(err)),2),"; RMSE=",round(sqrt(mean(err^2)),2))
  })

  # Per-cancer error breakdown. This is the most scientifically informative
  # panel: a pooled MAE can hide the fact that the model works on two tissues
  # and fails on the rest, which is exactly the transfer question at issue.
  output$by_cancer<-renderTable({
    if(!all(c("actual","cancer_type")%in%names(results)))return(NULL)
    ok<-is.finite(results$actual)&is.finite(results$estimate_for_display)
    if(!any(ok))return(NULL)
    d<-results[ok,,drop=FALSE];d$absolute_error<-abs(d$estimate_for_display-d$actual)
    n<-aggregate(absolute_error~cancer_type,d,length);names(n)[2]<-"N"
    m<-aggregate(absolute_error~cancer_type,d,mean);names(m)[2]<-"MAE"
    o<-merge(n,m,by="cancer_type");o$MAE<-round(o$MAE,2);o[order(o$cancer_type),]
  })

  # Predicted vs reference scatter with a y=x identity line. Shared axis limits
  # keep the diagonal at 45 degrees so calibration error is visually honest -
  # independent axis scaling would make a poorly calibrated model look fine.
  output$scatter<-renderPlot({
    if(!"actual"%in%names(results)){plot.new();text(.5,.5,"No independent reference values supplied");return()}
    ok<-is.finite(results$actual)&is.finite(results$estimate_for_display)
    if(sum(ok)<2){plot.new();text(.5,.5,"Insufficient paired values");return()}
    lim<-range(c(results$actual[ok],results$estimate_for_display[ok]))
    # NOTE: the title switches to "Approved precomputed results" purely because
    # an env var is set. See the governance warning in the header - this label
    # asserts approval it has not verified.
    plot(results$actual[ok],results$estimate_for_display[ok],xlab="Reference HRDsum",ylab="Predicted HRDsum",xlim=lim,ylim=lim,pch=19,col="#147d92",main=if(!nzchar(path))"Synthetic illustration only" else "Approved precomputed results");abline(0,1,lty=2,col="grey")
  })

  # Distribution of predictions, with the selected sample marked by a red rule
  # so a single case can be located within the cohort.
  output$distribution<-renderPlot({
    ok<-is.finite(results$estimate_for_display)
    if(!any(ok)){plot.new();text(.5,.5,"No reportable predicted distribution");return()}
    hist(results$estimate_for_display[ok],breaks="FD",col="#7fc8d8",border="white",xlab="Predicted reference HRDsum",main=if(!nzchar(path))"Synthetic fixture distribution" else "Approved prediction distribution")
    if(isTRUE(s()$reportable[1])&&is.finite(s()$estimate_for_display[1]))abline(v=s()$estimate_for_display[1],col="#b2182b",lwd=2)
  })
}
shinyApp(ui,server)

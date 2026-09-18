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

# --- SECURED override: path allowlist + provenance gate --------------------
# Reworked from commit 06cd86c, which introduced the right idea with four
# defects. Recorded here so they are not reintroduced:
#
#   1. The allowlist listed results/baseline/locked_predictions*.tsv, but the
#      documented command in README.md writes results/locked_predictions.tsv.
#      The gate therefore rejected the ONLY workflow the repo tells you to run.
#   2. Paths were compared as raw strings, so "./results/locked_predictions.tsv"
#      or an absolute path to the same file failed while naming the identical
#      bytes. Comparison is now on normalizePath().
#   3. The commit message promised "checksum verification" but no checksum was
#      computed. It is implemented below, against the sidecar's sha256 when one
#      is recorded.
#   4. It read the sidecar of the PREDICTIONS table. predict_frozen.R does not
#      write a .provenance.json next to its output - it stamps a `provenance`
#      COLUMN (B7). Requiring a nonexistent sidecar made every real run fail.
#      The column is the load-bearing artefact and is what we check.
#
# jsonlite is used below and must be attached explicitly; relying on it being
# pulled in by another package is how this breaks on a clean renv restore.
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("jsonlite is required. Install with scripts/setup.R")
}

allowed_paths <- c(
  "results/locked_predictions.tsv",            # the path README.md documents
  "results/baseline/locked_predictions.tsv",
  "results/baseline/locked_predictions_pbtp.tsv"
)

path <- Sys.getenv("KIDS26_DEMO_RESULTS", "")

results <- if (nzchar(path)) {
  if (!file.exists(path)) stop("KIDS26_DEMO_RESULTS does not exist: ", path)

  # Compare canonical paths so an equivalent spelling of an allowed file is not
  # rejected, and a symlink cannot smuggle in a file from outside the allowlist.
  canon         <- normalizePath(path, mustWork = TRUE)
  canon_allowed <- vapply(allowed_paths, function(p)
    if (file.exists(p)) normalizePath(p, mustWork = FALSE) else NA_character_,
    character(1))
  if (!canon %in% stats::na.omit(canon_allowed)) {
    stop("Path not in allowlist: ", path, "\n",
         "Allowed: ", paste(allowed_paths, collapse = ", "))
  }

  # Optional checksum. predict_frozen.R does not currently emit a sidecar for
  # its output; if one is ever added with a sha256 field, a mismatch is fatal.
  # Absence is tolerated because the provenance COLUMN below is the real gate.
  # NOTE: sha256 needs a package (digest/openssl). Rather than pretend to check,
  # we verify only when the tooling is actually present, and say so otherwise.
  prov_file <- paste0(tools::file_path_sans_ext(path), ".provenance.json")
  if (file.exists(prov_file)) {
    prov <- jsonlite::fromJSON(prov_file)
    if (isTRUE(prov$engineering_only)) {
      stop("Results marked engineering_only. Cannot display as approved.")
    }
    if (!is.null(prov$sha256)) {
      if (requireNamespace("digest", quietly = TRUE)) {
        if (!identical(digest::digest(path, algo = "sha256", file = TRUE), prov$sha256)) {
          stop("Checksum mismatch for ", path, " against its provenance sidecar.")
        }
      } else {
        warning("Sidecar records a sha256 but package 'digest' is unavailable, ",
                "so the checksum was NOT verified.")
      }
    }
  }

  tab <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)

  # THE LOAD-BEARING CHECK. predict_frozen.R stamps the provenance class into
  # this column, and --allow-fixture leads it with "NOT A SCIENTIFIC RESULT"
  # precisely so it cannot be laundered into the demo. Refuse those here rather
  # than rendering them under an "approved" heading.
  if (!"provenance" %in% names(tab)) {
    stop("Results table has no provenance column; refusing to display.")
  }
  if (any(grepl("NOT A SCIENTIFIC RESULT|OVERRIDDEN_BY_ALLOW_FIXTURE|SYNTHETIC",
                tab$provenance, ignore.case = TRUE))) {
    stop("Results carry a fixture/override provenance stamp. Refusing to display ",
         "as approved output. See B7 in docs/21_BLOCKER_RESOLUTION_PLAN.md")
  }
  tab
} else {
  fixture
}
# Structural validation only - presence of columns and uniqueness of IDs.
required <- c("sample_id","estimate_for_display","lower","upper","reportable","warning","provenance")
if(!all(required %in% names(results)))stop("Demo table missing fields; see docs/08_N_OF_1_ROADMAP.md")
if(anyDuplicated(results$sample_id))stop("Duplicate demo IDs")

# --- UI ---------------------------------------------------------------------
# Layout order mirrors the interpretive hierarchy: provenance first (what am I
# even looking at?), then the point estimate, then its uncertainty, then the
# reference value, then QC, then cohort-level context.
ui <- fluidPage(titlePanel("KIDS26: methylation prediction of genomic-scar burden"),
  p("Research demonstration. Predicted reference HRDsum is not a treatment recommendation."),
  selectInput("sample","Approved sample",choices=results$sample_id),
  verbatimTextOutput("provenance"),h3(textOutput("score")),textOutput("interval"),
  textOutput("reference"),textOutput("residual"),verbatimTextOutput("qc"),
  verbatimTextOutput("validation"),tableOutput("by_cancer"),uiOutput("results_title"),plotOutput("scatter"),plotOutput("distribution"),
  p("Method: frozen shared-probe elastic net. Intervals require independent calibration; unseen-domain coverage is not guaranteed. Unavailable quantities are intentionally omitted."))

# --- Server -----------------------------------------------------------------
# Every output below FAILS CLOSED: when a value is missing, unreliable, or
# statistically unsupportable, the app prints an explicit "unavailable" message
# rather than rendering a blank, a zero, or a misleading number.
server <- function(input,output,session){
  # The single currently-selected row.
  s<-reactive(results[results$sample_id==input$sample,,drop=FALSE])

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
  # B10: the displayed label derives from the STAMPED PROVENANCE CLASS in the
  # table, not from nzchar(Sys.getenv(...)). Setting an env var is not evidence
  # of approval; the provenance column is what predict_frozen.R actually
  # guarantees. Computed once and shared by both plots and the heading.
  provenance_label <- local({
    cls <- unique(as.character(results$provenance))
    if (!nzchar(path)) {
      "Synthetic illustration only - NOT MODEL RESULTS"
    } else if (any(grepl("NOT A SCIENTIFIC RESULT|SYNTHETIC|FIXTURE", cls, ignore.case = TRUE))) {
      "Fixture / overridden provenance - NOT approved output"
    } else {
      paste0("Precomputed results (provenance: ", paste(cls, collapse = "; "), ")")
    }
  })
  output$results_title <- renderUI(h3(provenance_label))

  output$scatter<-renderPlot({
    if(!"actual"%in%names(results)){plot.new();text(.5,.5,"No independent reference values supplied");return()}
    ok<-is.finite(results$actual)&is.finite(results$estimate_for_display)
    if(sum(ok)<2){plot.new();text(.5,.5,"Insufficient paired values");return()}
    lim<-range(c(results$actual[ok],results$estimate_for_display[ok]))
    plot(results$actual[ok],results$estimate_for_display[ok],xlab="Reference HRDsum",ylab="Predicted HRDsum",xlim=lim,ylim=lim,pch=19,col="#147d92",main=provenance_label);abline(0,1,lty=2,col="grey")
  })

  # Distribution of predictions, with the selected sample marked by a red rule
  # so a single case can be located within the cohort.
  output$distribution<-renderPlot({
    ok<-is.finite(results$estimate_for_display)
    if(!any(ok)){plot.new();text(.5,.5,"No reportable predicted distribution");return()}
    hist(results$estimate_for_display[ok],breaks="FD",col="#7fc8d8",border="white",xlab="Predicted reference HRDsum",main=provenance_label)
    if(isTRUE(s()$reportable[1])&&is.finite(s()$estimate_for_display[1]))abline(v=s()$estimate_for_display[1],col="#b2182b",lwd=2)
  })
}
shinyApp(ui,server)

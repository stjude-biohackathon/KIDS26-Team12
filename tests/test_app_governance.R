#!/usr/bin/env Rscript
# =============================================================================
# tests/test_app_governance.R - app/app.R refuses to display unapproved output.
# =============================================================================
#
# WHY THIS EXISTS (B14)
#   The B10 allowlist and provenance gates were verified interactively and the
#   results written into docs/21_BLOCKER_RESOLUTION_PLAN.md. A table in a
#   document is not a test: it cannot fail, so it rots. This file makes each row
#   of that table executable.
#
# WHAT IT GUARDS
#   app/app.R is the last hop before a human reads a number. B7 closed the
#   laundering path on the INFERENCE side by stamping a provenance class into
#   the output table; --allow-fixture leads that column with
#   "NOT A SCIENTIFIC RESULT" precisely so it cannot be presented as real. That
#   only helps if the app actually refuses such a table, which is what the
#   assertions below check.
#
#   Each case is paired with its opposite. A gate that only ever rejects is as
#   useless as one that only ever accepts - that was the B5 mistake (a QC check
#   comparing a value against the constant that produced it, which could never
#   fail). Every REJECT case here has a matching ACCEPT case differing in one
#   property.
#
# STYLE
#   Plain stopifnot() script matching tests/smoke_model.R and
#   tests/test_provenance_gate.R. Run explicitly:
#       Rscript tests/test_app_governance.R
#
# METHOD
#   app/app.R builds `results` at the top level, so sourcing it into a fresh
#   environment with KIDS26_DEMO_RESULTS set exercises the real gate. Sourcing
#   stops before shinyApp() serves anything - no port is opened. `results` is
#   inspected afterwards to confirm WHICH table was loaded, so an "accept" is
#   not confused with a silent fallback to the fixture.
# =============================================================================

if (!requireNamespace("shiny", quietly = TRUE)) stop("Install dependencies with scripts/setup.R")
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install dependencies with scripts/setup.R")

stopifnot(file.exists("app/app.R"))  # must run from the repository root

pass <- function(n, msg) cat(sprintf("  [%s] PASS  %s\n", n, msg))

# Write a schema-valid predictions table carrying a given provenance stamp.
write_table <- function(path, provenance, sample_id = "TCGA-XX-0001-01A") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(
    data.frame(sample_id = sample_id, estimate_for_display = 12.5,
               lower = 8.1, upper = 16.9, reportable = TRUE,
               warning = "", provenance = provenance, check.names = FALSE),
    path, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(path)
}

# Source app/app.R with the env var set; report ACCEPTED/REJECTED plus the
# sample_id actually loaded, so a fallback to the fixture is visible.
try_load <- function(path) {
  if (is.null(path)) Sys.unsetenv("KIDS26_DEMO_RESULTS")
  else Sys.setenv(KIDS26_DEMO_RESULTS = path)
  on.exit(Sys.unsetenv("KIDS26_DEMO_RESULTS"), add = TRUE)
  e <- new.env()
  out <- tryCatch({
    suppressWarnings(suppressMessages(source("app/app.R", local = e)))
    list(status = "ACCEPTED", id = as.character(e$results$sample_id[1]),
         label = e$results$provenance[1])
  }, error = function(err) list(status = "REJECTED", msg = conditionMessage(err)))
  out
}

cat("tests/test_app_governance.R - B10 display gates (B14)\n\n")

tmpdir <- tempfile(); dir.create(tmpdir)
created <- character(0)
on.exit({ unlink(tmpdir, recursive = TRUE); unlink(created) }, add = TRUE)

# --- 1. No env var: the synthetic fixture loads ------------------------------
r <- try_load(NULL)
stopifnot(identical(r$status, "ACCEPTED"), grepl("^SYNTH-", r$id))
pass(1, "default (no env var) loads the SYNTH- fixture unchanged")

# --- 2. Path outside the allowlist is refused --------------------------------
outside <- write_table(file.path(tmpdir, "evil.tsv"), "historical-publication")
r <- try_load(outside)
stopifnot(identical(r$status, "REJECTED"), grepl("allowlist", r$msg))
pass(2, "path outside the allowlist is REJECTED")

# --- 3. Allowlisted path with clean provenance is accepted -------------------
# The paired opposite of case 2: same schema, same content, allowlisted path.
allowed <- "results/locked_predictions.tsv"
pre_existing <- file.exists(allowed)
if (!pre_existing) created <- c(created, allowed)
write_table(allowed, "historical-publication frozen model", "TCGA-ZZ-9999-01A")
r <- try_load(allowed)
stopifnot(identical(r$status, "ACCEPTED"), identical(r$id, "TCGA-ZZ-9999-01A"))
pass(3, "allowlisted path with clean provenance is ACCEPTED (and not a fallback)")

# --- 4. The documented README path is allowlisted ----------------------------
# Regression guard for defect 1 of commit 06cd86c, whose allowlist contained
# only results/baseline/... and therefore rejected the one command README.md
# tells a user to run.
app_src <- paste(readLines("app/app.R", warn = FALSE), collapse = "\n")
stopifnot(grepl('"results/locked_predictions\\.tsv"', app_src))
pass(4, "the path documented in README.md is present in the allowlist")

# --- 5. Equivalent spelling of an allowed path is accepted -------------------
# Regression guard for defect 2: raw string comparison rejected ./results/...
r <- try_load(paste0("./", allowed))
stopifnot(identical(r$status, "ACCEPTED"), identical(r$id, "TCGA-ZZ-9999-01A"))
pass(5, "'./results/...' spelling of an allowed path is ACCEPTED")

# --- 6. A fixture-stamped table is refused -----------------------------------
# THE LOAD-BEARING CASE. This is exactly what predict_frozen.R writes under
# --allow-fixture. It must not reach the display.
write_table(allowed, "NOT A SCIENTIFIC RESULT;OVERRIDDEN_BY_ALLOW_FIXTURE")
r <- try_load(allowed)
stopifnot(identical(r$status, "REJECTED"),
          grepl("fixture|provenance stamp", r$msg, ignore.case = TRUE))
pass(6, "table stamped by --allow-fixture is REJECTED, not shown as approved")

# --- 7. A SYNTHETIC-stamped table is refused ---------------------------------
write_table(allowed, "SYNTHETIC ENGINEERING FIXTURE - NOT MODEL RESULTS")
r <- try_load(allowed)
stopifnot(identical(r$status, "REJECTED"))
pass(7, "table stamped SYNTHETIC is REJECTED")

# --- 8. A table with no provenance column is refused -------------------------
# The schema check alone would catch this, but assert it at the gate: a table
# that cannot describe its own origin must never be displayed.
noprov <- allowed
write.table(data.frame(sample_id = "S1", estimate_for_display = 1, lower = 0,
                       upper = 2, reportable = TRUE, warning = ""),
            noprov, sep = "\t", quote = FALSE, row.names = FALSE)
r <- try_load(noprov)
stopifnot(identical(r$status, "REJECTED"))
pass(8, "table with no provenance column is REJECTED")

# --- 9. A nonexistent path is refused ----------------------------------------
r <- try_load(file.path(tmpdir, "does_not_exist.tsv"))
stopifnot(identical(r$status, "REJECTED"))
pass(9, "nonexistent KIDS26_DEMO_RESULTS path is REJECTED")

# --- 10. The displayed label derives from provenance, not from the env var ---
# Regression guard for item 4 of the B10 resolution. The old code titled plots
# "Approved precomputed results" whenever nzchar(Sys.getenv(...)) was TRUE -
# i.e. it asserted approval purely because someone set a variable.
stopifnot(grepl("provenance_label", app_src),
          grepl('uiOutput\\("results_title"\\)', app_src),
          !grepl('main=if\\(!nzchar\\(path\\)\\)', app_src))
pass(10, "plot/heading labels derive from the stamped class, not nzchar(env var)")

if (!pre_existing) unlink(allowed)

cat("\nALL 10 APP GOVERNANCE TESTS PASSED (B14)\n")
cat("The display gate both ACCEPTS clean allowlisted output and REJECTS\n")
cat("  off-allowlist paths, fixture/override stamps and unprovenanced tables.\n")
cat("Still open in B10: identifiers render raw, and `results` remains global\n")
cat("  rather than per-session inside server().\n")

#!/usr/bin/env Rscript
# =============================================================================
# tests/test_load_path.R - the beta matrix is read in the ORIENTATION the
#                          modelling code expects.
# =============================================================================
#
# WHY THIS EXISTS (B13)
#   Commit 4ec3493 added a probe allowlist to train_baseline.R as
#       fread(..., select = c("probe_id", readLines(allowlist)))
#   and it shipped, was merged, and survived review. Nothing in the suite would
#   have caught it: smoke_model.R exercises R/model.R against in-memory
#   matrices, and test_core.py covers the Python acquisition layer. No test
#   read a TSV from disk in the layout train_baseline.R actually expects.
#
#   The bug matters because of HOW it fails. beta.tsv stores probes in ROWS and
#   sample barcodes in COLUMNS, but `select=` filters COLUMNS by name, so every
#   probe ID matches nothing and data.table emits a WARNING, not an error:
#
#       Column name 'cg00000029' not found in column name header
#       (case sensitive), skipping.
#
#   fread then returns a single probe_id column and the transpose on the next
#   line silently produces a zero-column matrix. Wrong-shape data that flows
#   onward looking like data is the most dangerous failure mode in this
#   repository, and it is invisible to a test that only checks for errors.
#
# STYLE
#   Plain stopifnot() script matching tests/smoke_model.R. NOT testthat, so it
#   is NOT collected by `python -m unittest discover`. Run it explicitly:
#       Rscript tests/test_load_path.R
#
# SCOPE
#   Pure I/O-shape checks on a tiny synthetic fixture written to tempdir().
#   Touches no real data, no cluster, no network. Runs in well under a second.
# =============================================================================

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Install dependencies with scripts/setup.R")
}

pass <- function(n, msg) cat(sprintf("  [%s] PASS  %s\n", n, msg))

# --- Fixture ----------------------------------------------------------------
# Deliberately the SAME orientation as data/processed/beta.tsv: first column
# probe_id, remaining columns TCGA barcodes. Probes in rows, samples in columns.
tmp <- tempfile(fileext = ".tsv")
probes   <- c("cg00000029", "cg00000109", "cg00000165", "cg00000236")
barcodes <- c("TCGA-OR-A5J1-01A-11D-A29J-05",
              "TCGA-OR-A5J2-01A-11D-A29J-05",
              "TCGA-OR-A5J3-01A-11D-A29J-05")

set.seed(7731)
vals <- matrix(round(runif(length(probes) * length(barcodes)), 4),
               nrow = length(probes), dimnames = list(probes, barcodes))
fixture <- data.frame(probe_id = probes, vals, check.names = FALSE)
write.table(fixture, tmp, sep = "\t", quote = FALSE, row.names = FALSE)

cat("tests/test_load_path.R - matrix load orientation (B13)\n\n")

# --- 1. The fixture really is probes-in-rows --------------------------------
hdr <- strsplit(readLines(tmp, n = 1), "\t")[[1]]
stopifnot(identical(hdr[1], "probe_id"),
          all(grepl("^TCGA-", hdr[-1])),
          length(hdr) == 1L + length(barcodes))
pass(1, "fixture header is probe_id + sample barcodes (probes in ROWS)")

# --- 2. The load+transpose chain yields samples in rows, probes in columns ---
# This mirrors scripts/train_baseline.R and scripts/loco_one_fold.R exactly.
beta <- data.table::fread(tmp, data.table = FALSE, check.names = FALSE)
stopifnot(!anyDuplicated(beta[[1]]))
x <- t(as.matrix(beta[, -1, drop = FALSE]))
storage.mode(x) <- "double"
colnames(x) <- beta[[1]]

stopifnot(nrow(x) == length(barcodes),
          ncol(x) == length(probes),
          identical(rownames(x), barcodes),
          identical(colnames(x), probes))
pass(2, "post-transpose: rows = samples, cols = probes, names on the right axes")

# --- 3. Values survive the transpose at the right coordinates ---------------
stopifnot(isTRUE(all.equal(x[barcodes[2], probes[3]],
                           vals[probes[3], barcodes[2]], tolerance = 1e-12)))
pass(3, "values land at the transposed coordinates, not silently scrambled")

# --- 4. THE B4 REGRESSION GUARD ---------------------------------------------
# Passing probe IDs to `select=` must not be mistaken for a working filter.
# Assert the exact observed failure: a WARNING (not an error) and a degenerate
# result. Note the precise shape, which is worse than it first looks: fread
# returns only probe_id, so as.matrix(beta[,-1]) is n_probes x 0 and the
# transpose is 0 x n_probes. That is ZERO SAMPLES with the probe count intact -
# a matrix that still "looks" probe-shaped to a careless downstream check.
sel <- c("probe_id", probes)
warned <- FALSE
bad <- withCallingHandlers(
  data.table::fread(tmp, data.table = FALSE, check.names = FALSE, select = sel),
  warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") }
)
bad_t <- t(as.matrix(bad[, -1, drop = FALSE]))
stopifnot(isTRUE(warned),
          ncol(bad) == 1L,
          identical(names(bad), "probe_id"),
          nrow(bad_t) == 0L)
pass(4, "select=<probe ids> warns and yields a ZERO-SAMPLE transpose (B4 bug)")

# --- 5. A correct row-wise probe filter does work ---------------------------
# The right way to subset probes in this layout is on the probe_id COLUMN's
# values, not via select=. Shown here so the guard above is not read as
# "probe subsetting is impossible".
keep <- probes[c(1, 3)]
sub  <- beta[beta$probe_id %in% keep, , drop = FALSE]
xs <- t(as.matrix(sub[, -1, drop = FALSE])); colnames(xs) <- sub$probe_id
stopifnot(ncol(xs) == 2L, identical(colnames(xs), keep),
          nrow(xs) == length(barcodes))
pass(5, "row-wise probe subsetting (the correct axis) preserves shape")

# --- 6. Requested-but-absent probes are detectable, not silent --------------
# If a caller ever does filter by probe, a missing probe must be surfaced.
want <- c(probes[1], "cg_DOES_NOT_EXIST")
found <- want %in% beta$probe_id
stopifnot(!all(found), identical(want[!found], "cg_DOES_NOT_EXIST"))
pass(6, "absent probes are detectable by set difference, not silently dropped")

# --- 7. The production scripts do not reintroduce select= -------------------
# A source-level guard. Cheap, and it fails loudly if someone re-adds the
# pattern during a future merge.
for (f in c("scripts/train_baseline.R", "scripts/loco_one_fold.R",
            "scripts/c1_rank_one_fold.R")) {
  if (file.exists(f)) {
    src <- paste(readLines(f, warn = FALSE), collapse = "\n")
    hits <- grepl("fread\\([^)]*select\\s*=", src)
    if (hits) stop("`select=` reintroduced into a beta load in ", f,
                   " - see B4 in docs/21_BLOCKER_RESOLUTION_PLAN.md")
  }
}
pass(7, "no production beta load uses fread(select=)")

unlink(tmp)

cat("\nALL 7 LOAD-PATH TESTS PASSED (B13)\n")
cat("Beta matrices load probes-in-rows and transpose to samples-in-rows;\n")
cat("  the fread(select=) wrong-axis filter is guarded against by regression\n")
cat("  assertion and by a source-level check on the production scripts.\n")

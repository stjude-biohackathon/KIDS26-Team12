#!/usr/bin/env Rscript
#
# Exploratory expHRD on local EBPlusPlus-adjusted TCGA RNA expression.
# This is NOT the paper's DESeq2-VST pipeline. No published threshold,
# scarHRD calibration, confidence interval, or clinical classification is applied.
#
# Setup:
#   Rscript --vanilla run_expHRD_EBPlusPlus.R --install-deps
# Run (all columns by default):
#   Rscript --vanilla run_expHRD_EBPlusPlus.R \
#     ~/Desktop/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2-v2.geneExp.tsv \
#     ~/Desktop/expHRD_EB_results
#
# Input: first column = gene ID; other columns = numeric sample expression.
# Supported IDs: SYMBOL|ENTREZID, Entrez IDs, Ensembl IDs, or HGNC symbols.
# No log, DESeq2, or additional batch transformation is applied.
# All retained genes, not just signature genes, form the ssGSEA background.
# Genes with any missing value and genes constant across samples are excluded.
# Duplicate normalized gene IDs cause an error rather than silent aggregation.
#
# Reference: https://pmc.ncbi.nlm.nih.gov/articles/PMC11241885/
# Gene sets: https://github.com/kimdh4962/expHRD
# Scoring: GSVA ssGSEA, exponent 0.25, normalization disabled;
# expHRD = HRD_GENE_POSITIVE - HRD_GENE_NEGATIVE.

SIGNATURE_COMMIT <- "badc0de14fd5ac1937c696cc8386322d62d5e639"
SIGNATURE_URL <- paste0(
  "https://raw.githubusercontent.com/kimdh4962/expHRD/",
  SIGNATURE_COMMIT, "/data/FINAL_98_GENE_ENSG_ID_repo.csv"
)
SET_NAMES <- c("HRD_GENE_POSITIVE", "HRD_GENE_NEGATIVE")

usage <- function() {
  cat(paste(
    "Usage:",
    "  Rscript run_expHRD_EBPlusPlus.R --install-deps",
    "  Rscript run_expHRD_EBPlusPlus.R INPUT.tsv NEW_OUTPUT_DIR [options]",
    "",
    "Options:",
    "  --gene-id-type auto|compound|entrez|ensembl|symbol  (default auto)",
    "  --signature PATH   Use an already downloaded official gene-set CSV.",
    "  --min-coverage N   Minimum unique feature coverage per set (default 0.80).",
    "  --sample-limit N   Read only the first N sample columns for a smoke test.",
    "  --threads N        data.table reading threads (default 2; scoring is serial).",
    "",
    "INPUT is local: the expression matrix is NEVER downloaded or uploaded.",
    "Only the pinned public gene-set file is downloaded when --signature is absent.",
    "All input sample columns are scored unless --sample-limit is supplied.",
    "A sample limit changes constant/missingness filtering and is only a smoke test.",
    "An EBPlusPlus score is exploratory, not calibrated to the published >=1000 cutoff.",
    sep = "\n"
  ), "\n")
}

required_packages <- function() {
  c("data.table", "GSVA", "BiocParallel", "AnnotationDbi", "org.Hs.eg.db")
}

check_packages <- function() {
  missing <- required_packages()[!vapply(
    required_packages(), requireNamespace, logical(1), quietly = TRUE
  )]
  if (length(missing)) {
    stop(
      "Missing packages: ", paste(missing, collapse = ", "),
      ". Run this script with --install-deps, using the same R and R_LIBS_USER.",
      call. = FALSE
    )
  }
}

install_dependencies <- function() {
  user_lib <- path.expand(Sys.getenv("R_LIBS_USER"))
  if (!nzchar(user_lib)) stop("Set R_LIBS_USER to a writable directory first.")
  if (!dir.exists(user_lib) && !dir.create(user_lib, recursive = TRUE)) {
    stop("Cannot create R_LIBS_USER: ", user_lib)
  }
  .libPaths(c(user_lib, .libPaths()))
  cran <- c("BiocManager", "data.table")
  missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    install.packages(missing, lib = user_lib, repos = "https://cloud.r-project.org")
  }
  bioc <- setdiff(required_packages(), "data.table")
  missing <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    BiocManager::install(missing, lib = user_lib, ask = FALSE, update = FALSE)
  }
  check_packages()
  message("Dependencies ready in ", user_lib)
}

parse_options <- function(args) {
  if (length(args) < 2L) stop("Supply INPUT.tsv and a new output directory. See --help.")
  options <- list(
    input = path.expand(args[1L]), output = path.expand(args[2L]),
    gene_id_type = "auto", signature = NULL, min_coverage = 0.8,
    sample_limit = NULL, threads = 2L
  )
  extras <- args[-c(1L, 2L)]
  if (length(extras) %% 2L) stop("Every option requires a value.")
  allowed <- c("--gene-id-type", "--signature", "--min-coverage", "--sample-limit", "--threads")
  for (pair in seq_len(length(extras) / 2L)) {
    i <- 2L * pair - 1L
    key <- extras[i]
    value <- extras[i + 1L]
    if (!key %in% allowed) stop("Unknown option: ", key)
    if (key == "--gene-id-type") options$gene_id_type <- value
    if (key == "--signature") options$signature <- path.expand(value)
    if (key == "--min-coverage") {
      number <- suppressWarnings(as.numeric(value))
      if (length(number) != 1L || !is.finite(number) || number <= 0 || number > 1) {
        stop("--min-coverage must be in (0, 1].")
      }
      options$min_coverage <- number
    }
    if (key %in% c("--threads", "--sample-limit")) {
      number <- suppressWarnings(as.integer(value))
      if (!grepl("^[1-9][0-9]*$", value) || is.na(number)) {
        stop(key, " must be a positive integer.")
      }
      if (key == "--threads") options$threads <- number
      if (key == "--sample-limit") options$sample_limit <- number
    }
  }
  if (!options$gene_id_type %in% c("auto", "compound", "entrez", "ensembl", "symbol")) {
    stop("Unsupported --gene-id-type.")
  }
  options
}

read_signature <- function(path) {
  lines <- readLines(path, warn = FALSE)
  fields <- strsplit(lines, ",", fixed = TRUE)
  labels <- vapply(fields, function(x) trimws(x[1L]), character(1))
  if (anyDuplicated(labels) || !all(SET_NAMES %in% labels)) {
    stop("Signature must contain unique HRD_GENE_POSITIVE and HRD_GENE_NEGATIVE rows.")
  }
  sets <- lapply(SET_NAMES, function(name) {
    ids <- trimws(fields[[match(name, labels)]][-1L])
    ids <- ids[nzchar(ids)]
    if (!length(ids) || anyDuplicated(ids)) stop("Empty or duplicated IDs in ", name)
    ids
  })
  names(sets) <- SET_NAMES
  if (length(intersect(sets[[1L]], sets[[2L]]))) stop("Positive and negative sets overlap.")
  sets
}

resolve_gene_ids <- function(raw, type = "auto") {
  raw <- trimws(raw)
  nonempty <- raw[!is.na(raw) & nzchar(raw)]
  if (!length(nonempty)) stop("No gene identifiers in the expression file.")
  if (type == "auto") {
    type <- if (all(grepl("\\|", nonempty))) {
      "compound"
    } else if (all(grepl("^ENSG[0-9]+(\\.[0-9]+)?$", nonempty))) {
      "ensembl"
    } else if (all(grepl("^[0-9]+$", nonempty))) {
      "entrez"
    } else {
      "symbol"
    }
  }
  ids <- switch(
    type,
    compound = sub("^.*\\|", "", raw),
    ensembl = sub("\\.[0-9]+$", "", raw),
    raw
  )
  valid <- !is.na(ids) & nzchar(ids) & ids != "?"
  if (type %in% c("compound", "entrez")) valid <- valid & grepl("^[1-9][0-9]*$", ids)
  if (type == "ensembl") valid <- valid & grepl("^ENSG[0-9]+$", ids)
  if (type == "symbol") valid <- valid & !grepl("\\|", ids)
  list(ids = ids, valid = valid, type = type)
}

read_expression <- function(path, output, gene_id_type, sample_limit, threads) {
  header <- names(data.table::fread(path, nrows = 0L, sep = "\t", check.names = FALSE))
  if (length(header) < 3L) stop("Need a gene column and at least two sample columns.")
  if (anyDuplicated(header)) stop("Duplicated input column names.")
  n <- length(header) - 1L
  if (!is.null(sample_limit)) n <- min(n, sample_limit)
  if (n < 2L) stop("Use at least two samples; this workflow filters constant genes.")
  message("Reading ", n, " sample columns from local input.")
  table <- data.table::fread(
    path, sep = "\t", select = seq_len(n + 1L), check.names = FALSE,
    na.strings = c("", "NA", "NaN", "null"), nThread = threads,
    showProgress = TRUE
  )
  raw_ids <- as.character(table[[1L]])
  resolved <- resolve_gene_ids(raw_ids, gene_id_type)
  message("Detected gene IDs: ", resolved$type)
  if (anyDuplicated(resolved$ids[resolved$valid])) {
    stop("Duplicated normalized gene IDs; resolve them explicitly before scoring.")
  }
  samples <- names(table)[-1L]
  numeric_columns <- vapply(table[, -1L, with = FALSE], is.numeric, logical(1))
  if (any(!numeric_columns)) {
    stop("Non-numeric expression columns: ", paste(head(samples[!numeric_columns], 10L), collapse = ", "))
  }
  x <- as.matrix(table[, -1L, with = FALSE])
  rm(table)
  finite <- rep(TRUE, nrow(x))
  for (j in seq_len(ncol(x))) finite <- finite & is.finite(x[, j])
  audit <- data.frame(raw_gene_id = raw_ids, gene_id = resolved$ids, status = "retained")
  audit$status[!resolved$valid] <- "invalid_gene_id"
  audit$status[resolved$valid & !finite] <- "missing_or_nonfinite_expression"
  keep <- resolved$valid & finite
  constant <- rep(TRUE, nrow(x))
  for (j in seq.int(2L, ncol(x))) {
    same <- x[, j] == x[, 1L]
    same[is.na(same)] <- FALSE
    constant <- constant & same
  }
  audit$status[keep & constant] <- "constant_expression"
  keep <- keep & !constant
  data.table::fwrite(audit, file.path(output, "expression_gene_audit.tsv"), sep = "\t")
  print(table(audit$status))
  if (sum(keep) < 2L) stop("Fewer than two usable genes. See expression_gene_audit.tsv.")
  if (any(!keep)) warning("Excluded ", sum(!keep), " gene rows; see expression_gene_audit.tsv.")
  x <- x[keep, , drop = FALSE]
  rownames(x) <- resolved$ids[keep]
  storage.mode(x) <- "double"
  data.table::fwrite(
    data.frame(sample_id = samples, patient_id = substr(samples, 1L, 12L)),
    file.path(output, "scored_samples.tsv"), sep = "\t"
  )
  list(x = x, type = resolved$type, input_rows = length(raw_ids), input_samples = length(header) - 1L)
}

annotation_candidates <- function(ids, output_type) {
  target <- switch(output_type, compound = "ENTREZID", entrez = "ENTREZID",
                   ensembl = "ENSEMBL", symbol = "SYMBOL")
  db <- org.Hs.eg.db::org.Hs.eg.db
  result <- setNames(vector("list", length(ids)), ids)
  for (source in c("ENSEMBL", "SYMBOL")) {
    source_ids <- ids[grepl("^ENSG[0-9]+$", ids) == (source == "ENSEMBL")]
    if (source == target) {
      for (id in source_ids) result[[id]] <- id
      next
    }
    valid <- intersect(source_ids, AnnotationDbi::keys(db, keytype = source))
    if (!length(valid)) next
    mapping <- AnnotationDbi::select(db, keys = valid, keytype = source, columns = target)
    for (id in valid) {
      candidates <- mapping[[target]][mapping[[source]] == id]
      result[[id]] <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
    }
  }
  result
}

match_signature <- function(sets, candidates, universe, output, min_coverage) {
  rows <- list()
  matched <- setNames(vector("list", length(sets)), names(sets))
  for (name in names(sets)) {
    for (id in sets[[name]]) {
      possible <- unique(candidates[[id]])
      expressed <- intersect(possible, universe)
      status <- if (length(expressed) == 1L) "matched" else if (length(expressed) > 1L) {
        "ambiguous_expressed_mapping"
      } else if (length(possible)) "not_in_expression_universe" else "unmapped_annotation"
      rows[[length(rows) + 1L]] <- data.frame(
        gene_set = name, signature_id = id, candidate_ids = paste(possible, collapse = ";"),
        expressed_ids = paste(expressed, collapse = ";"), status = status
      )
      if (status == "matched") matched[[name]] <- c(matched[[name]], expressed)
    }
    matched[[name]] <- unique(matched[[name]])
  }
  audit <- data.table::rbindlist(rows)
  coverage <- data.frame(
    gene_set = names(sets),
    signature_genes = lengths(sets),
    matched_signature_genes = vapply(names(sets), function(name) {
      sum(audit$gene_set == name & audit$status == "matched")
    }, integer(1)),
    unique_expression_features = lengths(matched),
    coverage = lengths(matched) / lengths(sets)
  )
  data.table::fwrite(audit, file.path(output, "signature_mapping_audit.tsv"), sep = "\t")
  data.table::fwrite(coverage, file.path(output, "signature_coverage.tsv"), sep = "\t")
  print(coverage, row.names = FALSE)
  if (any(coverage$coverage < min_coverage)) {
    stop("Insufficient gene-set coverage. Review audits; default 80% is a safety gate, not a validated cutoff.")
  }
  if (any(coverage$coverage < 1)) {
    warning("Incomplete signature coverage: scores use a reduced signature, not the full published signature.")
  }
  if (length(intersect(matched[[1L]], matched[[2L]]))) {
    stop("Positive and negative gene sets overlap after ID mapping.")
  }
  if (any(lengths(matched) >= length(universe))) {
    stop("ssGSEA needs background genes outside each signature. Supply the full transcriptome.")
  }
  list(sets = matched, coverage = coverage)
}

score_ssgsea <- function(x, sets) {
  if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
    parameter <- GSVA::ssgseaParam(
      x, sets, minSize = 1L, maxSize = Inf, alpha = 0.25, normalize = FALSE
    )
    scores <- GSVA::gsva(
      parameter, verbose = TRUE, BPPARAM = BiocParallel::SerialParam()
    )
  } else {
    scores <- GSVA::gsva(
      x, sets, method = "ssgsea", tau = 0.25, ssgsea.norm = FALSE,
      min.sz = 1L, max.sz = Inf, parallel.sz = 1L, verbose = TRUE
    )
  }
  if (!all(SET_NAMES %in% rownames(scores)) ||
      !identical(colnames(scores), colnames(x)) || any(!is.finite(scores))) {
    stop("Unexpected or nonfinite GSVA output.")
  }
  scores[SET_NAMES, , drop = FALSE]
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!length(args) || identical(args, "--help")) {
    usage()
    return(invisible(NULL))
  }
  if (identical(args, "--install-deps")) {
    install_dependencies()
    return(invisible(NULL))
  }
  options <- parse_options(args)
  check_packages()
  if (!file.exists(options$input)) stop("Local expression file not found: ", options$input)
  if (file.exists(options$output)) stop("Output already exists; choose a new directory: ", options$output)
  if (!is.null(options$signature) && !file.exists(options$signature)) {
    stop("Signature file not found: ", options$signature)
  }
  if (!dir.create(options$output, recursive = TRUE)) stop("Cannot create output directory.")
  warning("Exploratory EBPlusPlus scoring: no paper cutoff, scarHRD calibration, or clinical interpretation.")
  signature <- file.path(options$output, "official_expHRD_gene_sets.csv")
  if (is.null(options$signature)) {
    old_timeout <- getOption("timeout")
    options(timeout = max(300, old_timeout))
    on.exit(options(timeout = old_timeout), add = TRUE)
    status <- download.file(SIGNATURE_URL, signature, mode = "wb", method = "libcurl")
    if (status != 0L) stop("Gene-set download failed.")
  } else if (!file.copy(options$signature, signature, overwrite = FALSE)) {
    stop("Could not copy supplied gene-set file.")
  }
  sets <- read_signature(signature)
  expression <- read_expression(
    options$input, options$output, options$gene_id_type, options$sample_limit, options$threads
  )
  x <- expression$x
  candidates <- annotation_candidates(unique(unlist(sets, use.names = FALSE)), expression$type)
  mapping <- match_signature(
    sets, candidates, rownames(x), options$output, options$min_coverage
  )
  message("Scoring ", ncol(x), " samples against ", nrow(x), " background genes.")
  scores <- score_ssgsea(x, mapping$sets)
  results <- data.frame(
    sample_id = colnames(x), patient_id = substr(colnames(x), 1L, 12L),
    positive_ssGSEA = unname(scores[SET_NAMES[1L], ]),
    negative_ssGSEA = unname(scores[SET_NAMES[2L], ]),
    expHRD_exploratory = unname(scores[SET_NAMES[1L], ] - scores[SET_NAMES[2L], ]),
    input_preprocessing = "EBPlusPlus_adjusted_as_supplied",
    interpretation = "exploratory_not_calibrated"
  )
  data.table::fwrite(results, file.path(options$output, "expHRD_EBPlusPlus_scores.tsv"), sep = "\t")
  saveRDS(list(
    options = options, gene_id_type = expression$type,
    input_rows = expression$input_rows, input_samples = expression$input_samples,
    background_genes = rownames(x), scored_samples = colnames(x),
    signature_url = if (is.null(options$signature)) SIGNATURE_URL else NA_character_,
    signature_commit = if (is.null(options$signature)) SIGNATURE_COMMIT else NA_character_,
    input_md5 = tools::md5sum(c(options$input, signature)),
    coverage = mapping$coverage, mapped_gene_sets = mapping$sets,
    ssGSEA = list(alpha = 0.25, normalize = FALSE), session_info = sessionInfo()
  ), file.path(options$output, "run_metadata.rds"))
  writeLines(c(
    "Completed exploratory EBPlusPlus expHRD scoring.",
    "No published cutoff or predicted scarHRD calibration was applied.",
    "Background: all valid, finite, nonconstant genes in this input cohort.",
    "Inspect signature_coverage.tsv and signature_mapping_audit.tsv before interpretation.",
    "Scores may change with gene universe, annotation version, missing genes and GSVA version.",
    "No cancer-type or tumor/normal filter was applied: all requested sample columns were scored."
  ), file.path(options$output, "RUN_COMPLETE.txt"))
  print(head(results))
  message("Completed: ", normalizePath(options$output))
  invisible(results)
}

if (sys.nframe() == 0L) main()

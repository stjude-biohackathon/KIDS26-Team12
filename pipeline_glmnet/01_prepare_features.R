#!/usr/bin/env Rscript
# =============================================================================
# 01_prepare_features.R  —  KIDS26 Team 12
# Build the development feature matrix WITHOUT loading the 27 GB beta.tsv.
#
# * Development partition only; locked_CNS (GBM/LGG) is never read.
# * Streams the matrix in chunks. Label-free prefilter: probes with <= --max_na
#   missing and the --top_k highest variance are kept (default 50,000). The
#   protocol's 5,000-feature selection happens later, INSIDE each training fold.
# * Validation carried over from train_baseline.R: provenance JSON + allowlist
#   hash, beta range [0,1], duplicate patients, HRD component sum, QC annotation.
#
# Usage:
#   Rscript 01_prepare_features.R --beta=/path/beta.tsv --samples=/path/master_samples.tsv \
#       --provenance=/path/beta.provenance.json --out=features
# Subset prototype: add --cancers=BRCA,UCEC,BLCA,STAD,LUSC,HNSC --out=features_subset
# Smoke test: add --max_rows=20000 --top_k=3000 --out=features_smoke
# =============================================================================
suppressPackageStartupMessages(library(data.table))
arg <- function(name, default = NULL) {
  a <- commandArgs(trailingOnly = TRUE)
  hit <- grep(paste0("^--", name, "="), a, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
beta_path    <- arg("beta")
samples_path <- arg("samples")
prov_path    <- arg("provenance", if (!is.null(beta_path)) sub("\\.[^.]+$", ".provenance.json", beta_path))
out_dir      <- arg("out", "features")
top_k        <- as.integer(arg("top_k", "50000"))
max_na       <- as.numeric(arg("max_na", "0.05"))
chunk_lines  <- as.integer(arg("chunk", "4000"))
max_rows     <- as.numeric(arg("max_rows", "Inf"))
locked       <- strsplit(arg("locked", "GBM,LGG"), ",")[[1]]
expected_qc  <- strsplit(arg("expected_qc", "published_450K_no_exclusion"), ",")[[1]]
cancers_arg  <- arg("cancers", NULL)   # e.g. --cancers=BRCA,UCEC,BLCA,STAD,LUSC,HNSC (prototype subset)
stopifnot(!is.null(beta_path), !is.null(samples_path), file.exists(beta_path), file.exists(samples_path))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
t0 <- Sys.time()

# ---- provenance (blocker #4 in train_baseline.R) --------------------------------
if (!file.exists(prov_path)) stop("Missing matrix provenance JSON: ", prov_path)
prov_txt <- paste(readLines(prov_path, warn = FALSE), collapse = " ")
if (grepl('"engineering_only"\\s*:\\s*true', prov_txt)) stop("Engineering-only matrix; not for scientific runs.")
if (!grepl("allowlist", prov_txt)) stop("Provenance JSON has no probe allowlist hash field.")
msg("Provenance OK: ", prov_path)

# ---- samples + metadata checks ---------------------------------------------------
con <- file(beta_path, "r")
hdr <- strsplit(readLines(con, n = 1L), "\t", fixed = TRUE)[[1]]
meta <- fread(samples_path, sep = "\t", na.strings = c("", "NA"))
stopifnot(all(c("sample_id", "patient_id", "cancer_type", "HRDsum", "HRD_LOH", "LST", "TAI") %in% names(meta)))
meta[, sample_id := as.character(sample_id)]

if (any(!hdr[-1] %in% meta$sample_id)) stop("Beta columns without sample metadata: ", sum(!hdr[-1] %in% meta$sample_id))
dev <- meta[!cancer_type %in% locked & sample_id %in% hdr]
if (!is.null(cancers_arg)) {
  want <- strsplit(cancers_arg, ",")[[1]]
  if (any(want %in% locked)) stop("Locked cancer types requested: ", paste(intersect(want, locked), collapse = ","))
  missing_ct <- setdiff(want, dev$cancer_type)
  if (length(missing_ct)) stop("Requested cancer types not in development set: ", paste(missing_ct, collapse = ","))
  dev <- dev[cancer_type %in% want]
  msg("SUBSET MODE: restricted to ", paste(want, collapse = ","))
}
if ("partition" %in% names(dev) && any(dev$partition != "development"))
  stop("Non-development partition among non-CNS samples - check master table.")
if (anyDuplicated(dev$sample_id) || anyDuplicated(dev$patient_id)) stop("Duplicate patient/specimen")
if (any(!is.finite(dev$HRDsum)) || any(dev$HRDsum < 0)) stop("Invalid HRDsum target")
if (!all(abs(dev$HRDsum - dev$HRD_LOH - dev$LST - dev$TAI) < 1e-6)) stop("HRD label sum mismatch")
if ("quality_annotation" %in% names(dev)) {
  qa <- table(dev$quality_annotation, useNA = "ifany")
  msg("quality_annotation values: ", paste(names(qa), qa, sep = "=", collapse = "; "))
  if (anyNA(dev$quality_annotation) || !all(dev$quality_annotation %in% expected_qc))
    stop("Unexpected quality_annotation values (expected: ", paste(expected_qc, collapse = ","),
         "). If these are acceptable, pass --expected_qc=value1,value2")
}
msg("Development samples: ", nrow(dev), " across ", uniqueN(dev$cancer_type),
    " cancer types (locked, not read: ", paste(locked, collapse = "/"), ")")
sel <- c(1L, match(dev$sample_id, hdr))

# ---- stream matrix, keep running top_k by variance ------------------------------
pool_list <- list(); stat_list <- list(); pool_rows <- 0L
n_seen <- 0L; n_drop <- 0L; n_out_of_range <- 0
prune <- function() {
  P <- do.call(rbind, pool_list); S <- rbindlist(stat_list)
  keep <- order(S$variance, decreasing = TRUE)[seq_len(min(top_k, nrow(S)))]
  pool_list <<- list(P[keep, , drop = FALSE]); stat_list <<- list(S[keep]); pool_rows <<- length(keep)
}
repeat {
  n_read <- if (is.finite(max_rows)) min(chunk_lines, max_rows - n_seen) else chunk_lines
  if (n_read <= 0) break
  lines <- readLines(con, n = n_read)
  if (!length(lines)) break
  dt <- fread(text = lines, sep = "\t", header = FALSE, select = sel,
              na.strings = c("NA", "NaN", ""), showProgress = FALSE)
  m <- as.matrix(dt[, -1, with = FALSE]); storage.mode(m) <- "double"; rownames(m) <- dt[[1]]
  n_seen <- n_seen + nrow(m)
  n_out_of_range <- n_out_of_range + sum(m < 0 | m > 1, na.rm = TRUE)
  na_frac <- rowMeans(is.na(m)); mu <- rowMeans(m, na.rm = TRUE)
  v <- rowMeans(m^2, na.rm = TRUE) - mu^2
  ok <- na_frac <= max_na & is.finite(v) & v > 0
  n_drop <- n_drop + sum(!ok)
  pool_list[[length(pool_list) + 1L]] <- m[ok, , drop = FALSE]
  stat_list[[length(stat_list) + 1L]] <- data.table(probe_id = rownames(m)[ok], variance = v[ok], na_frac = na_frac[ok])
  pool_rows <- pool_rows + sum(ok)
  if (pool_rows > 2L * top_k) prune()
  if (n_seen %% (chunk_lines * 10L) == 0L) msg(n_seen, " probes scanned")
}
close(con)
if (n_out_of_range > 0) stop("Found ", n_out_of_range, " beta values outside [0,1]")
if (!length(pool_list)) stop("No probes passed the filters")
prune(); pool <- pool_list[[1]]; pool_stats <- stat_list[[1]]; rm(pool_list, stat_list)
if (anyDuplicated(pool_stats$probe_id)) stop("Duplicate probe IDs")
msg("Scanned ", n_seen, " probes; failed NA/variance filter: ", n_drop, "; kept: ", nrow(pool))

X <- t(pool); rm(pool); invisible(gc())
rownames(X) <- dev$sample_id
stopifnot(identical(rownames(X), dev$sample_id), identical(colnames(X), pool_stats$probe_id))

saveRDS(X, file.path(out_dir, "X_dev.rds"), compress = FALSE)
fwrite(dev, file.path(out_dir, "meta_dev.tsv"), sep = "\t")
fwrite(pool_stats, file.path(out_dir, "probes_prefilter.tsv"), sep = "\t")
invisible(file.copy(prov_path, file.path(out_dir, "beta.provenance.json"), overwrite = TRUE))
writeLines(c(sprintf("n_samples\t%d", nrow(X)), sprintf("n_probes\t%d", ncol(X)),
             sprintf("probes_scanned\t%d", n_seen), sprintf("top_k\t%d", top_k),
             sprintf("max_na\t%g", max_na), sprintf("cancers\t%s", paste(sort(unique(dev$cancer_type)), collapse = ",")), sprintf("beta\t%s", normalizePath(beta_path)),
             sprintf("R\t%s", R.version.string), sprintf("created\t%s", format(Sys.time()))),
           file.path(out_dir, "features_info.tsv"))
msg("Done in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min -> ", out_dir)

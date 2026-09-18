#!/usr/bin/env Rscript
# =============================================================================
# 03_predict_frozen.R  —  memory-safe replacement for scripts/predict_frozen.R
# Streams the beta matrix and keeps ONLY the frozen model's features (<= 5,000
# rows), so a 27 GB file never has to fit in memory.
#
#   Rscript 03_predict_frozen.R --model=results/run01/frozen_nonCNS.rds \
#       --beta=/path/beta.tsv --out=predictions.tsv [--samples=ids.txt]
#
# !! Do NOT run on locked_CNS (GBM/LGG) until development results are reviewed
#    and the model is formally locked. !!
# =============================================================================
suppressPackageStartupMessages(library(data.table))
arg <- function(name, default = NULL) {
  a <- commandArgs(trailingOnly = TRUE)
  hit <- grep(paste0("^--", name, "="), a, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(arg("model_r", file.path(script_dir, "R", "model.R")))
model_path <- arg("model"); beta_path <- arg("beta"); out_path <- arg("out")
stopifnot(!is.null(model_path), !is.null(beta_path), !is.null(out_path))

b <- readRDS(model_path)
feats <- b$preprocess$features
con <- file(beta_path, "r")
hdr <- strsplit(readLines(con, n = 1L), "\t", fixed = TRUE)[[1]]
cols <- seq_along(hdr)[-1]
if (!is.null(arg("samples"))) {
  want <- readLines(arg("samples")); cols <- which(hdr %in% want)
  if (!length(cols)) stop("None of the requested samples are in the beta header")
}
keep_lines <- character(0)
repeat {
  lines <- readLines(con, n = 4000L)
  if (!length(lines)) break
  pid <- sub("\t.*$", "", lines)
  keep_lines <- c(keep_lines, lines[pid %in% feats])
}
close(con)
if (!length(keep_lines)) stop("No model features found in the beta matrix")
dt <- fread(text = keep_lines, sep = "\t", header = FALSE, select = c(1L, cols), na.strings = c("NA", "NaN", ""))
x <- t(as.matrix(dt[, -1, with = FALSE])); storage.mode(x) <- "double"
colnames(x) <- dt[[1]]; rownames(x) <- hdr[cols]
if (anyDuplicated(colnames(x))) stop("Duplicate probes")
if (any(is.finite(x) & (x < 0 | x > 1))) stop("Invalid beta values")
cat("Features found:", ncol(x), "of", length(feats), "\n")

p <- predict_en(b, x); q <- if (is.null(b$interval_q)) Inf else b$interval_q
p$lower <- p$predicted_reference_HRDsum - q; p$upper <- p$predicted_reference_HRDsum + q
p$interval_status <- if (is.finite(q)) "empirical_no_domain_shift_guarantee" else "unavailable_insufficient_calibration"
p$warning <- ifelse(p$qc_fail, "excessive_missing_features",
             ifelse(p$ood, "outside_training_distance_reference", "research_prediction_unvalidated_domain"))
p$estimate_for_display <- ifelse(p$reportable, p$predicted_reference_HRDsum, NA_real_)
p$lower[!p$reportable | !is.finite(q)] <- NA_real_; p$upper[!p$reportable | !is.finite(q)] <- NA_real_
p$model_md5 <- unname(tools::md5sum(model_path))
p$provenance <- "frozen_research_model;domain_validity_requires_review"
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
fwrite(p, out_path, sep = "\t")
cat("Wrote", nrow(p), "predictions to", out_path, "\n")

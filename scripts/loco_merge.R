#!/usr/bin/env Rscript
# =============================================================================
# scripts/loco_merge.R - combine per-fold outputs into cohort-level results.
# =============================================================================
#
# PURPOSE
#   scripts/loco_one_fold.R writes one set of files per outer fold. This script
#   gathers them and produces exactly the files scripts/train_baseline.R would
#   have written from its serial loop, plus the pooled confounding controls.
#
#   Splitting merge from fit matters because the pooled analyses are the ones
#   that actually detect tissue and purity confounding, and they cannot run
#   until every fold is finished.
#
# USAGE
#   Rscript scripts/loco_merge.R <out_dir> <master.tsv>
#
# IT WILL REFUSE TO RUN if any fold is missing. A merged result computed from 27
# of 30 folds would silently be a different (and better-looking) analysis than
# the one specified, because failed folds are rarely failing at random.
# =============================================================================

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("Usage: Rscript scripts/loco_merge.R <out_dir> <master.tsv>")
out_dir <- args[1]; meta_path <- args[2]
fold_dir <- file.path(out_dir, "folds")
if (!dir.exists(fold_dir)) stop("No folds directory at ", fold_dir)

source("R/model.R")

meta <- read.delim(meta_path, check.names=FALSE, stringsAsFactors=FALSE)
# DEFECT 3: the hardcoded list and the `partition` column must agree exactly.
assert_partition_matches_cns(meta$cancer_type, meta$partition,
                             context = basename(meta_path))
cns <- meta$cancer_type %in% c("GBM","LGG")
types <- sort(unique(meta$cancer_type[!cns]))

# --- Completeness check -----------------------------------------------------
# DEFECT 4: this used to glob predictions_*.tsv ALONE. loco_one_fold.R writes
# predictions first and metrics/nullpanel afterwards, so a fold that crashed
# between those writes passed the gate and was merged as if complete - with its
# metrics row simply absent from the macro average. All three file families are
# now required, and their per-fold sets must be identical.
fold_set <- function(prefix) {
  sort(sub(paste0("^", prefix, "_"), "",
           sub("\\.tsv$", "", basename(Sys.glob(file.path(fold_dir, paste0(prefix, "_*.tsv")))))))
}
present  <- fold_set("predictions")
have_mt  <- fold_set("metrics")
have_np  <- fold_set("nullpanel")
missing <- setdiff(types, present)
if (length(missing)) {
  stop("Missing folds: ", paste(missing, collapse=", "),
       "\nRerun those array indices before merging. Refusing to merge a partial cohort.")
}
incomplete <- c(setdiff(present, have_mt), setdiff(present, have_np))
if (length(incomplete)) {
  stop("Folds with predictions but missing metrics_/nullpanel_ output: ",
       paste(sort(unique(incomplete)), collapse=", "),
       "\nThose folds crashed part-way through writing. Rerun them; refusing to merge.")
}
if (!identical(present, have_mt) || !identical(present, have_np)) {
  stop("Per-fold file sets disagree: predictions=", length(present),
       " metrics=", length(have_mt), " nullpanel=", length(have_np),
       "\nRefusing to merge a cohort whose artefacts do not line up fold for fold.")
}
cat("all", length(types), "folds present (predictions + metrics + nullpanel)\n\n")

read_all <- function(pattern) {
  do.call(rbind, lapply(Sys.glob(file.path(fold_dir, pattern)),
                        function(f) read.delim(f, check.names=FALSE, stringsAsFactors=FALSE)))
}
pooled    <- read_all("predictions_*.tsv")
mt        <- read_all("metrics_*.tsv")
all_nulls <- read_all("nullpanel_*.tsv")

# DEFECT 4 (row counts): file presence is necessary but not sufficient. Every
# fold contributes exactly one metrics row and one nullpanel row, and metrics$n
# is the count of that fold's predictions, so the two must reconcile. A fold
# whose metrics file was written from a different (e.g. re-run, differently
# masked) prediction set shows up here and nowhere else.
if (nrow(mt) != length(present) || nrow(all_nulls) != length(present)) {
  stop(sprintf(paste0("Row counts do not match fold counts: %d folds, %d metrics rows, ",
                      "%d nullpanel rows. Refusing to merge."),
               length(present), nrow(mt), nrow(all_nulls)))
}
if (all(c("cancer_type","n") %in% names(mt))) {
  pred_n <- table(pooled$cancer_type)
  mt_n   <- stats::setNames(as.integer(mt$n), as.character(mt$cancer_type))
  common <- intersect(names(pred_n), names(mt_n))
  disagree <- common[as.integer(pred_n[common]) != mt_n[common]]
  if (!identical(sort(names(pred_n)), sort(names(mt_n))) || length(disagree)) {
    stop("metrics_*.tsv sample counts disagree with predictions_*.tsv for: ",
         paste(c(disagree, setdiff(union(names(pred_n),names(mt_n)), common)), collapse=", "),
         "\nA metrics row written from a different prediction set is not mergeable.")
  }
}

# DEFECT 4 (protocol homogeneity): folds fitted on different targets or with
# different probe-ranking strategies are not one experiment, and pooling them
# produces a number that describes no run that was ever performed. Both columns
# are written per fold by loco_one_fold.R precisely so this can be checked.
for (col in c("target_transform","feature_rank")) {
  if (!col %in% names(mt)) {
    stop("metrics_*.tsv has no `", col, "` column, so folds cannot be shown to share ",
         "one protocol. Rerun the folds with the current scripts/loco_one_fold.R.")
  }
  vals <- unique(as.character(mt[[col]]))
  if (length(vals) != 1L) {
    stop("Folds disagree on ", col, ": ", paste(sort(vals), collapse=", "),
         "\nRefusing to merge folds run under different protocols.")
  }
}
cat(sprintf("protocol: target_transform=%s  feature_rank=%s (identical across all %d folds)\n\n",
            mt$target_transform[1], mt$feature_rank[1], nrow(mt)))

write.table(pooled, file.path(out_dir,"loco_predictions.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(mt, file.path(out_dir,"loco_metrics.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(all_nulls, file.path(out_dir,"loco_null_panel_by_cancer.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
writeLines(paste("macro_MAE", mean(mt$MAE)), file.path(out_dir,"macro_metrics.txt"))

# Surface any fold whose tuning hit a path endpoint. One or two is unremarkable;
# a majority means the lambda path needs widening before the result is trusted.
if ("lambda_at_boundary" %in% names(mt)) {
  nb <- sum(as.logical(mt$lambda_at_boundary), na.rm=TRUE)
  cat(sprintf("lambda at path boundary in %d of %d folds\n", nb, nrow(mt)))
  if (nb > nrow(mt)/2) cat("  WARNING: most folds hit a boundary; widen the lambda path.\n")
}

# ===========================================================================
# POOLED CONFOUNDING CONTROLS
# ===========================================================================
# Per-fold, every sample shares one cancer type, so the tissue-mean null is just
# that fold's mean. Pooling restores between-tissue variation, which is the only
# setting where the tissue confound is visible. See docs/22.
pooled_null <- null_panel(pooled$actual, pooled$predicted_reference_HRDsum,
                          pooled$cancer_type, mean(meta$HRDsum[!cns]))
write.table(pooled_null, file.path(out_dir,"pooled_null_panel.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

pooled_within <- within_tissue_metrics(pooled$actual, pooled$predicted_reference_HRDsum, pooled$cancer_type)
write.table(pooled_within, file.path(out_dir,"pooled_within_tissue_metrics.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

pooled_perm <- permutation_test_within_tissue(pooled$actual, pooled$predicted_reference_HRDsum,
                                              pooled$cancer_type, n_perm=1000L)
write.table(pooled_perm, file.path(out_dir,"pooled_within_tissue_permutation.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# --- Lineage confounding: tissue-identity R^2 (DEFECT 4) --------------------
# R^2 of a quantity on factor(cancer_type) - "how much of this is just knowing
# which tissue it is". This is the 56.2% headline from run 01, and it was only
# ever computed ad hoc inside scripts/c1_zeroshot_percentile.R, so a standard
# LOCO run did not ship it. It does now, on every run.
#
# BOTH rows matter and neither should be quoted alone. The PREDICTION's R^2 is
# the confounding number; the OBSERVED label's R^2 is the reference showing how
# much tissue structure HRDsum genuinely has. A prediction R^2 far above the
# label's means the model is more lineage-determined than the biology is.
tissue_r2 <- data.frame(
  quantity = c("predicted_reference_HRDsum", "actual_HRDsum"),
  r2_tissue_identity = c(tissue_identity_r2(pooled$predicted_reference_HRDsum, pooled$cancer_type),
                         tissue_identity_r2(pooled$actual, pooled$cancer_type)),
  n = c(sum(is.finite(pooled$predicted_reference_HRDsum)), sum(is.finite(pooled$actual))),
  n_tissues = length(unique(pooled$cancer_type)),
  note = c("fraction of the MODEL OUTPUT explained by tissue identity alone (lineage confounding)",
           "fraction of the OBSERVED label explained by tissue identity alone (the reference point)"),
  stringsAsFactors = FALSE)
write.table(tissue_r2, file.path(out_dir,"pooled_tissue_identity_r2.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
cat(sprintf("TISSUE  R2(predicted ~ factor(cancer_type))=%.4f  R2(observed ~ factor(cancer_type))=%.4f\n",
            tissue_r2$r2_tissue_identity[1], tissue_r2$r2_tissue_identity[2]))

cat(sprintf("\nPOOLED  MAE_model=%.3f  MAE_tissue_mean_null=%.3f  skill_vs_tissue_mean=%.3f\n",
            pooled_null$MAE_model, pooled_null$MAE_tissue_mean_null, pooled_null$skill_vs_tissue_mean))
if (is.finite(pooled_null$skill_vs_tissue_mean) && pooled_null$skill_vs_tissue_mean <= 0) {
  cat("WARNING: the model does NOT beat the tissue-mean null. Consistent with the model\n",
      "having learned tissue lineage rather than HRD biology; do not present this as an HRD result.\n")
}

# --- Purity confounding -----------------------------------------------------
if ("purity" %in% names(meta)) {
  pur <- meta$purity[match(pooled$patient_id, meta$patient_id)]
  if (sum(is.finite(pur)) >= 30L) {
    legs <- purity_confound_legs(pooled$actual, pooled$predicted_reference_HRDsum, pur, pooled$cancer_type)
    write.table(legs, file.path(out_dir,"pooled_purity_confound_legs.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
    ps <- purity_stratified_metrics(pooled$actual, pooled$predicted_reference_HRDsum, pur,
                                    pooled$cancer_type, mean(meta$HRDsum[!cns]))
    if (!is.null(ps)) write.table(ps, file.path(out_dir,"pooled_purity_stratified.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
    pm <- purity_matched_subset(pooled$actual, pooled$predicted_reference_HRDsum, pur,
                                pooled$cancer_type, training_mean=mean(meta$HRDsum[!cns]))
    if (!is.null(pm)) write.table(pm, file.path(out_dir,"pooled_purity_matched.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
    cat(sprintf("PURITY  cor(pred,purity)_within=%.3f  cor(label,purity)_within=%.3f\n",
                legs$cor_pred_purity_within, legs$cor_label_purity_within))
  } else {
    cat("Too few finite purity values for stratified analysis; skipping.\n")
  }
} else {
  cat("No purity column in master table; purity confounding NOT assessed.\n")
}

cat("\nmerge complete ->", out_dir, "\n")

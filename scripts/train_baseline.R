# =============================================================================
# scripts/train_baseline.R - Development training and leave-one-cancer-out
#                            evaluation for the methylation -> HRDsum model.
# =============================================================================
#
# USAGE
#   Rscript scripts/train_baseline.R <beta.tsv> <master_samples.tsv> <output_dir>
#
#   e.g. Rscript scripts/train_baseline.R \
#          data/processed/beta.tsv \
#          data/processed/master_samples.tsv \
#          results/baseline
#
#   Must be run from the repository root: source("R/model.R") below is a
#   relative path with no fallback resolution.
#
# WHAT THIS SCRIPT DOES, IN ORDER
#   1. Load and validate the beta matrix and the master sample table.
#   2. Refuse to run on an engineering/smoke fixture (provenance gate).
#   3. Hold out GBM and LGG entirely - the "locked CNS" partition.
#   4. Leave-one-cancer-out over the remaining development cancers, writing
#      per-fold predictions, metrics and fitted bundles.
#   5. Reserve a stratified 20% calibration set, fit a final model on the rest,
#      derive a split-conformal interval, and freeze that model to disk.
#
# WHAT IT DELIBERATELY DOES NOT DO
#   It never evaluates on the locked CNS partition. Architecture and
#   hyperparameter decisions must be finalised BEFORE anything touches the
#   locked set, otherwise the locked set stops being a held-out test and
#   becomes just another tuning signal. See the note at the bottom.
#
# ---------------------------------------------------------------------------
# OPERATIONAL WARNING - THIS SCRIPT WILL NOT RUN AS-IS AT FULL SCALE
# ---------------------------------------------------------------------------
# Measured against the real data/processed/beta.tsv (336,480 probes x 7,707
# samples):
#
#   MEMORY: line 6's fread() materialises ~21 GB as a data.frame, and the
#   transpose chain on the x <- t(as.matrix(...)) line creates three to four
#   live copies before modelling starts - 60 GB or more at peak. `beta` is
#   never rm()'d. Submit to a large-memory LSF queue; do not run on a login
#   node (typically ~3 GB free).
#
#   TIME: fit_preprocess() in R/model.R uses apply() over ~336k columns,
#   benchmarked at ~2.8 min per call. Nested LOCO invokes it roughly 900 times
#   (30 inner folds x 30 outer folds) -> about 42 hours. Substituting
#   matrixStats::colMedians/colVars brings that to ~17 hours.
#
# Recommended fixes before a real run: convert the matrix to a float HDF5 or
# bigstatsr FBM backing (halves it to ~10 GB), swap apply() for matrixStats,
# and add rm(beta); gc() after the transpose.
#
# ---------------------------------------------------------------------------
# AUDIT NOTES ON INTERPRETING THE OUTPUT (2026-09-16)
# ---------------------------------------------------------------------------
#   * The ONLY null comparator produced is null_MAE, the training-mean
#     prediction. There is no permuted-label control and no cancer-type-mean
#     (tissue identity) control. HRDsum varies strongly by tissue, so without a
#     tissue null you cannot tell whether the model learned HRD biology or
#     learned to recognise the tumour type. Add these before presenting results.
#   * null_MAE uses the TRAINING mean as its baseline, while the R2 reported by
#     metrics() uses the HELD-OUT CANCER's own mean. They are different
#     baselines; do not compare them side by side without saying so.
#   * `selected_features` counts features ENTERING glmnet, not features with
#     non-zero coefficients. It is not a sparsity measure.
#   * The LOCO metrics and the frozen model come from DIFFERENT fits - the
#     frozen model excludes the 20% calibration reservation. macro_metrics.txt
#     does not describe the artefact that ships in frozen_nonCNS.rds.
#   * purity and ploidy are present in the master table and never used here,
#     despite being plausible confounders (methylation predicts tumour purity;
#     purity correlates with scar scores).
# =============================================================================

# Run from repository root: Rscript scripts/train_baseline.R beta.tsv master.tsv results/run01
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3) stop("Usage: Rscript scripts/train_baseline.R beta.tsv master.tsv output_dir")
# All modelling functions live in R/model.R. Relative path => run from repo root.
source("R/model.R")
if (!requireNamespace("data.table",quietly=TRUE)) stop("Install dependencies with scripts/setup.R")

# --- Load the beta matrix --------------------------------------------------
# On-disk layout is probes in ROWS, samples in COLUMNS, first column = probe_id.
# check.names=FALSE preserves TCGA barcodes verbatim; R would otherwise mangle
# the hyphens in "TCGA-A2-A0SV-01A-..." into dots.
# MEMORY: this is the ~21 GB allocation described in the header.
#
# B4 NOTE - DO NOT re-add a `select=` probe allowlist here. Commit 4ec3493 passed
#   select = c("probe_id", readLines("config/shared_autosomal_probes.txt"))
# which is wrong on two independent counts:
#
#   1. ORIENTATION. `select=` chooses COLUMNS by name, but this matrix stores
#      probes in ROWS and sample barcodes in COLUMNS (header is
#      "probe_id  TCGA-OR-A5J1-01A-...  TCGA-OR-A5J2-01A-..."). Every probe ID is
#      therefore "not found in column name header" and is SKIPPED WITH A WARNING,
#      not an error. fread returns a single probe_id column, and the next line
#      (t(as.matrix(beta[,-1]))) then transposes a zero-column frame. Silent
#      wrong-shape data, not a crash.
#
#   2. REDUNDANCY. prepare_beta.py already intersected the matrix with the
#      384,640-probe bridge at build time. data/processed/beta.provenance.json
#      records probe_allowlist_sha256=f359e43b... and n_probes=336480, and the
#      file on disk has 336,480 rows. There is nothing left to filter out.
#
# The genuine memory saving from that commit was rm(beta)+gc() below, which is
# retained. See docs/21_BLOCKER_RESOLUTION_PLAN.md (B4).
beta <- data.table::fread(args[1],data.table=FALSE,check.names=FALSE)
if (anyDuplicated(beta[[1]])) stop("Duplicate probe IDs")

# Transpose to the samples-in-rows orientation that R/model.R expects, then
# label columns with probe IDs so apply_preprocess() can match features BY NAME.
# storage.mode ensures a numeric matrix even if fread typed a column as integer
# or character. NOTE: if any cell were non-numeric this coercion emits a warning
# ("NAs introduced by coercion") rather than an error, and those NAs would then
# flow onward looking like legitimate missing data.
x <- t(as.matrix(beta[,-1,drop=FALSE])); storage.mode(x)<-"double";colnames(x)<-beta[[1]]
rm(beta); gc()

# Beta values are methylation proportions and must lie in [0,1]. The is.finite()
# conjunct means NA is tolerated (handled later by imputation) but an out-of-range
# real number is fatal. Caveat: Inf passes this check, because is.finite(Inf) is
# FALSE - see the related note in R/model.R's fit_preprocess().
if(any(is.finite(x)&(x<0|x>1))) stop("Expected beta values in [0,1]")

# --- Load and validate the sample metadata ---------------------------------
meta <- read.delim(args[2],check.names=FALSE,stringsAsFactors=FALSE)
# One row per patient AND per specimen. Duplicates would leak a patient across
# the fold boundary during cross-validation.
if(anyDuplicated(meta$sample_id)||anyDuplicated(meta$patient_id)) stop("Duplicate patient/specimen")
# Every matrix column must have metadata. The reverse is allowed: extra metadata
# rows (samples not in this matrix) are simply dropped by the match() below.
if(any(!rownames(x)%in%meta$sample_id)) stop("Missing sample metadata")
# Align metadata row order to matrix row order. EVERY subsequent meta$... lookup
# depends on this alignment holding, so it must come before any modelling.
meta <- meta[match(rownames(x),meta$sample_id),]

# Target sanity: HRDsum is a count of genomic scars, so it cannot be negative.
if(any(!is.finite(meta$HRDsum))||any(meta$HRDsum<0)) stop("Invalid target")
# Label integrity: HRDsum must equal the sum of its three components. This is
# the single most valuable check in the file - it catches column mis-mapping
# between HRD_LOH / LST / TAI, which would otherwise be invisible and would
# silently train the model against a scrambled target.
if(!all(abs(meta$HRDsum-meta$HRD_LOH-meta$LST-meta$TAI)<1e-6)) stop("Label sum mismatch")

# --- QC gate ----------------------------------------------------------------
# The previous version of this gate compared quality_annotation against the
# constant string build_master.py writes for every row, so it could never fail.
# A gate that cannot fail is worse than no gate: it puts a reassuring line in
# the audit trail while verifying nothing.
#
# The replacement asserts properties that can actually be violated. Each check
# is separate so a failure names the specific problem.

# 1. The column must exist and be fully populated.
if(!"quality_annotation" %in% names(meta)) stop("Missing quality_annotation column")
if(any(is.na(meta$quality_annotation)|!nzchar(trimws(as.character(meta$quality_annotation)))))
  stop("Blank quality_annotation for ",
       sum(is.na(meta$quality_annotation)|!nzchar(trimws(as.character(meta$quality_annotation)))),
       " specimens")

# 2. Every annotation must be one the pipeline knows how to interpret. An
#    unrecognised value means build_master.py changed and this gate was not
#    revisited - stop rather than guess what it meant.
known_qc <- c("published_450K_no_exclusion","qc_pass","qc_pass_reviewed")
if(!all(meta$quality_annotation %in% known_qc))
  stop("Unrecognised quality_annotation value(s): ",
       paste(sort(unique(meta$quality_annotation[!meta$quality_annotation %in% known_qc])),collapse=", "))

# 3. The real risk this gate guards against: silently training on a cohort
#    whose QC was never adjudicated per specimen. If every row carries the
#    identical blanket annotation, no per-sample QC decision was made, so warn
#    loudly and record it rather than presenting the run as QC-gated.
#    Set KIDS26_REQUIRE_ADJUDICATED_QC=1 to make this fatal.
qc_is_blanket <- length(unique(meta$quality_annotation))==1L
if(qc_is_blanket) {
  msg <- paste0("Blanket quality_annotation ('",meta$quality_annotation[1],
                "') on all ",nrow(meta)," specimens: no per-sample QC adjudication is recorded. ",
                "Array-level QC (detection p-value pass rate, sex concordance, duplicate audit) ",
                "has NOT gated this cohort. See docs/16_COHORT_QC.md.")
  if(nzchar(Sys.getenv("KIDS26_REQUIRE_ADJUDICATED_QC"))) stop(msg)
  warning(msg,call.=FALSE);message("WARNING: ",msg)
}

# 4. Partition bookkeeping must be internally consistent where present.
if("partition" %in% names(meta) && any(is.na(meta$partition)))
  stop("Missing partition assignment for ",sum(is.na(meta$partition))," specimens")

# --- Provenance gate: refuse to train on a smoke fixture --------------------
# scripts/prepare_beta.py writes a sidecar <matrix_basename>.provenance.json
# recording whether the matrix was built with a real probe allowlist or with the
# tiny --smoke-probes engineering subset. Without this gate it would be
# alarmingly easy to present fixture output as a scientific result.
# Never silently run a smoke-probe fixture as scientific evidence.
prov <- sub("\\.[^.]+$",".provenance.json",args[1])
if(!file.exists(prov)) stop("Missing matrix provenance JSON")
if(!requireNamespace("jsonlite",quietly=TRUE)) stop("Install jsonlite")
pr <- jsonlite::fromJSON(prov)
if(isTRUE(pr$engineering_only)) stop("Engineering-only matrix. Build a real technical probe allowlist first.")
if(is.null(pr$probe_allowlist_sha256)) stop("Missing technical allowlist hash")

out <- args[3];dir.create(out,recursive=TRUE,showWarnings=FALSE)

# --- Define the locked CNS partition ---------------------------------------
# GBM and LGG (642 samples) are withheld as the adult-CNS stand-in for the
# eventual pediatric high-grade glioma transfer question. They are never fitted
# and never scored in this script.
#
# NOTE: this re-derives the lock from a hardcoded cancer-type list rather than
# reading the `partition` column that build_master.py already wrote. The two
# agree today (verified: 642 locked_CNS rows, exactly GBM+LGG, zero anomalies),
# but they are two sources of truth that could drift. Reading partition and
# stopifnot()-ing agreement would be safer.
cns <- meta$cancer_type %in% c("GBM","LGG")

# Anti-fixture guards: refuse to make transfer claims from a toy cohort.
if(length(unique(meta$cancer_type[!cns]))<3) stop("Need >=3 development cancers; six-patient smoke subset is insufficient")
if(sum(!cns)<50) stop("Need >=50 development patients for this protocol; do not claim transfer from a tiny fixture")

# ===========================================================================
# PHASE 1: Leave-one-cancer-out (LOCO) development evaluation
# ===========================================================================
# Each iteration holds out one entire cancer type. This tests whether the model
# has learned something transferable about HRD biology, or has merely learned
# tissue-specific methylation signatures. Holding out random SAMPLES instead
# would let the model recognise the tissue and look far better than it is.
all_predictions <- list();all_metrics<-list();all_nulls<-list()
for (type in sort(unique(meta$cancer_type[!cns]))) {
 # Both masks carry !cns, so locked CNS samples appear in neither train nor test.
 tr <- !cns & meta$cancer_type!=type;te <- !cns & meta$cancer_type==type

 # fit_en() runs its own inner CV for hyperparameters using only these training
 # rows - the held-out cancer is invisible to hyperparameter selection.
 b <- fit_en(x[tr,,drop=FALSE],meta$HRDsum[tr],meta$patient_id[tr],meta$cancer_type[tr])

 p <- predict_en(b,x[te,,drop=FALSE]);p$actual<-meta$HRDsum[te];p$cancer_type<-type;p$patient_id<-meta$patient_id[te]
 # The null comparator: predict this fold's training mean for every sample.
 # Stored per-row so downstream analysis can compute any null metric it wants.
 p$null_prediction<-b$training_mean;p$split<-"development_LOCO"
 all_predictions[[type]]<-p

 # Metrics are computed over ALL held-out samples, including those flagged
 # non-reportable. That is the non-inflating choice: reporting only "confident"
 # predictions would be selective reporting. abstention_rate is tracked
 # separately so the trade-off stays visible.
 mm<-metrics(p$actual,p$predicted_reference_HRDsum);mm$cancer_type<-type;mm$train_n<-sum(tr);mm$selected_features<-length(b$preprocess$features);mm$alpha<-b$alpha;mm$lambda<-b$lambda;mm$null_MAE<-mean(abs(p$actual-b$training_mean));mm$abstention_rate<-mean(!p$reportable)
 all_metrics[[type]]<-mm

 # Within-fold null panel. Every sample here shares one cancer type, so the
 # tissue-mean null collapses to this fold's own mean - a strict bar, but it
 # says nothing about the between-tissue confound. The pooled analysis after
 # the loop is where that becomes visible.
 np<-null_panel(p$actual,p$predicted_reference_HRDsum,p$cancer_type,b$training_mean)
 np$cancer_type<-type;all_nulls[[type]]<-np

 # Persist the per-fold bundle and inner-fold assignments so the grid search and
 # fold structure remain auditable after the run.
 saveRDS(b,file.path(out,paste0("loco_",type,".rds")))
 write.table(b$inner_folds,file.path(out,paste0("inner_folds_",type,".tsv")),sep="\t",row.names=FALSE,quote=FALSE)
}

write.table(do.call(rbind,all_predictions),file.path(out,"loco_predictions.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
mt<-do.call(rbind,all_metrics);write.table(mt,file.path(out,"loco_metrics.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
# Macro average: every cancer counts equally regardless of sample count, so BRCA
# (n=743) cannot drown out CHOL (n=35). Propagates NA loudly if any fold failed.
writeLines(paste("macro_MAE",mean(mt$MAE)),file.path(out,"macro_metrics.txt"))

# ===========================================================================
# POOLED CONFOUNDING CONTROLS
# ===========================================================================
# Everything above is per-fold, where each fold is one cancer type and the
# tissue-mean null therefore collapses to that fold's own mean. Pooling all
# held-out predictions restores the between-tissue variation, which is the only
# way the tissue confound becomes visible.
#
# The decisive question: does the model beat an opponent that knows ONLY each
# sample's cancer type and that type's average HRDsum? If not, the model has
# learned tissue lineage rather than HRD biology, and no other result in this
# run should be presented as evidence of the latter. See docs/22.
pooled<-do.call(rbind,all_predictions)
write.table(do.call(rbind,all_nulls),file.path(out,"loco_null_panel_by_cancer.tsv"),sep="\t",row.names=FALSE,quote=FALSE)

pooled_null <- null_panel(pooled$actual,pooled$predicted_reference_HRDsum,
                          pooled$cancer_type,mean(meta$HRDsum[!cns]))
write.table(pooled_null,file.path(out,"pooled_null_panel.tsv"),sep="\t",row.names=FALSE,quote=FALSE)

# E2: the within-tissue estimand, with both outcome and prediction centered by
# cancer type. A model whose entire skill is tissue-level offsets scores ~0.
pooled_within <- within_tissue_metrics(pooled$actual,pooled$predicted_reference_HRDsum,pooled$cancer_type)
write.table(pooled_within,file.path(out,"pooled_within_tissue_metrics.tsv"),sep="\t",row.names=FALSE,quote=FALSE)

# Permutation null for E2. Predictions are held fixed and the outcome is
# shuffled within tissue, so this asks whether the within-tissue ordering
# carries information - NOT whether the pipeline could manufacture it from
# noise, which would require refitting per permutation.
pooled_perm <- permutation_test_within_tissue(pooled$actual,pooled$predicted_reference_HRDsum,
                                              pooled$cancer_type,n_perm=1000L)
write.table(pooled_perm,file.path(out,"pooled_within_tissue_permutation.tsv"),sep="\t",row.names=FALSE,quote=FALSE)

# Surface the decisive comparison in the run log so it cannot be overlooked.
message(sprintf("POOLED  MAE_model=%.3f  MAE_tissue_mean_null=%.3f  skill_vs_tissue_mean=%.3f",
                pooled_null$MAE_model,pooled_null$MAE_tissue_mean_null,pooled_null$skill_vs_tissue_mean))
if(is.finite(pooled_null$skill_vs_tissue_mean)&&pooled_null$skill_vs_tissue_mean<=0)
  message("WARNING: the model does NOT beat the tissue-mean null. Consistent with the model ",
          "having learned tissue lineage rather than HRD biology; do not present this as an HRD result.")

# --- Purity confounding -----------------------------------------------------
# Purity is excluded from the features but is NOT thereby adjusted for. It is
# readable from methylation AND shapes the ABSOLUTE-derived label, so a model
# can score well by inferring purity and exploiting its correlation with the
# outcome. These analyses measure whether that is happening. See docs/22.
if("purity" %in% names(meta)) {
  pur <- meta$purity[match(pooled$patient_id,meta$patient_id)]
  n_pur <- sum(is.finite(pur))
  message(sprintf("purity available for %d of %d pooled predictions",n_pur,nrow(pooled)))

  if(n_pur>=30L) {
    legs <- purity_confound_legs(pooled$actual,pooled$predicted_reference_HRDsum,pur,pooled$cancer_type)
    write.table(legs,file.path(out,"pooled_purity_confound_legs.tsv"),sep="\t",row.names=FALSE,quote=FALSE)

    ps <- purity_stratified_metrics(pooled$actual,pooled$predicted_reference_HRDsum,pur,
                                    pooled$cancer_type,mean(meta$HRDsum[!cns]))
    if(!is.null(ps)) write.table(ps,file.path(out,"pooled_purity_stratified.tsv"),sep="\t",row.names=FALSE,quote=FALSE)

    pm <- purity_matched_subset(pooled$actual,pooled$predicted_reference_HRDsum,pur,
                                pooled$cancer_type,training_mean=mean(meta$HRDsum[!cns]))
    if(!is.null(pm)) write.table(pm,file.path(out,"pooled_purity_matched.tsv"),sep="\t",row.names=FALSE,quote=FALSE)

    # Both legs must hold for the confound to operate; report them together so
    # the pathway is judged as a whole rather than from one suggestive number.
    message(sprintf("PURITY  cor(pred,purity)_within=%.3f  cor(label,purity)_within=%.3f",
                    legs$cor_pred_purity_within,legs$cor_label_purity_within))
    if(!is.null(ps)&&nrow(ps)>1L) {
      rng <- range(ps$skill_vs_tissue_mean[is.finite(ps$skill_vs_tissue_mean)])
      if(length(rng)==2L) message(sprintf("PURITY  skill_vs_tissue_mean across strata: %.3f to %.3f",rng[1],rng[2]))
    }
  } else {
    message("Too few finite purity values for stratified analysis; skipping.")
  }
} else {
  message("No purity column in master table; purity confounding NOT assessed.")
}

# ===========================================================================
# PHASE 2: Fit and freeze the final model, with a conformal interval
# ===========================================================================
# Deterministic calibration reservation within non-CNS cancers; no refit on calibration.
#
# Stratified 20% per cancer type, so the calibration set mirrors the development
# cohort's tissue composition rather than being dominated by the largest cancers.
# The seed makes the split reproducible.
#
# FIXED: sample(ii, k) where ii has length 1 is interpreted by R as
# sample(1:ii, k) - it returns a random integer in 1..ii rather than ii itself,
# which would point at an arbitrary row (possibly a locked_CNS one). It did not
# fire on the current cohort because the smallest development group is OV with
# n=10, but it would fire on any subset or re-partition producing a singleton
# group. sample.int(length(ii)) is length-safe: for length 1 it yields index 1,
# so ii[1] is returned as intended.
pick_calibration <- function(ii) {
  k <- max(1L,floor(length(ii)*0.2))
  ii[sample.int(length(ii),k)]
}
set.seed(260910);dev<-which(!cns)
cal<-unlist(lapply(split(dev,meta$cancer_type[dev]),pick_calibration),use.names=FALSE)

# Assert the reservation is sane before anything depends on it: calibration
# rows must be development rows, distinct, and must not exhaust any group.
stopifnot(!anyDuplicated(cal),all(cal %in% dev),!any(cns[cal]))
tr<-setdiff(dev,cal)
if(!length(tr)) stop("Calibration reservation consumed every development sample")

# Fit on development-minus-calibration only. The calibration samples must stay
# unseen or the conformal interval below loses its validity guarantee.
b<-fit_en(x[tr,,drop=FALSE],meta$HRDsum[tr],meta$patient_id[tr],meta$cancer_type[tr])

# Split-conformal: predict the untouched calibration set, then take the residual
# quantile as the interval half-width. No refit happens on `cal`.
cp<-predict_en(b,x[cal,,drop=FALSE]);b$interval_q<-conformal_q(meta$HRDsum[cal],cp$predicted_reference_HRDsum)

# Record exactly which patients were used for what, so any later question about
# contamination can be answered from the frozen object itself.
b$calibration_ids<-meta$patient_id[cal];b$training_ids<-meta$patient_id[tr];b$probe_provenance<-pr
b$interval_note<-"Split-conformal residual interval; no pediatric/domain-shift coverage guarantee"

# THE SHIPPING ARTEFACT. scripts/predict_frozen.R consumes this file.
saveRDS(b,file.path(out,"frozen_nonCNS.rds"))

# Export coefficients for interpretation. Rows with coefficient 0 were shrunk
# out by the elastic net; the non-zero count here is the real sparsity measure
# (unlike the `selected_features` column in loco_metrics.tsv).
cc<-as.matrix(coef(b$model,s=b$lambda))
write.table(data.frame(feature=rownames(cc),coefficient=as.numeric(cc[,1]),stringsAsFactors=FALSE),file.path(out,"frozen_coefficients.tsv"),sep="\t",row.names=FALSE,quote=FALSE)

# Full role assignment for every sample: locked_CNS / calibration / training.
# This is the audit trail proving the locked set was never fitted.
write.table(data.frame(sample_id=meta$sample_id,patient_id=meta$patient_id,cancer_type=meta$cancer_type,role=ifelse(cns,"locked_CNS",ifelse(seq_len(nrow(meta))%in%cal,"calibration","training"))),file.path(out,"final_partitions.tsv"),sep="\t",row.names=FALSE,quote=FALSE)

# Deliberately do not evaluate CNS here: architecture selection must precede locked testing.
# Scoring the locked set now would convert it from a held-out test into part of
# the development loop. Evaluate it exactly once, after the model is frozen, via
# scripts/predict_frozen.R.
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
cat("Development LOCO complete; frozen_nonCNS.rds saved. Review development and lock before external prediction.\n")

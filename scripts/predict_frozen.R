# =============================================================================
# scripts/predict_frozen.R - Apply a frozen model to new samples.
# =============================================================================
#
# USAGE
#   Rscript scripts/predict_frozen.R <model.rds> <beta.tsv> <output.tsv> \
#     [--allow-fixture]
#
#   e.g. Rscript scripts/predict_frozen.R \
#          results/baseline/frozen_nonCNS.rds \
#          data/processed/locked_test_beta.tsv \
#          results/locked_predictions.tsv
#
#   --allow-fixture  score a matrix that fails the provenance gate (engineering
#                    fixture, missing sidecar, or unverifiable probe allowlist).
#                    For engineering only. The output is stamped as
#                    non-scientific and the stamp cannot be suppressed.
#
#   Run from the repository root (source("R/model.R") is a relative path).
#
# PURPOSE AND DISCIPLINE
# ----------------------
# This is the inference half of the leakage boundary. It LOADS a model that was
# frozen by scripts/train_baseline.R and applies it unchanged. It performs no
# fitting, no tuning, and no re-estimation of any statistic - all preprocessing
# constants (medians, centres, scales), the coefficients, the OOD cutoff and the
# conformal half-width travel inside the .rds bundle.
#
# This script is what you run ONCE against the locked CNS partition after the
# model is frozen, and what you would run against approved external data. Every
# additional run against the same held-out set erodes its status as a held-out
# set, because the results inevitably inform the next decision.
#
# OUTPUT COLUMN GUIDE
# -------------------
#   predicted_reference_HRDsum  raw model output, ALWAYS populated even when QC
#                               fails. Analysts should use this one.
#   estimate_for_display        the same value, blanked to NA unless reportable.
#                               Display layers (the Shiny app) should use this.
#   lower / upper               conformal interval, NA unless reportable and a
#                               finite calibration quantile exists.
#   reportable                  passed QC AND inside the training distribution.
#   warning                     why a prediction is being withheld, if it is.
#   interval_status             whether the interval is usable at all.
#   model_md5                   checksum of the model file, for provenance.
#   matrix_provenance_class     what the scored matrix was verified to be.
#   probe_set_match             whether its probe allowlist matches the model's.
#
# KNOWN GAPS (audited 2026-09-16)
# --------------------------------
#   * PROVENANCE GATE - RESOLVED 2026-09-16 (blocker B7). The sidecar check that
#     train_baseline.R performs now runs here too, via R/provenance.R, before the
#     matrix is read. Fixtures, unlabelled matrices and unverifiable allowlists
#     are refused unless --allow-fixture is passed, and the resulting class is
#     written into every output row so that app/app.R - which validates
#     KIDS26_DEMO_RESULTS structurally only - displays what it actually got.
#   * There is no HRD_high_probability column. It previously existed and was
#     written as NA for every row; it was removed rather than left as a
#     placeholder, because an always-NA column implies a validated threshold
#     exists. The exploratory cut of 42 in config/analysis_protocol.json
#     describes the TCGA reference distribution and is not a pediatric clinical
#     threshold. A calibrated probability needs a separately fitted and
#     validated mapping.
#     Either wire it up or drop the column - an always-NA field in a
#     clinical-looking table invites misreading.
# =============================================================================

# Rscript scripts/predict_frozen.R model.rds beta.tsv output.tsv [--allow-fixture]
args<-commandArgs(trailingOnly=TRUE)
allow_fixture<-"--allow-fixture" %in% args;args<-args[args!="--allow-fixture"]
if(length(args)!=3)stop("Expected model, beta matrix, output [--allow-fixture]")
source("R/model.R");source("R/provenance.R")
if(!requireNamespace("data.table",quietly=TRUE))stop("Install data.table")

# --- Provenance gate (B7) --------------------------------------------------
# Runs BEFORE the matrix is read, so a fixture costs nothing and an unlabelled
# 28 GB matrix is not loaded only to be rejected afterwards. train_baseline.R
# applies the same gate at fit time; this is the inference half.
matrix_pr<-read_matrix_provenance(args[2])
gate<-gate_matrix_provenance(matrix_pr,allow_fixture=allow_fixture)

# Load the frozen bundle and the new beta matrix. Same on-disk orientation as
# training: probes in rows, samples in columns, first column = probe_id.
b<-readRDS(args[1]);d<-data.table::fread(args[2],data.table=FALSE,check.names=FALSE)
if(anyDuplicated(d[[1]]))stop("Duplicate probes")

# NOTE (merge 2026-09-17): a second, independently written provenance gate was
# added here in commit a6de701. It has been removed rather than kept alongside
# the R/provenance.R gate above, for three reasons:
#   1. It was UNREACHABLE. gate_matrix_provenance() runs before the matrix is
#      read and stops on the same conditions, so this block never executed.
#   2. It ran AFTER fread(), so on an unlabelled 28 GB matrix it would have
#      loaded the whole file before refusing it.
#   3. It checked only engineering_only, missing the no-allowlist-hash case,
#      and it did not stamp the output - so a fixture scored under
#      --allow-fixture would still have produced a clean-looking results table.
# The two path conventions were verified to agree before removal:
#   sub("\\.[^.]+$", ...) and tools::file_path_sans_ext() both map
#   data/processed/beta.tsv -> data/processed/beta.provenance.json
# Did this matrix come off the same probe allowlist the model was trained on?
# Non-fatal by design (see R/provenance.R): cross-platform transfer is the point
# of the project, and features align by name. But it is recorded per row.
probe_match<-compare_probe_allowlist(b$probe_provenance,matrix_pr)
if(probe_match!="matched")
  warning(sprintf("Probe allowlist vs training matrix: %s. Absent probes are imputed with training medians; check missing_fraction.",probe_match),call.=FALSE)

# Transpose to samples-in-rows and name the columns by probe ID. The naming
# matters: apply_preprocess() inside predict_en() aligns features BY NAME, so
# the incoming matrix does not need the same probe order, or even the same probe
# set, as the training matrix. Missing probes are imputed with stored training
# medians and counted toward the per-sample missing fraction.
x<-t(as.matrix(d[,-1,drop=FALSE]));storage.mode(x)<-"double";colnames(x)<-d[[1]]
if(any(is.finite(x)&(x<0|x>1)))stop("Invalid beta")

# Apply the frozen model. predict_en() returns the raw prediction plus the
# missing_fraction / ood / qc_fail / reportable flags.
p<-predict_en(b,x);q<-b$interval_q
# A model frozen without a calibration step has no interval half-width; Inf
# propagates to NA bounds below rather than fabricating a narrow interval.
if(is.null(q))q<-Inf

# Symmetric conformal interval around the point estimate.
p$lower<-p$predicted_reference_HRDsum-q;p$upper<-p$predicted_reference_HRDsum+q

# Be explicit that conformal coverage is marginal and assumes the new samples
# are exchangeable with the TCGA calibration set. For pediatric or
# different-platform data that assumption does not hold, so the guarantee does
# not transfer - hence "no_domain_shift_guarantee".
p$interval_status<-if(is.finite(q))"empirical_no_domain_shift_guarantee" else "unavailable_insufficient_calibration"

# Human-readable reason per row. Note the ordering: qc_fail (too many missing
# probes) takes precedence over ood (inside the assay but outside the training
# distribution). Even a passing sample gets a warning string, because no
# prediction from this research model is validated for a new domain.
p$warning<-ifelse(p$qc_fail,"excessive_missing_features",ifelse(p$ood,"outside_training_distance_reference","research_prediction_unvalidated_domain"))

# Two-column design, deliberate: the raw prediction is retained for analysis
# while estimate_for_display is blanked unless the sample is reportable. Display
# layers must read estimate_for_display, never the raw column.
p$estimate_for_display<-ifelse(p$reportable,p$predicted_reference_HRDsum,NA_real_)
# Bounds are withheld under the same rule, and additionally whenever the
# calibration quantile is not finite.
p$lower[!p$reportable|!is.finite(q)]<-NA_real_;p$upper[!p$reportable|!is.finite(q)]<-NA_real_

# NOTE: there is deliberately no HRD_high_probability column. An always-NA
# column invites a reader to assume a validated clinical threshold exists. The
# exploratory reference cut of 42 in config/analysis_protocol.json is a
# property of the TCGA reference distribution, not a pediatric clinical
# threshold, and a point regression prediction cannot be converted into a
# calibrated probability. Reintroduce this column only when a probability
# mapping has been fitted on separate calibration data and validated
# (docs/06_VALIDATION_STRATEGY.md, "Binary secondary").

# Provenance: record which model file produced these numbers so a result can
# always be traced back to a specific frozen artefact, and carry the matrix
# class through so the display layer cannot present a fixture as a result.
# app/app.R validates KIDS26_DEMO_RESULTS structurally only and prints
# `provenance` verbatim, so the class must live in that column, not beside it.
p$model_md5<-unname(tools::md5sum(args[1]))
p$matrix_provenance_class<-gate$class
p$probe_set_match<-probe_match
p$provenance<-if(isTRUE(gate$fatal))
  paste0("NOT A SCIENTIFIC RESULT - ",gate$class,"; scored under --allow-fixture") else
  paste0("frozen_research_model;domain_validity_requires_review;matrix=",gate$class,
         ";probe_set=",probe_match)

# A fixture run must not be mistaken for a reportable one downstream, whatever
# the per-sample QC happened to say.
if(isTRUE(gate$fatal)){p$reportable<-FALSE;p$estimate_for_display<-NA_real_
  p$lower<-NA_real_;p$upper<-NA_real_;p$warning<-"engineering_fixture_not_a_result"}

dir.create(dirname(args[3]),recursive=TRUE,showWarnings=FALSE)
write.table(p,args[3],sep="\t",row.names=FALSE,quote=FALSE)

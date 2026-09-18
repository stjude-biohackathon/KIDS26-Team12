# =============================================================================
# scripts/predict_frozen.R - Apply a frozen model to new samples.
# =============================================================================
#
# USAGE
#   Rscript scripts/predict_frozen.R <model.rds> <beta.tsv> <output.tsv> \
#     [--allow-fixture] [--metadata=<master_samples.tsv>] [--locked-evaluation]
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
#   --metadata=PATH  OPTIONAL master sample table. Supplying it turns on the
#                    locked-partition guard below and adds a `partition` column
#                    to the output. Omitting it keeps every existing invocation
#                    working exactly as before.
#
#   --locked-evaluation  the explicit, auditable opt-in required to score any
#                    sample whose partition is "locked_CNS". Equivalent env var:
#                    KIDS26_LOCKED_EVALUATION=1 (matching the
#                    KIDS26_REQUIRE_ADJUDICATED_QC convention in
#                    train_baseline.R). Without it, locked samples are REFUSED.
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
# THE LOCKED-PARTITION GUARD (defect 2, audited 2026-09-18)
# ---------------------------------------------------------
# "Run this ONCE" was a paragraph of prose with nothing enforcing it: this
# script had no idea which partition a sample belonged to, and left no trace of
# having scored the locked set. Two mechanisms now exist:
#
#   1. OPT-IN. With --metadata supplied, scoring any locked_CNS sample requires
#      --locked-evaluation (or KIDS26_LOCKED_EVALUATION=1). Without it the run
#      FAILS before any prediction is written. The metadata argument is optional
#      so existing invocations keep working; the guard is only as strong as the
#      habit of passing it, which is why the locked run MUST pass it.
#   2. LEDGER. Every locked evaluation appends one row to
#      results/LOCKED_EVALUATION_LEDGER.tsv recording the timestamp, git commit,
#      model artefact sha256, input matrix sha256/provenance, sample counts and
#      output sha256. The file is APPEND-ONLY and never rewritten, so a second
#      scoring of the locked set is visible as a second row rather than being
#      silently indistinguishable from the first.
#
# The guard only reads sample_id / cancer_type / partition from the metadata.
# It never reads, prints or writes locked HRDsum label values.
#
# OUTPUT COLUMN GUIDE
# -------------------
#   predicted_reference_HRDsum      AUTHORITATIVE. The C2 nonnegativity rule has
#                               been applied: pmax(prediction, 0). This is the
#                               column analysts and every downstream consumer
#                               should use.
#   predicted_reference_HRDsum_raw  the RAW model output BEFORE clipping (and
#                               before any target back-transform). Kept so the
#                               clip is auditable and exactly reversible - the
#                               authoritative column is a pure function of this
#                               one. Never report it: it can be negative, which
#                               HRDsum cannot.
#   clipping_rule               the rule that was applied, written per row so a
#                               reader of the table alone can tell.
#   estimate_for_display        the authoritative value, blanked to NA unless
#                               reportable. Display layers (the Shiny app)
#                               should use this.
#   lower / upper               conformal interval, NA unless reportable and a
#                               finite calibration quantile exists.
#   reportable                  passed QC AND inside the training distribution.
#   warning                     why a prediction is being withheld, if it is.
#   interval_status             whether the interval is usable at all.
#   model_md5                   checksum of the model file, for provenance.
#   matrix_provenance_class     what the scored matrix was verified to be.
#   probe_set_match             whether its probe allowlist matches the model's.
#   partition                   development / locked_CNS, only when --metadata
#                               was supplied.
#
# KNOWN GAPS (audited 2026-09-16)
# --------------------------------
#   * CLIPPING - RESOLVED 2026-09-18 (blocker C2). This script previously
#     emitted the elastic net's unbounded output directly, so it shipped
#     negative HRDsum predictions for roughly 3.2% of samples (docs/24), even
#     though C2 ADOPTED clipping at 0 (docs/21) and loco_one_fold.R had applied
#     the back-transform since run 01. The transform pair now comes from
#     R/calibration.R - the same get_transform() object the fold scripts use -
#     rather than being re-spelled here, so the two paths cannot drift.
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
#        [--metadata=<master.tsv>] [--locked-evaluation]
args<-commandArgs(trailingOnly=TRUE)
allow_fixture<-"--allow-fixture" %in% args;args<-args[args!="--allow-fixture"]
# Locked-evaluation opt-in. Accepted as a CLI flag OR as an env var, matching
# the KIDS26_REQUIRE_ADJUDICATED_QC convention in train_baseline.R, so it works
# from an LSF submission script without editing the command line.
locked_opt_in<-("--locked-evaluation" %in% args)||nzchar(Sys.getenv("KIDS26_LOCKED_EVALUATION"))
args<-args[args!="--locked-evaluation"]
# Optional metadata table. Supplying it enables the partition guard; omitting it
# leaves every pre-existing three-argument invocation working unchanged.
meta_arg<-grep("^--metadata=",args,value=TRUE)
args<-args[!grepl("^--metadata=",args)]
if(length(meta_arg)>1)stop("Pass --metadata= at most once")
meta_path<-if(length(meta_arg))sub("^--metadata=","",meta_arg) else NA_character_
if(length(args)!=3)stop("Expected model, beta matrix, output [--allow-fixture] [--metadata=<master.tsv>] [--locked-evaluation]")
source("R/model.R");source("R/provenance.R")
# C2 clipping lives in R/calibration.R as a matched (forward, inverse) pair, and
# is reached through get_transform() - the SAME object scripts/loco_one_fold.R
# back-transforms with. Re-spelling pmax(p, 0) here instead would recreate the
# two-sources-of-truth problem this repo already has with the CNS lock.
source("R/calibration.R")
if(!requireNamespace("data.table",quietly=TRUE))stop("Install data.table")
# glmnet must be loaded for its predict() S3 METHOD to be registered. fit_en()
# happens to do this at fit time, so every training script gets it for free and
# the omission here was invisible: this inference-only path died with
# "no applicable method for 'predict' applied to an object of class 'elnet'"
# the moment it was run standalone. Found while testing the clipping fix.
if(!requireNamespace("glmnet",quietly=TRUE))stop("Install glmnet via scripts/setup.R")

# Append-only ledger of every locked evaluation (defect 2). Path is fixed rather
# than an argument: a ledger you can redirect is a ledger you can avoid.
LEDGER_PATH<-"results/LOCKED_EVALUATION_LEDGER.tsv"
LEDGER_COLS<-c("timestamp_utc","git_commit","git_dirty","model_path","model_sha256",
               "input_matrix","input_matrix_sha256","matrix_provenance_class",
               "matrix_sidecar_sha256","n_samples_scored","n_locked_CNS_scored",
               "output_path","output_sha256","clipping_rule","opt_in_source","user","host")

# sha256 of a file. digest is already a declared dependency (scripts/setup.R);
# fall back to the system tool, and record the failure rather than a blank if
# neither is available - "unavailable" in the ledger is honest, "" is not.
file_sha256<-function(path){
  if(is.null(path)||is.na(path)||!nzchar(path)||!file.exists(path))return("file_absent")
  if(requireNamespace("digest",quietly=TRUE))
    return(digest::digest(file=path,algo="sha256"))
  out<-suppressWarnings(try(system2("sha256sum",shQuote(path),stdout=TRUE,stderr=FALSE),silent=TRUE))
  if(inherits(out,"try-error")||!length(out))return("sha256_unavailable")
  sub("\\s.*$","",out[1])
}

git_field<-function(args){
  out<-suppressWarnings(try(system2("git",args,stdout=TRUE,stderr=FALSE),silent=TRUE))
  if(inherits(out,"try-error")||!length(out))return(NA_character_)
  paste(out,collapse=" ")
}

# --- Provenance gate (B7) --------------------------------------------------
# Runs BEFORE the matrix is read, so a fixture costs nothing and an unlabelled
# 28 GB matrix is not loaded only to be rejected afterwards. train_baseline.R
# applies the same gate at fit time; this is the inference half.
matrix_pr<-read_matrix_provenance(args[2])
gate<-gate_matrix_provenance(matrix_pr,allow_fixture=allow_fixture)

# --- Locked-partition metadata (defect 2) ----------------------------------
# Read BEFORE the (potentially 28 GB) matrix so a governance refusal costs
# nothing. Only the three governance columns are touched; HRDsum and the other
# label columns are dropped immediately and are never read, printed or written.
meta<-NULL
if(!is.na(meta_path)){
  if(!file.exists(meta_path))stop("--metadata file not found: ",meta_path)
  meta_full<-read.delim(meta_path,check.names=FALSE,stringsAsFactors=FALSE)
  need<-c("sample_id","cancer_type","partition")
  if(!all(need %in% names(meta_full)))
    stop("--metadata table must have columns: ",paste(need,collapse=", "),
         ". Rebuild it with scripts/build_master.py.")
  meta<-meta_full[,need,drop=FALSE];rm(meta_full)
  # Defect 3 again, on the inference side: the hardcoded CNS list and the
  # `partition` column must agree before either is trusted to define the lock.
  assert_partition_matches_cns(meta$cancer_type,meta$partition,
                               context=basename(meta_path))
}

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

# --- THE LOCKED-PARTITION GUARD (defect 2) ---------------------------------
# Runs before predict_en(), so a refused run writes nothing at all. Samples the
# metadata does not describe are an error rather than an implicit "development":
# an unknown sample is exactly the case where the guard must not assume safety.
scored_partition<-rep(NA_character_,nrow(x))
n_locked<-0L
if(!is.null(meta)){
  idx<-match(rownames(x),meta$sample_id)
  if(anyNA(idx))
    stop(sprintf(paste0("%d of %d samples in the matrix are absent from --metadata (%s), so their ",
                        "partition is unknown. Refusing to score samples whose lock status cannot be ",
                        "established."),sum(is.na(idx)),nrow(x),basename(meta_path)))
  scored_partition<-meta$partition[idx]
  locked<-scored_partition==CNS_LOCKED_LABEL
  n_locked<-sum(locked)
  if(n_locked>0L&&!locked_opt_in)
    stop(sprintf(paste0("REFUSING TO SCORE THE LOCKED PARTITION.\n",
                        "  %d of %d samples have partition == '%s' (cancer types: %s).\n",
                        "  This set is the project's one-shot external test. Scoring it is a\n",
                        "  decision, not a default, and every run erodes its held-out status.\n",
                        "  To proceed deliberately, pass --locked-evaluation (or set\n",
                        "  KIDS26_LOCKED_EVALUATION=1). The run will then be recorded in\n",
                        "  %s."),
                 n_locked,nrow(x),CNS_LOCKED_LABEL,
                 paste(sort(unique(meta$cancer_type[idx][locked])),collapse=", "),
                 LEDGER_PATH),call.=FALSE)
  if(n_locked>0L)
    message(sprintf("LOCKED EVALUATION: scoring %d locked_CNS samples under explicit opt-in. This will be appended to %s.",
                    n_locked,LEDGER_PATH))
}else if(locked_opt_in){
  # Opting in without metadata cannot be honoured: nothing here knows which
  # samples are locked. Say so rather than letting the flag imply a check ran.
  stop("--locked-evaluation was given but no --metadata=<master.tsv>, so partition ",
       "cannot be determined and no locked evaluation can be verified or recorded.",
       call.=FALSE)
}

# Apply the frozen model. predict_en() returns the raw prediction plus the
# missing_fraction / ood / qc_fail / reportable flags.
p<-predict_en(b,x);q<-b$interval_q
# A model frozen without a calibration step has no interval half-width; Inf
# propagates to NA bounds below rather than fabricating a narrow interval.
if(is.null(q))q<-Inf

# --- C2: back-transform and clip (defect 1) ---------------------------------
# predict_en() returns the elastic net's UNBOUNDED output. HRDsum = HRD_LOH +
# LST + TAI is a count-like scar burden bounded below at zero, and this script
# was shipping negative values for it.
#
# Consistency with scripts/loco_one_fold.R is the point: that script does
#     p$predicted_reference_HRDsum <- tf$inverse(p$predicted_reference_HRDsum)
# immediately after predict_en(), using a transform pair from R/calibration.R.
# The identical call is made here, with the transform taken from the frozen
# bundle's own `target_transform` when it recorded one - so a model fitted on
# log1p is back-transformed with expm1 here exactly as it is in the folds -
# falling back to "clip" for bundles frozen before that field existed (run 01
# was fitted on the identity scale, whose back-transform IS the clip). Both
# "clip" and "log1p" inverses end in pmax(., 0), so the adopted C2 rule holds
# either way.
#
# The pre-clip value is preserved in its own column. Clipping is deliberately
# the kind of change that must remain auditable: a reader can recompute the
# authoritative column from the raw one and confirm nothing else moved, and the
# raw column is what any future re-analysis of the clip decision needs.
bundle_transform<-if(is.null(b$target_transform))"clip" else b$target_transform
tf<-get_transform(bundle_transform)
p$predicted_reference_HRDsum_raw<-p$predicted_reference_HRDsum
p$predicted_reference_HRDsum<-tf$inverse(p$predicted_reference_HRDsum_raw)
# Belt and braces, and the executable statement of the adopted rule: whatever
# the bundle's transform was, the shipped column is non-negative.
p$predicted_reference_HRDsum<-pmax(p$predicted_reference_HRDsum,0)
n_clipped<-sum(is.finite(p$predicted_reference_HRDsum_raw)&p$predicted_reference_HRDsum_raw<0)
clipping_rule<-sprintf("C2_ADOPTED:pmax(inverse_%s(raw),0)",bundle_transform)
p$clipping_rule<-clipping_rule
if(n_clipped>0L)
  message(sprintf("C2 clipping: %d of %d raw predictions were negative and were clipped to 0 (rule: %s). Pre-clip values retained in predicted_reference_HRDsum_raw.",
                  n_clipped,nrow(p),clipping_rule))

# Symmetric conformal interval around the point estimate. Built from the CLIPPED
# (authoritative) value, so the interval is centred on the number that is
# actually reported; the lower bound is then clipped for the same reason the
# point estimate is - a negative HRDsum bound is not a statement about HRD.
p$lower<-pmax(p$predicted_reference_HRDsum-q,0);p$upper<-p$predicted_reference_HRDsum+q

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
# Partition travels with every row when it is known, so a reader of the table
# alone can tell whether they are looking at locked output.
p$partition<-if(is.null(meta))"not_supplied" else scored_partition
p$provenance<-if(isTRUE(gate$fatal))
  paste0("NOT A SCIENTIFIC RESULT - ",gate$class,"; scored under --allow-fixture") else
  paste0("frozen_research_model;domain_validity_requires_review;matrix=",gate$class,
         ";probe_set=",probe_match,";",clipping_rule,
         if(n_locked>0L)";LOCKED_EVALUATION_RECORDED" else "")

# A fixture run must not be mistaken for a reportable one downstream, whatever
# the per-sample QC happened to say.
if(isTRUE(gate$fatal)){p$reportable<-FALSE;p$estimate_for_display<-NA_real_
  p$lower<-NA_real_;p$upper<-NA_real_;p$warning<-"engineering_fixture_not_a_result"}

dir.create(dirname(args[3]),recursive=TRUE,showWarnings=FALSE)
write.table(p,args[3],sep="\t",row.names=FALSE,quote=FALSE)

# =============================================================================
# APPEND-ONLY LEDGER OF LOCKED EVALUATIONS (defect 2)
# =============================================================================
# Written AFTER the output exists, so the output's own sha256 can be recorded
# and the ledger row describes an artefact that is actually on disk.
#
# APPEND-ONLY is the whole design. Nothing here reads, rewrites or deduplicates
# the existing file: a second scoring of the locked set appends a second row.
# Repeated evaluation of a one-shot test set is a real scientific cost, and the
# purpose of this file is to make that cost visible rather than to prevent it -
# a guard that can be argued past silently is worse than one that leaves a
# receipt. Note that `results/` is gitignored, so this ledger lives on the
# filesystem alongside the run; quote its rows in any write-up of a locked
# evaluation.
if(n_locked>0L){
  row<-data.frame(
    timestamp_utc           = format(as.POSIXct(Sys.time(),tz="UTC"),"%Y-%m-%dT%H:%M:%SZ"),
    git_commit              = git_field(c("rev-parse","HEAD")),
    git_dirty               = if(is.na(git_field(c("status","--porcelain")))) NA_character_
                              else as.character(nzchar(git_field(c("status","--porcelain")))),
    model_path              = args[1],
    model_sha256            = file_sha256(args[1]),
    input_matrix            = args[2],
    # The matrix sha256 is recorded from the sidecar when it carries one: the
    # locked matrix is up to 28 GB and re-hashing it here would add tens of
    # minutes to every run. The sidecar hash IS the provenance record, and
    # matrix_sidecar_sha256 names where the value came from.
    input_matrix_sha256     = if(!is.null(matrix_pr)&&!is.null(matrix_pr$sha256))
                                as.character(matrix_pr$sha256) else file_sha256(args[2]),
    matrix_provenance_class = gate$class,
    matrix_sidecar_sha256   = file_sha256(provenance_path(args[2])),
    n_samples_scored        = nrow(p),
    n_locked_CNS_scored     = n_locked,
    output_path             = args[3],
    output_sha256           = file_sha256(args[3]),
    clipping_rule           = clipping_rule,
    opt_in_source           = if("--locked-evaluation" %in% commandArgs(trailingOnly=TRUE))
                                "cli:--locked-evaluation" else "env:KIDS26_LOCKED_EVALUATION",
    user                    = Sys.info()[["user"]],
    host                    = Sys.info()[["nodename"]],
    stringsAsFactors        = FALSE)
  row<-row[,LEDGER_COLS,drop=FALSE]
  dir.create(dirname(LEDGER_PATH),recursive=TRUE,showWarnings=FALSE)
  # Header only when the file does not yet exist; thereafter pure append.
  new_file<-!file.exists(LEDGER_PATH)
  write.table(row,LEDGER_PATH,sep="\t",row.names=FALSE,quote=FALSE,
              col.names=new_file,append=!new_file)
  n_prior<-length(readLines(LEDGER_PATH,warn=FALSE))-2L
  message(sprintf("Locked evaluation appended to %s (%d locked samples).",LEDGER_PATH,n_locked))
  if(n_prior>0L)
    warning(sprintf(paste0("%s already contained %d earlier locked evaluation(s). The locked set is ",
                           "a ONE-SHOT external test; this is at least evaluation number %d and its ",
                           "held-out status is correspondingly weaker."),
                    LEDGER_PATH,n_prior,n_prior+1L),call.=FALSE)
}

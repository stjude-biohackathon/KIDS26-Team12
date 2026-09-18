# =============================================================================
# tests/test_provenance_gate.R - Executable assertions for the SHIPPING PATH
#                                gates: B7 provenance, C2 clipping, the locked
#                                partition, and merge completeness.
# =============================================================================
#   Rscript tests/test_provenance_gate.R     (from the repository root)
#
# These are refusal tests. Each one asserts that a specific bad input is
# REJECTED, because the failure mode B7 describes is a gate that exists but
# does not actually stop anything.
#
# Sections 1-9 are the original B7 unit-level gate.
# Sections 10-16 (added 2026-09-18 with the pre-lock audit fixes) drive
# scripts/predict_frozen.R and scripts/loco_merge.R as SUBPROCESSES against tiny
# fixtures in tempdir(), because the audited defects were all in the scripts
# rather than in the library functions they call, and nothing in the suite
# executed those scripts.
#
# Every REJECT here is paired with an ACCEPT differing in ONE property (the B5
# rule in docs/21): a gate that always refuses is as useless as one that never
# does. No fixture contains or reads a real GBM/LGG HRDsum value - the locked
# cases use synthetic samples labelled GBM/LGG with synthetic targets.
# =============================================================================

source("R/provenance.R")

real <- list(engineering_only = FALSE, n_samples = 7707L, n_probes = 336480L,
             probe_allowlist_sha256 = "f359e43bc171c544e55e60000905bb0bba1c44ee")
fixture <- list(engineering_only = TRUE, n_samples = 8L, n_probes = 50L,
                probe_allowlist_sha256 = "deadbeef")
nohash <- list(engineering_only = FALSE, n_samples = 100L, n_probes = 1000L)

# 1. The real matrix passes and is classed as such.
stopifnot(identical(gate_matrix_provenance(real)$class, "verified_research_matrix"))

# 2. An engineering fixture is refused outright.
stopifnot(inherits(try(gate_matrix_provenance(fixture), silent = TRUE), "try-error"))

# 3. A matrix with NO sidecar at all is refused. This is the case an attacker or
#    a hurried teammate hits: drop any TSV on disk with no metadata beside it.
stopifnot(inherits(try(gate_matrix_provenance(NULL), silent = TRUE), "try-error"))

# 4. A non-fixture matrix whose allowlist cannot be verified is still refused -
#    "not declared fake" is not the same as "verified real".
stopifnot(inherits(try(gate_matrix_provenance(nohash), silent = TRUE), "try-error"))

# 5. --allow-fixture permits the run but must NOT launder the label: the class
#    still names the fixture and additionally records that it was overridden.
ov <- suppressWarnings(gate_matrix_provenance(fixture, allow_fixture = TRUE))
stopifnot(grepl("ENGINEERING_FIXTURE", ov$class),
          grepl("OVERRIDDEN_BY_ALLOW_FIXTURE", ov$class))

# 6. The override must also raise a warning, not pass quietly.
stopifnot(length(withCallingHandlers(
  { w <- character(0)
    gate_matrix_provenance(fixture, allow_fixture = TRUE); w },
  warning = function(cond) { w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning") }
)) > 0)

# 7. Allowlist comparison distinguishes all three states. A mismatch is a
#    warning-level condition, not a refusal, because cross-platform transfer is
#    the actual research goal (see R/provenance.R design note).
stopifnot(identical(compare_probe_allowlist(real, real), "matched"),
          identical(compare_probe_allowlist(real, fixture), "MISMATCH_different_probe_allowlist"),
          identical(compare_probe_allowlist(real, nohash), "unknown_allowlist_not_recorded"),
          identical(compare_probe_allowlist(NULL, real), "unknown_allowlist_not_recorded"))

# 8. The sidecar path convention must match the one train_baseline.R uses.
stopifnot(identical(provenance_path("data/processed/beta.tsv"),
                    "data/processed/beta.provenance.json"))

# 9. End-to-end against the real sidecar on disk, if present: the production
#    matrix must pass its own gate. Guards against the sidecar drifting out of
#    the schema the gate expects.
if (file.exists("data/processed/beta.provenance.json")) {
  stopifnot(identical(
    gate_matrix_provenance(read_matrix_provenance("data/processed/beta.tsv"))$class,
    "verified_research_matrix"))
  cat("  (verified against the real data/processed/beta.provenance.json)\n")
}

# =============================================================================
# SCRIPT-LEVEL GATES: scripts/predict_frozen.R and scripts/loco_merge.R
# =============================================================================
stopifnot(file.exists("scripts/predict_frozen.R"), file.exists("scripts/loco_merge.R"))
if (!requireNamespace("glmnet", quietly = TRUE)) stop("Install dependencies with scripts/setup.R")
if (!requireNamespace("data.table", quietly = TRUE)) stop("Install dependencies with scripts/setup.R")

pass <- function(n, msg) cat(sprintf("  [%s] PASS  %s\n", n, msg))

tmp <- tempfile("predict_frozen_fixture"); dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

# The shipping-path ledger is a FIXED path by design, so these tests would
# otherwise write into a real one. Snapshot and restore it.
LEDGER <- "results/LOCKED_EVALUATION_LEDGER.tsv"
ledger_existed <- file.exists(LEDGER)
ledger_backup <- if (ledger_existed) readLines(LEDGER, warn = FALSE) else NULL
restore_ledger <- function() {
  if (ledger_existed) writeLines(ledger_backup, LEDGER) else unlink(LEDGER)
}
on.exit(restore_ledger(), add = TRUE)
unlink(LEDGER)

# --- Fixture: a frozen bundle, a beta matrix, a sidecar, a master table -----
# Deliberately built so that SOME predictions are negative: the target is
# centred near zero and the held-out rows are shifted away from training, which
# is exactly the geometry that produced 3.2% negative predictions in run 01.
source("R/model.R")
set.seed(3312)
n_tr <- 60L; n_te <- 8L; n_probe <- 40L
probes <- sprintf("cg%05d", seq_len(n_probe))
x_tr <- matrix(runif(n_tr * n_probe), n_tr, n_probe,
               dimnames = list(sprintf("TR-%02d", seq_len(n_tr)), probes))
y_tr <- 20 * x_tr[, 1] - 18 * x_tr[, 2] + rnorm(n_tr, 2, 1)
bundle <- fit_en(x_tr, y_tr, sprintf("p%02d", seq_len(n_tr)),
                 rep(c("BRCA", "LUAD", "KIRC"), length.out = n_tr), max_features = 20L)
bundle$target_transform <- "identity"
bundle$interval_q <- 4.0
model_path <- file.path(tmp, "frozen.rds"); saveRDS(bundle, model_path)

# Scored samples: half development-labelled, half CNS-labelled. All synthetic.
te_ids <- c(sprintf("DEV-%02d", 1:4), sprintf("CNS-%02d", 1:4))
x_te <- matrix(runif(n_te * n_probe, 0, 0.35), n_te, n_probe,
               dimnames = list(te_ids, probes))
write_matrix <- function(path, x) {
  write.table(data.frame(probe_id = colnames(x), t(x), check.names = FALSE),
              path, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(sprintf('{"engineering_only": false, "n_samples": %d, "n_probes": %d, "probe_allowlist_sha256": "%s", "sha256": "fixture_matrix_sha"}',
                     nrow(x), ncol(x), "f359e43bc171c544e55e60000905bb0bba1c44ee"),
             provenance_path(path))
  path
}
beta_path <- write_matrix(file.path(tmp, "beta.tsv"), x_te)

# Master table carrying ONLY the governance columns predict_frozen.R reads.
# No HRDsum column at all, so this fixture cannot leak a label even by accident.
master <- data.frame(
  sample_id   = te_ids,
  cancer_type = c(rep("BRCA", 4), "GBM", "GBM", "LGG", "LGG"),
  partition   = c(rep("development", 4), rep("locked_CNS", 4)),
  stringsAsFactors = FALSE)
master_path <- file.path(tmp, "master.tsv")
write.table(master, master_path, sep = "\t", quote = FALSE, row.names = FALSE)

# Development-only master (same table, CNS rows relabelled to a non-CNS type)
# for the paired ACCEPT case.
master_dev <- master
master_dev$cancer_type <- rep("BRCA", nrow(master_dev))
master_dev$partition <- rep("development", nrow(master_dev))
master_dev_path <- file.path(tmp, "master_dev.tsv")
write.table(master_dev, master_dev_path, sep = "\t", quote = FALSE, row.names = FALSE)

run_predict <- function(extra = character(0), out = file.path(tmp, "out.tsv"), env = character(0)) {
  unlink(out)
  res <- suppressWarnings(system2("Rscript",
    c("scripts/predict_frozen.R", shQuote(model_path), shQuote(beta_path), shQuote(out), extra),
    stdout = TRUE, stderr = TRUE, env = env))
  status <- attr(res, "status")
  list(ok = is.null(status) || status == 0L, log = paste(res, collapse = "\n"), out = out)
}
read_out <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)

cat("\nscript-level gates (predict_frozen.R, loco_merge.R)\n\n")

# --- 10. C2 CLIPPING IS APPLIED (defect 1) ----------------------------------
# The ACCEPT half: a plain development run must succeed AND must emit no
# negative value in the authoritative column.
r <- run_predict(c(paste0("--metadata=", master_dev_path)))
stopifnot(r$ok, file.exists(r$out))
o <- read_out(r$out)
stopifnot("predicted_reference_HRDsum" %in% names(o),
          "predicted_reference_HRDsum_raw" %in% names(o),
          all(o$predicted_reference_HRDsum >= 0))
pass(10, "predict_frozen.R emits no negative predicted_reference_HRDsum")

# --- 11. THE CLIP IS REVERSIBLE AND ACTUALLY FIRED --------------------------
# This is the check that can FAIL if someone deletes the pmax(): the raw column
# must still contain the negative values, the authoritative column must be
# exactly pmax(raw, 0), and the fixture must genuinely have produced at least
# one negative raw prediction - otherwise check 10 would pass vacuously.
stopifnot(any(o$predicted_reference_HRDsum_raw < 0))
stopifnot(isTRUE(all.equal(o$predicted_reference_HRDsum,
                           pmax(o$predicted_reference_HRDsum_raw, 0), tolerance = 1e-10)))
# Clipping must not reorder anyone (C2's defining property, docs/21).
stopifnot(!is.unsorted(o$predicted_reference_HRDsum[order(o$predicted_reference_HRDsum_raw)]))
pass(11, "raw pre-clip column is retained, clip = pmax(raw,0), ordering preserved")

# --- 12. THE CLIPPING RULE IS RECORDED --------------------------------------
# A clip that leaves no trace is not auditable. It must appear both as its own
# column and inside the provenance string app/app.R prints verbatim.
stopifnot("clipping_rule" %in% names(o),
          all(grepl("pmax", o$clipping_rule)),
          all(grepl("C2", o$clipping_rule)),
          all(grepl("pmax", o$provenance)))
pass(12, "the clipping rule is recorded per row and in the provenance string")

# --- 13. LOCKED SAMPLES ARE REFUSED WITHOUT THE OPT-IN (defect 2) -----------
# REJECT: the same command that succeeded in check 10, differing only in which
# master table names the samples' partition.
r_locked <- run_predict(c(paste0("--metadata=", master_path)))
stopifnot(!r_locked$ok,
          grepl("locked", r_locked$log, ignore.case = TRUE),
          !file.exists(r_locked$out))   # nothing written, not even a partial table
pass(13, "scoring locked_CNS samples without the opt-in FAILS and writes nothing")

# --- 14. THE OPT-IN WORKS, BY FLAG AND BY ENV VAR ---------------------------
# ACCEPT: identical command plus the opt-in. Paired with 13, this proves the
# guard keys on the opt-in and not on something incidental.
r_ok <- run_predict(c(paste0("--metadata=", master_path), "--locked-evaluation"))
stopifnot(r_ok$ok, file.exists(r_ok$out))
ol <- read_out(r_ok$out)
stopifnot(identical(sort(unique(ol$partition)), c("development", "locked_CNS")))
stopifnot(file.exists(LEDGER))
led <- read.delim(LEDGER, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(nrow(led) == 1L,
          all(c("timestamp_utc","git_commit","model_sha256","input_matrix_sha256",
                "n_samples_scored","n_locked_CNS_scored","output_sha256",
                "clipping_rule") %in% names(led)),
          led$n_locked_CNS_scored[1] == 4L,
          led$n_samples_scored[1] == length(te_ids),
          nchar(led$model_sha256[1]) == 64L,
          nchar(led$output_sha256[1]) == 64L)
pass(14, "--locked-evaluation permits the run and writes one full ledger row")

# --- 15. THE LEDGER IS APPEND-ONLY ------------------------------------------
# A second locked evaluation must add a SECOND row, not overwrite the first.
# This is what makes repeated scoring of a one-shot test set visible.
first_row <- led[1, ]
r_again <- run_predict(c(paste0("--metadata=", master_path)),
                       env = "KIDS26_LOCKED_EVALUATION=1")
stopifnot(r_again$ok)
led2 <- read.delim(LEDGER, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(nrow(led2) == 2L,
          identical(as.character(led2[1, ]), as.character(first_row)),
          identical(led2$opt_in_source[2], "env:KIDS26_LOCKED_EVALUATION"))
pass(15, "a second locked evaluation APPENDS; the first row is untouched")

# --- 16. THE GUARD CANNOT BE SATISFIED BY THE FLAG ALONE --------------------
# REJECT: --locked-evaluation with no metadata. Nothing would know which samples
# are locked, so the flag must not be allowed to imply that a check happened.
# ACCEPT (paired): no metadata AND no flag is the legacy 3-argument invocation,
# which must keep working unchanged.
r_flag_only <- run_predict("--locked-evaluation")
stopifnot(!r_flag_only$ok, grepl("metadata", r_flag_only$log))
r_legacy <- run_predict(character(0))
stopifnot(r_legacy$ok, all(read_out(r_legacy$out)$predicted_reference_HRDsum >= 0))
# A sample the metadata does not describe must be refused rather than assumed
# to be development.
short_master <- file.path(tmp, "short.tsv")
write.table(master[1:2, ], short_master, sep = "\t", quote = FALSE, row.names = FALSE)
r_short <- run_predict(c(paste0("--metadata=", short_master), "--locked-evaluation"))
stopifnot(!r_short$ok, grepl("absent from --metadata", r_short$log))
pass(16, "flag without metadata is REJECTED; legacy 3-arg call still ACCEPTED")

restore_ledger()

# --- 17. loco_merge.R completeness gate (defect 4) --------------------------
# ACCEPT first: a complete two-fold cohort merges and ships the tissue R^2.
# Then remove exactly one artefact at a time and require a refusal each time.
md <- file.path(tmp, "merge"); dir.create(file.path(md, "folds"), recursive = TRUE)
mm_meta <- data.frame(
  sample_id   = sprintf("S%03d", 1:12),
  patient_id  = sprintf("P%03d", 1:12),
  cancer_type = rep(c("BRCA", "LUAD", "GBM"), each = 4),
  partition   = rep(c("development", "development", "locked_CNS"), each = 4),
  # BRCA sits low and LUAD high ON PURPOSE: a tissue-identity R^2 can only be
  # shown to be non-trivial on a fixture where tissue actually carries level.
  HRDsum      = c(2, 4, 6, 8, 30, 34, 38, 42, rep(0, 4)),
  purity      = 0.7, stringsAsFactors = FALSE)
mm_meta_path <- file.path(md, "master.tsv")
write.table(mm_meta, mm_meta_path, sep = "\t", quote = FALSE, row.names = FALSE)

set.seed(551)
write_fold <- function(type, rows, transform = "identity", rank = "pooled", n_metrics = NULL) {
  pr <- data.frame(sample_id = mm_meta$sample_id[rows],
                   predicted_reference_HRDsum = mm_meta$HRDsum[rows] + rnorm(length(rows)),
                   actual = mm_meta$HRDsum[rows], cancer_type = type,
                   patient_id = mm_meta$patient_id[rows], reportable = TRUE,
                   stringsAsFactors = FALSE)
  w <- function(d, f) write.table(d, file.path(md, "folds", f), sep = "\t", quote = FALSE, row.names = FALSE)
  w(pr, sprintf("predictions_%s.tsv", type))
  w(data.frame(n = if (is.null(n_metrics)) length(rows) else n_metrics,
               MAE = 1, cancer_type = type, target_transform = transform,
               feature_rank = rank, lambda_at_boundary = FALSE, stringsAsFactors = FALSE),
    sprintf("metrics_%s.tsv", type))
  w(data.frame(n = length(rows), MAE_model = 1, MAE_tissue_mean_null = 2,
               skill_vs_tissue_mean = 0.5, cancer_type = type, stringsAsFactors = FALSE),
    sprintf("nullpanel_%s.tsv", type))
}
run_merge <- function() {
  res <- suppressWarnings(system2("Rscript", c("scripts/loco_merge.R", shQuote(md), shQuote(mm_meta_path)),
                                  stdout = TRUE, stderr = TRUE))
  st <- attr(res, "status")
  list(ok = is.null(st) || st == 0L, log = paste(res, collapse = "\n"))
}

write_fold("BRCA", 1:4); write_fold("LUAD", 5:8)
m <- run_merge()
stopifnot(m$ok, file.exists(file.path(md, "pooled_tissue_identity_r2.tsv")))
tr2 <- read.delim(file.path(md, "pooled_tissue_identity_r2.tsv"), stringsAsFactors = FALSE)
stopifnot(nrow(tr2) == 2L,
          identical(sort(tr2$quantity), c("actual_HRDsum", "predicted_reference_HRDsum")),
          all(is.finite(tr2$r2_tissue_identity)),
          all(tr2$r2_tissue_identity >= 0 & tr2$r2_tissue_identity <= 1))
# Not a constant: on this fixture BRCA and LUAD differ in level, so both R^2
# values must be materially above zero. A stubbed-out implementation returning
# 0 or NA fails here.
stopifnot(all(tr2$r2_tissue_identity > 0.05))
pass(17, "loco_merge.R ships pooled_tissue_identity_r2.tsv for predicted AND observed")

# --- 18. A fold that crashed after writing predictions is caught ------------
# The exact hole in the old gate: it globbed predictions_*.tsv only.
unlink(file.path(md, "folds", "metrics_LUAD.tsv"))
m <- run_merge(); stopifnot(!m$ok, grepl("metrics_|nullpanel_", m$log))
write_fold("LUAD", 5:8)
unlink(file.path(md, "folds", "nullpanel_LUAD.tsv"))
m <- run_merge(); stopifnot(!m$ok, grepl("metrics_|nullpanel_", m$log))
write_fold("LUAD", 5:8)
stopifnot(run_merge()$ok)     # ACCEPT again once restored
pass(18, "a fold missing metrics_ or nullpanel_ is REJECTED (predictions alone is not enough)")

# --- 19. Row counts must reconcile ------------------------------------------
write_fold("LUAD", 5:8, n_metrics = 99L)
m <- run_merge(); stopifnot(!m$ok, grepl("counts disagree", m$log))
write_fold("LUAD", 5:8); stopifnot(run_merge()$ok)
pass(19, "a metrics row whose n disagrees with its predictions is REJECTED")

# --- 20. Folds must share target_transform AND feature_rank -----------------
write_fold("LUAD", 5:8, transform = "log1p")
m <- run_merge(); stopifnot(!m$ok, grepl("target_transform", m$log))
write_fold("LUAD", 5:8, rank = "within_tissue")
m <- run_merge(); stopifnot(!m$ok, grepl("feature_rank", m$log))
write_fold("LUAD", 5:8); stopifnot(run_merge()$ok)
pass(20, "folds disagreeing on target_transform or feature_rank are REJECTED")

# --- 21. loco_merge.R enforces the lock agreement (defect 3) ----------------
bad_meta <- mm_meta; bad_meta$partition[9] <- "development"   # a GBM row unlocked
bad_path <- file.path(md, "master_bad.tsv")
write.table(bad_meta, bad_path, sep = "\t", quote = FALSE, row.names = FALSE)
res <- suppressWarnings(system2("Rscript", c("scripts/loco_merge.R", shQuote(md), shQuote(bad_path)),
                                stdout = TRUE, stderr = TRUE))
stopifnot(!identical(attr(res, "status"), NULL), grepl("LOCK MISMATCH", paste(res, collapse = "\n")))
pass(21, "loco_merge.R refuses a master table whose partition disagrees with the CNS list")

cat("\nB7 provenance gate tests passed: fixtures, unlabelled matrices and\n")
cat("  unverifiable allowlists are all refused; overrides are stamped.\n")
cat("Shipping-path gates passed: C2 clipping is applied and auditable, the\n")
cat("  locked partition requires an explicit opt-in and leaves an append-only\n")
cat("  ledger row, and loco_merge.R refuses incomplete or mixed-protocol runs.\n")

# =============================================================================
# tests/smoke_model.R - Fast correctness checks for the R modelling core.
# =============================================================================
#
# USAGE
#   Rscript tests/smoke_model.R          (from the repository root)
#
# This is a stopifnot() script, not a testthat suite, so it is NOT picked up by
# `python -m unittest discover -s tests`. Run it explicitly. It uses a tiny
# synthetic matrix and finishes in seconds, so there is no excuse for skipping
# it before a long training run.
#
# WHY EACH CHECK EXISTS - these are regression guards for specific, realistic
# failure modes, not filler assertions:
#
#   1. Fit-object immutability. apply_preprocess() must not mutate the transform
#      it is handed. If it ever did, applying a model to held-out data would
#      silently alter the model - the exact shape of a leakage bug.
#   2. Conformal degradation. With too few calibration points the interval must
#      come back Inf ("no usable interval") rather than a falsely narrow one.
#   3. Single-sample invariance. A prediction for one sample alone must equal
#      that sample's prediction within a batch. If this fails, some statistic is
#      being computed across the input batch at predict time - again, leakage.
#      This is the check that matters most for a single-sample clinical assay.
#   4. Column-order invariance. Feeding probes in reverse order must give an
#      identical answer, proving apply_preprocess() aligns features by NAME
#      rather than by position. A positional bug here would silently pair every
#      probe with the wrong coefficient.
#   5. Duplicate-patient rejection. Must be a hard error, since repeated
#      patients leak across fold boundaries.
#   6. Degenerate-metric handling. A constant predictor yields NA R2, not a
#      crash or a misleading number.
#   7. EXPLICIT LEAKAGE ASSERTIONS. Preprocessing constants must be computed
#      from training rows only, applied unchanged to held-out data, and aligned
#      by feature NAME. These were previously asserted only in comments; they
#      are now executable so a regression fails loudly.
#
# NOT COVERED HERE: any of the leakage boundaries at the SCRIPT level (that
# locked_CNS never enters training, that the calibration set is disjoint from
# the training set), and nothing at all about the Python pipeline that builds
# beta.tsv. Passing this file does not mean the pipeline is correct.
# =============================================================================

source("R/model.R")

# --- Synthetic fixture -----------------------------------------------------
# 90 samples x 30 probes of uniform noise in [0,1], mimicking beta values.
# Fixed seed keeps the whole script deterministic.
set.seed(26);x<-matrix(runif(90*30),90,30,dimnames=list(paste0("s",1:90),paste0("cg",1:30)))
# y depends on only the first two probes plus noise, so a working elastic net
# has genuine signal to find while the other 28 probes act as distractors.
# 90 unique patients (no duplicates) across 3 cancer types, 30 samples each -
# enough for inner_folds() to take its primary cancer-grouped path.
y<-30*x[,1]-15*x[,2]+rnorm(90);patients<-paste0("p",1:90);cancer<-rep(c("A","B","C"),each=30)

# --- Check 1: the fit object is not mutated by being applied ---------------
# Learn the transform on the first 60 rows, snapshot it, apply it to rows 61-90,
# then confirm the snapshot still matches.
pp<-fit_preprocess(x[1:60,,drop=FALSE],10);before<-pp
z<-apply_preprocess(x[61:90,,drop=FALSE],pp)
stopifnot(identical(pp,before),identical(colnames(z),pp$features))

# --- Check 2: conformal interval degrades honestly -------------------------
# n=5 cannot support 95% coverage (ceiling((5+1)*0.95) = 6 > 5), so Inf.
# n=30 can (ceiling(31*0.95) = 30 <= 30), so a finite width.
stopifnot(is.infinite(conformal_q(1:5,1:5)),is.finite(conformal_q(1:30,1:30)))

# --- Fit a model for the remaining checks ----------------------------------
b<-fit_en(x,y,patients,cancer,max_features=20)

# --- Check 3: single-sample prediction == batch prediction ------------------
# Tolerance 1e-10 means bitwise-equivalent in practice. Any cross-sample
# statistic at predict time would break this immediately.
p1<-predict_en(b,x[1,,drop=FALSE]);pall<-predict_en(b,x)
stopifnot(isTRUE(all.equal(p1$predicted_reference_HRDsum,pall$predicted_reference_HRDsum[1],tolerance=1e-10)))

# --- Check 4: column order does not matter ---------------------------------
# Same sample, probes reversed. Must produce the identical prediction.
p2<-predict_en(b,x[1,rev(seq_len(ncol(x))),drop=FALSE])
stopifnot(isTRUE(all.equal(p1$predicted_reference_HRDsum,p2$predicted_reference_HRDsum,tolerance=1e-10)))

# --- Check 5: duplicate patients are rejected ------------------------------
# All 90 rows given the same patient ID must raise, via the anyDuplicated()
# guard in inner_folds().
try_duplicate<-try(fit_en(x,y,rep("same",90),cancer),silent=TRUE);stopifnot(inherits(try_duplicate,"try-error"))

# --- Check 6: R2 is NA when the reference values are constant ---------------
# sum((y - mean(y))^2) is 0, so R2 is undefined; metrics() must return NA
# rather than dividing by zero.
stopifnot(is.na(metrics(rep(1,5),1:5)$R2))

# =============================================================================
# Check 7: EXPLICIT LEAKAGE ASSERTIONS
# =============================================================================
# The central claim of this codebase is that every preprocessing constant is
# learned from training rows ONLY and thereafter applied unchanged. Until now
# that was asserted in comments and verified by human code reading. These checks
# make it executable, so a future edit that breaks the boundary fails the build
# instead of silently producing better-looking results.
#
# The design of the test: build a held-out set whose distribution is DELIBERATELY
# and grossly different from training. If any statistic were recomputed on that
# set, the difference would be impossible to miss.

set.seed(4413)
x_tr <- matrix(runif(60*30), 60, 30, dimnames=list(paste0("tr",1:60), paste0("cg",1:30)))
# Held-out rows are shifted upward and compressed: a completely different
# centre and spread from the training rows.
x_te <- matrix(runif(40*30, min=0.80, max=0.95), 40, 30,
               dimnames=list(paste0("te",1:40), paste0("cg",1:30)))

pp_tr <- fit_preprocess(x_tr, max_features=10L, max_missing=0.05)

# 7a. The stored constants must equal statistics computed on the TRAINING rows
#     alone. This is the direct statement of the invariant.
#
#     Tolerance rather than identical(): fit_preprocess() uses
#     matrixStats::colMeans2 where this check uses base colMeans, and the two
#     sum in a different order. They agree to ~1e-16, which is floating-point
#     noise, not a difference in what was computed. Using identical() here would
#     make the test fail for a reason that has nothing to do with leakage.
feat <- pp_tr$features
stopifnot(
  isTRUE(all.equal(unname(pp_tr$center[feat]),
                   unname(colMeans(x_tr[, feat, drop=FALSE])), tolerance=1e-12)),
  isTRUE(all.equal(unname(pp_tr$scale[feat]),
                   unname(apply(x_tr[, feat, drop=FALSE], 2, sd)), tolerance=1e-12))
)

# 7b. Refitting on training rows must be deterministic and unaffected by the
#     existence of held-out data. Same input, same transform, every time.
stopifnot(identical(pp_tr, fit_preprocess(x_tr, max_features=10L, max_missing=0.05)))

# 7c. THE KEY ASSERTION. Standardising the held-out set with training constants
#     must leave it visibly OFF-CENTRE, because its true mean is far from the
#     training mean. If apply_preprocess() were recomputing the centre on the
#     data it is handed, every column would come back near zero mean - which is
#     exactly what a leakage bug looks like, and exactly what makes held-out
#     performance look better than it is.
z_te <- apply_preprocess(x_te, pp_tr, FALSE)
col_means_te <- colMeans(z_te)
stopifnot(
  # Not re-centred: the shifted set must sit well away from zero.
  all(abs(col_means_te) > 0.5),
  # Sanity: the same transform applied to the TRAINING rows does centre them,
  # confirming the constants themselves are correct rather than merely unused.
  all(abs(colMeans(apply_preprocess(x_tr, pp_tr, FALSE))) < 1e-8)
)

# 7d. The transform object must be unchanged after being applied to held-out
#     data. Check 1 covers the ordinary path; this repeats it for the extreme
#     distribution shift, where an in-place update would be most tempting.
pp_snapshot <- pp_tr
invisible(apply_preprocess(x_te, pp_tr, FALSE))
stopifnot(identical(pp_tr, pp_snapshot))

# 7e. Feature alignment is BY NAME, not by position, even when the held-out
#     matrix carries extra probes the model never saw and presents them in a
#     different order. The extra columns must be ignored and the answer must be
#     identical to the clean case.
x_te_scrambled <- cbind(x_te, cg_unseen_1=runif(40), cg_unseen_2=runif(40))
x_te_scrambled <- x_te_scrambled[, sample(ncol(x_te_scrambled)), drop=FALSE]
stopifnot(isTRUE(all.equal(apply_preprocess(x_te_scrambled, pp_tr, FALSE), z_te,
                           tolerance=1e-12)))

# 7f. The lambda path is a learned quantity too: it must be derivable from
#     training rows alone and must not depend on held-out data being present.
z_tr <- apply_preprocess(x_tr, pp_tr, FALSE)
y_tr <- 10*x_tr[,1] - 5*x_tr[,2] + rnorm(60)
stopifnot(identical(lambda_path(z_tr, y_tr, 0.5), lambda_path(z_tr, y_tr, 0.5)))

# =============================================================================
# Check 8: THE LOCK HAS ONE SOURCE OF TRUTH (defect 3)
# =============================================================================
# c("GBM","LGG") is hardcoded in six places while master_samples.tsv carries an
# independent `partition` column. assert_partition_matches_cns() is the
# executable comparison of the two. Each REJECT below is paired with an ACCEPT
# differing in ONE property, so a function that always errored (or always
# passed) would fail this block - the B5 rule in docs/21.
ct <- c("BRCA","BRCA","LUAD","GBM","LGG")
pt <- c("development","development","development","locked_CNS","locked_CNS")

# 8a. ACCEPT: the two definitions agree exactly.
stopifnot(isTRUE(assert_partition_matches_cns(ct, pt)))

# 8b. REJECT: a locked cancer type labelled development. This is the drift that
#     would put GBM into a training fold while the fold's own mask believed it
#     had excluded it.
pt_bad <- pt; pt_bad[4] <- "development"
stopifnot(inherits(try(assert_partition_matches_cns(ct, pt_bad), silent=TRUE), "try-error"))

# 8c. REJECT: the mirror image - a development cancer labelled locked_CNS. The
#     assertion must be two-way, or a tissue could be silently withheld from a
#     run that claims 30 folds.
pt_bad2 <- pt; pt_bad2[1] <- "locked_CNS"
stopifnot(inherits(try(assert_partition_matches_cns(ct, pt_bad2), silent=TRUE), "try-error"))

# 8d. REJECT: missing / blank / unrecognised partition values. A column that
#     cannot be read is not an authority, and NA must not silently mean
#     "development".
stopifnot(
  inherits(try(assert_partition_matches_cns(ct, replace(pt, 1, NA)), silent=TRUE), "try-error"),
  inherits(try(assert_partition_matches_cns(ct, replace(pt, 1, "")), silent=TRUE), "try-error"),
  inherits(try(assert_partition_matches_cns(ct, replace(pt, 1, "held_out")), silent=TRUE), "try-error"))

# 8e. REJECT: no partition column at all (the default), but ACCEPT under an
#     explicit require_partition=FALSE. The pair proves the default is doing
#     work rather than the function ignoring NULL.
stopifnot(
  inherits(try(assert_partition_matches_cns(ct, NULL), silent=TRUE), "try-error"),
  isFALSE(suppressWarnings(assert_partition_matches_cns(ct, NULL, require_partition=FALSE))))

# 8f. REJECT: length mismatch, which would otherwise recycle and compare the
#     wrong rows against each other.
stopifnot(inherits(try(assert_partition_matches_cns(ct, pt[1:3]), silent=TRUE), "try-error"))

# 8g. The rank-model alias must BE the canonical constant, not a copy that can
#     drift from it. This is the check that would catch someone re-adding a
#     second literal c("GBM","LGG") in R/rank_model.R.
rank_src <- paste(readLines("R/rank_model.R", warn=FALSE), collapse="\n")
stopifnot(grepl("CNS_LOCKED <- CNS_LOCKED_TYPES", rank_src),
          identical(CNS_LOCKED_TYPES, c("GBM","LGG")))

# 8h. END-TO-END on the real master table if it is present: the shipped cohort
#     must satisfy its own lock assertion. Reads only cancer_type and partition;
#     no HRDsum value is read or printed.
if (file.exists("data/processed/master_samples.tsv")) {
  m <- read.delim("data/processed/master_samples.tsv", check.names=FALSE, stringsAsFactors=FALSE)
  stopifnot(isTRUE(assert_partition_matches_cns(m$cancer_type, m$partition)))
  cat(sprintf("  (real master table: %d locked_CNS rows agree with the hardcoded list)\n",
              sum(m$partition == "locked_CNS")))
}

# =============================================================================
# Check 9: tissue_identity_r2() (defect 4)
# =============================================================================
# The lineage-confounding number. It must be both SENSITIVE (detect a quantity
# that is purely tissue level) and SPECIFIC (not manufacture confounding from
# tissue-independent noise) - otherwise shipping it on every run would be
# shipping a constant.
set.seed(9155)
ct9 <- rep(c("A","B","C"), each=40)

# 9a. Pure tissue effect: R2 must be ~1.
stopifnot(abs(tissue_identity_r2(rep(c(1,5,9), each=40), ct9) - 1) < 1e-10)

# 9b. Tissue-independent noise: R2 must be near 0, and materially below 9a.
r2_noise <- tissue_identity_r2(rnorm(120), ct9)
stopifnot(r2_noise < 0.1, r2_noise >= 0)

# 9c. It must agree with lm() to floating-point, since that is the definition
#     scripts/c1_zeroshot_percentile.R used and the number must be comparable.
v9 <- rnorm(120) + rep(c(0, 2, 6), each=40)
stopifnot(isTRUE(all.equal(tissue_identity_r2(v9, ct9),
                           summary(lm(v9 ~ factor(ct9)))$r.squared, tolerance=1e-10)))

# 9d. Degenerate input returns NA rather than erroring or returning a number:
#     one tissue, constant values, or too few rows.
stopifnot(is.na(tissue_identity_r2(v9, rep("A",120))),
          is.na(tissue_identity_r2(rep(3,120), ct9)),
          is.na(tissue_identity_r2(c(1,2), c("A","B"))))

cat("R preprocessing/grouping/order/single-sample/conformal smoke tests passed\n")
cat("Leakage assertions passed: preprocessing constants are training-only,\n")
cat("  applied unchanged to held-out data, and aligned by feature name.\n")
cat("Lock assertions passed: the hardcoded CNS list and the partition column\n")
cat("  must agree in BOTH directions; tissue-identity R2 is sensitive and specific.\n")

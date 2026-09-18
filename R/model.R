# =============================================================================
# R/model.R - Fold-safe elastic-net helpers for predicting reference HRDsum
#             from DNA methylation beta values.
# =============================================================================
#
# WHAT THIS FILE IS
# -----------------
# The modelling core of the project. It is deliberately a library of plain
# functions with no side effects: nothing here reads a file, writes a file, or
# looks at the command line. Callers (scripts/train_baseline.R and
# scripts/predict_frozen.R) own all I/O. That separation is what makes the
# leakage discipline below auditable.
#
# THE CENTRAL INVARIANT - READ THIS BEFORE EDITING ANYTHING
# ---------------------------------------------------------
# Every quantity that is *learned from data* must be estimated from training
# rows only, stored inside a fit object, and then *applied* unchanged to
# held-out rows. In this file the learned quantities are:
#
#   1. which probes survive the missingness filter   (fit_preprocess: keep)
#   2. the per-probe median used for imputation      (fit_preprocess: med)
#   3. which probes survive variance ranking         (fit_preprocess: ii)
#   4. the per-probe centre and scale                (fit_preprocess: mu, s)
#   5. the elastic-net coefficients                  (fit_en: mod)
#   6. the out-of-distribution distance cutoff       (fit_en: ood_cut)
#
# All six live in the list returned by fit_preprocess()/fit_en(). apply_preprocess()
# and predict_en() only ever *consume* them. If you find yourself computing a
# mean, median, variance, or quantile inside apply_preprocess() or predict_en(),
# you have introduced test-set leakage and the reported performance becomes
# meaningless. Do not do it.
#
# A subtle corollary: glmnet is always called with standardize=FALSE. glmnet's
# default is standardize=TRUE, which would re-standardise using whatever matrix
# it was handed - including validation rows - silently undoing invariant (4).
# The FALSE is load-bearing, not stylistic.
#
# WHAT IS BEING PREDICTED, AND WHAT THAT DOES NOT MEAN
# -----------------------------------------------------
# The target is HRDsum = HRD_LOH + LST + TAI, a genomic-scar burden score
# derived from allele-aware data (SNP6 arrays). Methylation arrays cannot
# observe loss of heterozygosity or telomeric allelic imbalance directly, so
# this model is a *predictor of an independently measured reference score*, not
# a re-derivation of it. Nothing produced here is a clinical HRD assay.
#
# KNOWN LIMITATIONS OF THIS IMPLEMENTATION (audited 2026-09-16)
# --------------------------------------------------------------
#   * RESOLVED 2026-09-16: fit_preprocess() now uses matrixStats compiled
#     reductions instead of apply(z, 2, ...) over ~336k columns. See the
#     benchmark note in that function.
#   * The lambda grid is five hardcoded values spanning four decades and is not
#     anchored to glmnet's data-derived lambda.max. There is no check that the
#     selected lambda is interior to the grid, so a boundary optimum would be
#     accepted silently.
#   * Tissue and permutation controls now exist at the bottom of this file
#     (tissue_mean_null, permute_within_tissue, permutation_test_within_tissue,
#     null_panel, within_tissue_metrics) and are wired into train_baseline.R.
#     They are necessary but not sufficient: the permutation test reuses the
#     fitted predictions rather than refitting on permuted labels, so it asks
#     "do these predictions track the outcome within tissue better than chance"
#     and NOT "could the whole pipeline manufacture this from noise". The
#     refit-based version is the stronger control and is still outstanding.
#   * Purity and ploidy remain excluded from the feature set. Excluding a
#     confounder does not adjust for it; see docs/22 for the stratified and
#     purity-matched analyses that do.
# =============================================================================


# -----------------------------------------------------------------------------
# within_tissue_var(): pooled WITHIN-TISSUE variance per probe, blocked.
# -----------------------------------------------------------------------------
# MOVED HERE FROM R/rank_model.R so that BOTH the tissue-relative model (Phase C)
# and the ABSOLUTE-target model (V3-abs) can share one tested implementation.
# R/rank_model.R sources this file, so its callers are unchanged.
#
# Training rows only - the caller subsets before calling. Non-finite cells are
# filled with the supplied per-probe medians (the same training medians the
# preprocessing uses), so this never invents values and never looks at y.
#
# Identity used: within_SS = sum_t (Q_t - S_t^2 / n_t), and
# var_within = within_SS / (N - T), where S_t = colSums(x_t) and
# Q_t = colSums(x_t^2). Columns are processed in blocks so peak memory is
# O(N * block_size), independent of the total probe count. No residualised copy
# of the matrix is ever created.
within_tissue_var <- function(x, cancer, med = NULL, block_size = 20000L) {
  stopifnot(is.matrix(x), !is.null(colnames(x)), nrow(x) == length(cancer))
  cancer <- as.character(cancer)
  groups <- sort(unique(cancer))
  idx <- lapply(groups, function(g) which(cancer == g))
  nt <- vapply(idx, length, integer(1))
  N <- nrow(x); Tn <- length(groups)
  if (N - Tn < 1L) stop("Need more training rows than tissues for within-tissue variance")

  p <- ncol(x)
  out <- numeric(p); names(out) <- colnames(x)
  starts <- seq.int(1L, p, by = block_size)

  for (st in starts) {
    cols <- st:min(st + block_size - 1L, p)
    xb <- x[, cols, drop = FALSE]
    storage.mode(xb) <- "double"
    nf <- !is.finite(xb)
    if (any(nf)) {
      if (is.null(med)) stop("within_tissue_var(): non-finite values need `med`")
      mb <- med[colnames(xb)]
      if (any(!is.finite(mb))) stop("within_tissue_var(): non-finite median supplied")
      # Fill by column, touching only the affected columns.
      for (j in which(matrixStatsOrBase_colAnys(nf))) xb[nf[, j], j] <- mb[j]
    }
    ss <- numeric(length(cols))
    for (k in seq_along(idx)) {
      i <- idx[[k]]
      if (!length(i)) next
      xt <- xb[i, , drop = FALSE]
      S <- colSums(xt)
      Q <- colSums(xt * xt)
      ss <- ss + (Q - (S * S) / nt[k])
      rm(xt, S, Q)
    }
    out[cols] <- ss / (N - Tn)
    rm(xb, nf, ss)
  }
  # Numerical guard: the identity can return tiny negatives for constant columns.
  out[out < 0 & out > -1e-8] <- 0
  out
}

# Tiny helper so within_tissue_var() can use the compiled reduction when
# matrixStats is available, exactly as fit_preprocess() does, without depending
# on it.
matrixStatsOrBase_colAnys <- function(m) {
  if (requireNamespace("matrixStats", quietly = TRUE)) matrixStats::colAnys(m)
  else apply(m, 2, any)
}


# =============================================================================
# V4-lineage-penalized: supervised within-tissue HRD meta-association ranking.
# docs/27_V3_PREREGISTRATION.md section 6, declared before any V3 result existed.
# =============================================================================
#
# WHAT THIS IS
#   The V3-abs filter is unsupervised: it ranks probes by how much they VARY
#   inside tissues, never asking whether that variation tracks HRD. V4 replaces
#   it with a SUPERVISED statistic that is still explicitly de-lineaged:
#
#     score_j = |z_meta_j| * c_j  -  lambda_pen * l_j
#
#   z_meta_j  fixed-effect meta-analysis, across training tissues, of the
#             within-tissue Pearson correlation between probe j and HRDsum,
#             Fisher-z transformed and inverse-variance weighted (w_t = n_t-3).
#   c_j       weighted sign consistency in [0,1]: 1 when every training tissue
#             agrees on the direction, 0 when the weight splits evenly.
#   l_j       lineage discriminability, the ICC-like between/within variance
#             ratio from the SAME sums-of-squares decomposition used by
#             within_tissue_var(), robustly rescaled onto the units of the HRD
#             term so that a lambda_pen grid of order 1 is meaningful.
#
# LEAKAGE
#   Everything here is a function of the rows it is handed. fit_en() hands it
#   the outer-training block for the refit and each inner fold's own training
#   rows for tuning; there is no global precomputation and no reuse across
#   folds. Unlike "pooled" and "within_tissue" this statistic reads y, which is
#   exactly why the leakage assertion in tests/test_v4_lineage_penalized.R
#   corrupts the held-out LABELS as well as the held-out feature rows.
#
# WHY IT IS ONE PASS
#   Per column block and per tissue t we accumulate three compiled reductions:
#     S_t = colSums(x_t), Q_t = colSums(x_t^2), P_t = x_t' y_t^c.
#   From those alone:
#     SS_within(j)   = sum_t (Q_t - S_t^2/n_t)          <- the V3 decomposition
#     SS_total(j)    = sum_t Q_t - (sum_t S_t)^2 / N
#     r_{t,j}        = P_t / sqrt((Q_t - S_t^2/n_t) * SS(y_t))
#   so the lineage term costs NOTHING extra beyond the correlation sweep, and
#   the whole ranking is a single blocked pass over the matrix. Peak memory is
#   O(N * block_size) exactly as within_tissue_var(), which is why the same
#   block_size=20000 default is used.
#
# Standardising methylation and HRDsum within each tissue before correlating is
# a no-op for Pearson's r (it is location/scale invariant), so the correlation
# below IS the within-tissue standardised association docs/27 section 6 asks for.
lineage_meta_stats <- function(x, y, cancer, med = NULL, block_size = 20000L,
                               min_n = 4L, r_clamp = 1 - 1e-6) {
  stopifnot(is.matrix(x), !is.null(colnames(x)), nrow(x) == length(cancer),
            length(y) == nrow(x))
  if (!all(is.finite(y))) stop("lineage_meta_stats(): y must be finite")
  cancer <- as.character(cancer)
  groups <- sort(unique(cancer))
  idx <- lapply(groups, function(g) which(cancer == g))
  nt  <- vapply(idx, length, integer(1))
  N <- nrow(x); Tn <- length(groups)
  if (N - Tn < 1L) stop("Need more training rows than tissues for the meta statistic")

  # A Fisher-z has sampling variance 1/(n_t-3), so a tissue with n_t <= 3 has no
  # usable weight at all. Such tissues still contribute to the variance
  # decomposition (they are real rows) but carry zero weight in the meta term.
  w <- as.numeric(pmax(nt - 3L, 0L)); w[nt < min_n] <- 0
  sw <- sum(w)
  if (sw <= 0) stop("lineage_meta_stats(): no training tissue has enough samples for a weight")

  # Centre y WITHIN each tissue. This is the "standardise HRDsum within tissue"
  # step; the scale half cancels in the correlation.
  yc  <- lapply(idx, function(i) y[i] - mean(y[i]))
  ssy <- vapply(yc, function(v) sum(v * v), numeric(1))

  p <- ncol(x)
  z_meta <- numeric(p); cons <- numeric(p); i2 <- numeric(p)
  lin <- numeric(p); nonconst <- logical(p)
  starts <- seq.int(1L, p, by = block_size)
  dfQ <- max(sum(w > 0) - 1L, 0L)

  for (st in starts) {
    cols <- st:min(st + block_size - 1L, p)
    m <- length(cols)
    xb <- x[, cols, drop = FALSE]
    storage.mode(xb) <- "double"
    nf <- !is.finite(xb)
    if (any(nf)) {
      if (is.null(med)) stop("lineage_meta_stats(): non-finite values need `med`")
      mb <- med[colnames(xb)]
      if (any(!is.finite(mb))) stop("lineage_meta_stats(): non-finite median supplied")
      for (j in which(matrixStatsOrBase_colAnys(nf))) xb[nf[, j], j] <- mb[j]
    }

    zt   <- matrix(0, Tn, m)     # Fisher-z per tissue x probe for this block
    Stot <- numeric(m); Qtot <- numeric(m); wSS <- numeric(m)
    for (k in seq_along(idx)) {
      i <- idx[[k]]
      if (!length(i)) next
      xt <- xb[i, , drop = FALSE]
      S <- colSums(xt); Q <- colSums(xt * xt)
      Stot <- Stot + S; Qtot <- Qtot + Q
      ssx <- Q - (S * S) / nt[k]
      ssx[ssx < 0 & ssx > -1e-8] <- 0
      wSS <- wSS + ssx
      if (w[k] > 0 && ssy[k] > 0) {
        # crossprod is a single BLAS call: the cross-product of every probe in
        # the block with this tissue's centred y.
        P <- as.numeric(crossprod(xt, yc[[k]]))
        den <- sqrt(ssx * ssy[k])
        r <- ifelse(den > 0, P / den, 0)
        # Clamp before atanh so a probe that is a perfect linear function of y
        # inside one tissue yields a large-but-finite z instead of Inf.
        r <- pmin(pmax(r, -r_clamp), r_clamp)
        zt[k, ] <- atanh(r)
      }
      rm(xt, S, Q, ssx)
    }

    totSS <- Qtot - (Stot * Stot) / N
    totSS[totSS < 0 & totSS > -1e-8] <- 0
    nonconst[cols] <- totSS > 0

    wz   <- colSums(w * zt)                    # w recycles down each column
    z_meta[cols] <- wz / sqrt(sw)              # fixed-effect combined z
    cons[cols]   <- abs(colSums(w * sign(zt))) / sw
    zbar <- wz / sw
    Qc   <- colSums(w * (zt - rep(zbar, each = Tn))^2)   # Cochran's Q
    i2[cols] <- if (dfQ > 0) pmax(0, (Qc - dfQ) / pmax(Qc, .Machine$double.eps)) else 0
    lin[cols] <- ifelse(totSS > 0, pmin(1, pmax(0, 1 - wSS / totSS)), 0)
    rm(xb, nf, zt, Stot, Qtot, wSS, totSS, wz, zbar, Qc)
  }

  nm <- colnames(x)
  names(z_meta) <- nm; names(cons) <- nm; names(i2) <- nm
  names(lin) <- nm; names(nonconst) <- nm
  list(z_meta = z_meta, consistency = cons, heterogeneity_i2 = i2,
       lineage = lin, nonconstant = nonconst,
       tissues = groups, n_per_tissue = nt, weights = w)
}


# -----------------------------------------------------------------------------
# lineage_score(): turn the meta statistics into the ranking score.
# -----------------------------------------------------------------------------
# score_j = |z_meta_j| * c_j - lambda_pen * l_j_standardised
#
# The two terms live on different scales: |z_meta| * c is a weighted z, l is a
# variance ratio bounded in [0,1]. docs/27 section 6 requires l to be put on a
# comparable scale; we do that ROBUSTLY (median / MAD) so a handful of extreme
# probes cannot set the exchange rate, and then multiply by the MAD of the HRD
# term so that lambda_pen = 1 means "one robust SD of lineage costs one robust
# SD of HRD evidence". Subtracting the median of l shifts every score by the
# same constant and therefore cannot change the ranking; only the scale matters.
#
# lambda_pen = 0 reduces exactly to the unpenalised meta-association ranking,
# which is the anchor point of the grid.
lineage_score <- function(st, lambda_pen) {
  stopifnot(is.numeric(lambda_pen), length(lambda_pen) == 1L, is.finite(lambda_pen))
  h <- abs(st$z_meta) * st$consistency
  l <- st$lineage
  rscale <- function(v) {
    s <- stats::mad(v)
    if (!is.finite(s) || s <= 0) s <- stats::sd(v)
    if (!is.finite(s) || s <= 0) s <- 1
    s
  }
  l_std <- ((l - stats::median(l)) / rscale(l)) * rscale(h)
  v <- h - lambda_pen * l_std
  names(v) <- names(st$z_meta)
  v
}


# -----------------------------------------------------------------------------
# fit_preprocess(): LEARN the preprocessing transform from training rows only.
# -----------------------------------------------------------------------------
# Input : x            numeric matrix, samples in ROWS, probes in COLUMNS,
#                      with colnames set to probe IDs (e.g. "cg00000029").
#         max_features how many probes to keep after variance ranking.
#         max_missing  a probe is dropped if more than this fraction of training
#                      samples have a non-finite value for it.
#         feature_rank which unsupervised statistic ranks probes:
#                        "pooled"        - total variance across all training
#                                          rows. THE DEFAULT and the pre-V3
#                                          behaviour, bit-for-bit.
#                        "within_tissue" - pooled within-tissue variance, i.e.
#                                          variance after removing each tissue's
#                                          own mean. Requires `cancer`.
#                        "lineage_penalized" - V4. SUPERVISED: within-tissue
#                                          HRD meta-association minus a lineage
#                                          penalty. Requires `cancer` AND `y`.
#         cancer       tissue label per ROW of x. Required for "within_tissue"
#                      and "lineage_penalized"; ignored otherwise.
#         y            target per ROW of x. Required ONLY for
#                      "lineage_penalized"; ignored (and must be ignored)
#                      otherwise, so the two unsupervised paths stay
#                      bit-identical whether or not a caller passes it.
#         lambda_pen   lineage penalty weight for "lineage_penalized". Chosen in
#                      the inner folds by fit_en(), never by hand here.
#         block_size   column block width for the within-tissue sweep.
#
# Output: a list describing the transform. It contains NO sample data - only
#         per-probe summary statistics - so it is safe to save and ship.
#
# This function must only ever be handed TRAINING rows. Every caller below
# subsets x before calling it.
#
# STRUCTURE (introduced for V4, behaviour-preserving)
#   The body is split into two helpers:
#     preprocess_base()     missingness filter + training medians + imputation.
#                           This is the EXPENSIVE part (~73 s at full width) and
#                           it does not depend on the ranking statistic at all.
#     preprocess_finalize() ranking + top-max_features + centre/scale.
#   fit_preprocess() is exactly base-then-finalize, so every existing call is
#   bit-identical to before the split. fit_en()'s V4 path calls the two halves
#   separately so that a grid of lambda_pen values can share ONE base pass and
#   ONE meta-statistic pass over the fold's rows. That is a within-fold cache of
#   a quantity computed from that fold's own rows - it is NOT global
#   precomputation and nothing is ever reused across folds.
preprocess_base <- function(x, max_missing=0.05) {
  stopifnot(is.matrix(x), !is.null(colnames(x)), !anyDuplicated(colnames(x)))

  # ---------------------------------------------------------------------------
  # PERFORMANCE NOTE (for collaborators reading this for the first time)
  # ---------------------------------------------------------------------------
  # This function previously used apply(z, 2, median) and friends. apply() is
  # convenient but slow on a wide matrix: it walks the columns in interpreted R
  # code, allocating a fresh vector for every one of the ~336,000 columns.
  #
  # matrixStats does the identical arithmetic in compiled C with one pass over
  # memory. The NUMBERS DO NOT CHANGE - only the time taken. We verified this
  # explicitly: on input with no Inf values the two implementations agree to
  # floating-point exactness (max absolute difference ~1e-16).
  #
  # has_ms lets the file still run if matrixStats is not installed, falling back
  # to the original base-R calls. Every reduction below is written as
  #     if (has_ms) <fast version> else <base version>
  # so the two paths stay visibly side by side and can be compared by eye.
  # ---------------------------------------------------------------------------
  has_ms <- requireNamespace("matrixStats", quietly=TRUE)

  # ---------------------------------------------------------------------------
  # LEARNED QUANTITY 1 of 6: the missingness filter (which probes to keep).
  # ---------------------------------------------------------------------------
  # A probe is unusable if too many training samples lack a value for it.
  #
  # is.finite() is FALSE for NA, NaN, Inf and -Inf, so `nonfinite` is a
  # TRUE/FALSE matrix of "this cell is unusable". Taking the column MEAN of a
  # logical matrix gives the FRACTION of samples missing that probe, because R
  # treats TRUE as 1 and FALSE as 0.
  #
  # We store `nonfinite` in a variable instead of recomputing !is.finite(z)
  # later. At this scale that matrix is ~2.6 GB, so computing it twice would be
  # both slow and a real memory risk.
  nonfinite <- !is.finite(x)
  keep <- (if (has_ms) matrixStats::colMeans2(nonfinite) else colMeans(nonfinite)) <= max_missing
  if (!any(keep)) stop("No features pass train-only missingness")

  # Subset BOTH the data and the missingness mask so their columns stay aligned.
  # If these ever fell out of sync, we would impute using the wrong probe's
  # median - a silent, hard-to-detect corruption.
  z <- x[,keep,drop=FALSE]; nonfinite <- nonfinite[,keep,drop=FALSE]

  # ---------------------------------------------------------------------------
  # LEARNED QUANTITY 2 of 6: per-probe imputation medians (training rows only).
  # ---------------------------------------------------------------------------
  # Why median and not mean: beta values live in [0,1] and are usually BIMODAL -
  # a probe is typically either methylated (~0.9) or not (~0.1), with few
  # samples in between. The mean of such a distribution lands near 0.5, a value
  # almost no real sample has. The median lands on one of the two real modes.
  #
  # Why we set non-finite cells to NA first: median(na.rm=TRUE) removes NA and
  # NaN but NOT Inf. In the old code, a column containing Inf was counted as
  # missing by the filter above while still contributing Inf to its own median,
  # so that column's median could come back as Inf. Converting everything
  # unusable to NA up front makes the two steps agree on what "missing" means.
  any_missing <- any(nonfinite)
  if (any_missing) z[nonfinite] <- NA_real_
  med <- if (has_ms) matrixStats::colMedians(z,na.rm=TRUE) else apply(z,2,median,na.rm=TRUE)

  # matrixStats returns a bare unnamed vector, so we reattach the probe IDs.
  # Downstream code looks these up BY NAME (med[colnames(z)]), so losing the
  # names would silently mis-align medians to probes.
  names(med) <- colnames(z)

  # Fill the gaps. Two efficiencies worth understanding:
  #   1. We skip the whole block when nothing is missing.
  #   2. colAnys() tells us WHICH columns actually contain a gap. After the
  #      missingness filter most columns are complete, so looping over only the
  #      affected ones avoids touching hundreds of thousands of clean columns.
  # The loop body reads: "in column j, set every missing cell to that column's
  # median". Note the median came from the OTHER (observed) samples in the
  # training set - never from held-out data.
  if (any_missing) {
    cols <- which(if (has_ms) matrixStats::colAnys(nonfinite) else apply(nonfinite,2,any))
    for (j in cols) z[nonfinite[,j],j] <- med[j]
  }
  list(z=z, med=med, max_missing=max_missing, has_ms=has_ms)
}


# -----------------------------------------------------------------------------
# preprocess_finalize(): rank, select the top max_features, centre and scale.
# -----------------------------------------------------------------------------
# `base` is the output of preprocess_base() for THIS fold's training rows.
# `stats` is an optional pre-computed lineage_meta_stats() object for the SAME
# rows, supplied only so a lambda_pen grid does not repay for the sweep.
preprocess_finalize <- function(base, max_features=5000L,
                                feature_rank=c("pooled","within_tissue","lineage_penalized"),
                                cancer=NULL, y=NULL, lambda_pen=0,
                                block_size=20000L, stats=NULL) {
  feature_rank <- match.arg(feature_rank)
  z <- base$z; med <- base$med; has_ms <- base$has_ms

  # ---------------------------------------------------------------------------
  # LEARNED QUANTITY 3 of 6: unsupervised feature selection by variance.
  # ---------------------------------------------------------------------------
  # We keep the max_features probes that vary MOST across training samples.
  #
  # "Unsupervised" is the important word: this ranking never looks at y (the
  # HRD score). That matters because supervised selection - picking probes
  # because they correlate with the outcome - computed on the full dataset is
  # one of the classic ways to leak test information into a model and produce
  # impressive but fictional performance. Variance is a property of the
  # predictors alone, so it cannot leak the label.
  #
  # Probes with zero variance are dropped: they carry no information, and
  # dividing by their (zero) standard deviation later would produce NaN.
  #
  # V3-abs: when feature_rank="within_tissue" the ranking statistic is the
  # POOLED WITHIN-TISSUE variance instead of the total variance. Everything else
  # - the missingness filter, the medians, the centre/scale, the tie-break - is
  # untouched. z here is the TRAINING block only (the caller subset x), already
  # imputed with TRAINING medians, so this statistic is training-fold-only by
  # construction. A probe that is constant across the training block has zero
  # within-tissue variance too, so the `v > 0` filter still removes it; a probe
  # that varies only BETWEEN tissues also falls out here, which is the whole
  # point of the V3-abs change.
  #
  # V4-lineage-penalized: the ranking becomes SUPERVISED. See the long note on
  # lineage_meta_stats() above. z is still the TRAINING block only and y is
  # still the TRAINING target only (the caller subsets both with the same mask),
  # so the statistic is training-fold-only by construction exactly as V3's is -
  # but because it reads y, the leakage test must corrupt held-out LABELS too.
  if (feature_rank == "lineage_penalized") {
    # V4. The score is NOT a variance, so `v > 0` is the wrong eligibility test
    # here: a perfectly good probe can score zero or negative once the lineage
    # penalty is applied. Eligibility is instead "has non-zero total variance in
    # the training block", which is the same set the pooled path would admit,
    # and the ORDERING is by score. The tie-break on probe name is identical.
    st <- if (is.null(stats)) lineage_meta_stats(z, y, cancer, med=med, block_size=block_size)
          else stats
    v <- lineage_score(st, lambda_pen)
    names(v) <- colnames(z); ii <- which(is.finite(v) & st$nonconstant)
  } else {
    if (feature_rank == "within_tissue") {
      v <- within_tissue_var(z, cancer, med=med, block_size=block_size)
    } else {
      v <- if (has_ms) matrixStats::colVars(z) else apply(z,2,var)
    }
    names(v) <- colnames(z); ii <- which(is.finite(v) & v>0)
  }

  # Sort by decreasing variance, breaking ties on probe NAME. The tie-break
  # makes selection fully deterministic: two probes with identical variance
  # always resolve the same way regardless of the column order in the input
  # file, so the same data always yields the same frozen feature list.
  ii <- ii[order(-v[ii],names(v)[ii])]
  ii <- head(ii,max_features)
  if (length(ii)<2) stop("Fewer than two nonconstant features")

  # ---------------------------------------------------------------------------
  # LEARNED QUANTITIES 4 and 5 of 6: the centre and scale for standardisation.
  # ---------------------------------------------------------------------------
  # Elastic net penalises all coefficients equally, which is only fair if the
  # predictors share a common scale. We therefore subtract each probe's mean and
  # divide by its standard deviation - but using TRAINING statistics only, which
  # are stored here and reapplied unchanged to any future sample.
  #
  # This is why glmnet is later called with standardize=FALSE: the data arrives
  # already standardised, and letting glmnet redo it on a held-out batch would
  # recompute statistics from that batch, breaking the leakage boundary.
  z <- z[,ii,drop=FALSE]
  mu <- if (has_ms) matrixStats::colMeans2(z) else colMeans(z)
  s  <- if (has_ms) matrixStats::colSds(z)    else apply(z,2,sd)
  names(mu) <- colnames(z); names(s) <- colnames(z)

  # Return the transform. Note med is subset to the finally-selected features so
  # that every stored vector has the same length and the same order as
  # $features. apply_preprocess() relies on that correspondence.
  #
  # This object contains ONLY per-probe summary statistics - no sample-level
  # data - which is what makes it safe to save, share and ship inside a model.
  out <- list(features=colnames(z), median=med[colnames(z)],center=mu,scale=s,
              max_missing=base$max_missing, feature_rank=feature_rank)
  # lambda_pen is recorded ONLY on the V4 path. Appending it unconditionally
  # would change the shape of every "pooled" transform ever produced and break
  # the bit-identity guarantee the V3 tests assert.
  if (feature_rank == "lineage_penalized") out$lambda_pen <- lambda_pen
  out
}


# -----------------------------------------------------------------------------
# fit_preprocess(): the public entry point. base + finalize, unchanged contract.
# -----------------------------------------------------------------------------
fit_preprocess <- function(x, max_features=5000L, max_missing=0.05,
                           feature_rank=c("pooled","within_tissue","lineage_penalized"),
                           cancer=NULL, y=NULL, lambda_pen=0, block_size=20000L) {
  feature_rank <- match.arg(feature_rank)
  if (feature_rank %in% c("within_tissue","lineage_penalized")) {
    if (is.null(cancer)) stop(sprintf("feature_rank='%s' requires `cancer`", feature_rank))
    stopifnot(length(cancer) == nrow(x))
  }
  # V4 is the first SUPERVISED ranker in this file, so it is the first that can
  # be silently wrong by falling back to an unsupervised statistic. Refuse
  # loudly instead. A missing y here would make an entire sentinel meaningless
  # while looking perfectly healthy in the logs.
  if (feature_rank == "lineage_penalized") {
    if (is.null(y)) stop("feature_rank='lineage_penalized' requires `y`")
    stopifnot(length(y) == nrow(x))
    if (!all(is.finite(y))) stop("feature_rank='lineage_penalized' requires finite `y`")
    stopifnot(is.numeric(lambda_pen), length(lambda_pen) == 1L, is.finite(lambda_pen))
  }
  base <- preprocess_base(x, max_missing)
  preprocess_finalize(base, max_features=max_features, feature_rank=feature_rank,
                      cancer=cancer, y=y, lambda_pen=lambda_pen,
                      block_size=block_size)
}


# -----------------------------------------------------------------------------
# apply_preprocess(): APPLY a previously learned transform. Learns nothing.
# -----------------------------------------------------------------------------
# This is the held-out side of the leakage boundary. It reads pp$median,
# pp$center and pp$scale and never recomputes them.
#
# reject_missing controls whether excessive missingness is a hard error (used by
# the smoke test) or a soft flag (used everywhere in production, where
# predict_en() reports qc_fail instead of aborting the whole run for one bad
# sample).
apply_preprocess <- function(x,pp, reject_missing=TRUE) {
  if (anyDuplicated(colnames(x))) stop("Duplicate input features")

  # Rebuild the feature matrix BY NAME, not by position. This is what makes the
  # model invariant to column order in the incoming file, and it is explicitly
  # tested in tests/smoke_model.R by feeding in a reversed matrix.
  # Probes the model expects but the new data lacks stay NA here and are filled
  # with the stored training median below - which is the correct behaviour for a
  # single-sample assay that may not cover every probe.
  z <- matrix(NA_real_,nrow(x),length(pp$features),dimnames=list(rownames(x),pp$features))
  common <- intersect(colnames(x),pp$features);z[,common] <- x[,common,drop=FALSE]

  # Per-SAMPLE missingness (rowMeans, not colMeans): "how much of what this
  # model needs is absent for this particular sample?" This drives the
  # per-sample QC flag rather than any feature-level decision.
  miss <- rowMeans(!is.finite(z))
  if (reject_missing && any(miss>pp$max_missing)) stop("Too many missing model features")

  # Impute with the TRAINING medians. pp$median[j] is positionally aligned with
  # pp$features[j] because fit_preprocess() subset it by colnames(z).
  for (j in seq_len(ncol(z))) z[!is.finite(z[,j]),j] <- pp$median[j]

  # Centre and scale with the TRAINING statistics. Note this is the only
  # standardisation that ever happens - glmnet is called with standardize=FALSE.
  z <- sweep(sweep(z,2,pp$center,"-"),2,pp$scale,"/")

  # Carry the missing fraction along as an attribute so predict_en() can flag
  # samples without needing to recompute it.
  attr(z,"missing_fraction") <- miss;z
}


# -----------------------------------------------------------------------------
# inner_folds(): assign INNER cross-validation folds for hyperparameter tuning.
# -----------------------------------------------------------------------------
# Called from inside fit_en(), which is itself called once per OUTER
# leave-one-cancer-out fold. So this produces the inner layer of a properly
# nested CV: hyperparameters are chosen without ever touching the outer held-out
# cancer.
inner_folds <- function(patient,cancer,seed=260910L) {
  # One specimen per patient is a hard requirement. Two samples from the same
  # patient in different folds would leak that patient's methylation profile
  # across the fold boundary.
  if (anyDuplicated(patient)) stop("Use one specimen per patient for this MVP")

  # PRIMARY PATH: group by cancer type, so each remaining cancer becomes its own
  # inner fold. Combined with the outer loop this is leave-one-cancer-out nested
  # inside leave-one-cancer-out. Hyperparameters are therefore selected for their
  # ability to generalise ACROSS TISSUES, which is the actual scientific question,
  # rather than for within-tissue fit.
  if (length(unique(cancer))>=3) return(match(cancer,sort(unique(cancer))))

  # FALLBACK PATH: random 3-fold over patients. Unreachable from
  # train_baseline.R, which already refuses to run with fewer than 3 development
  # cancers; it exists for the smoke test and for future small-cohort use.
  if (length(patient)<9) stop("Need >=9 training patients for patient-fold fallback")
  set.seed(seed); as.integer(sample(rep(1:3,length.out=length(patient))))
}


# -----------------------------------------------------------------------------
# lambda_path(): build a data-derived regularisation path.
# -----------------------------------------------------------------------------
# WHY THIS EXISTS
#   lambda controls how hard the elastic net shrinks coefficients toward zero.
#   Its meaningful scale is a property of the DATA, not a universal constant:
#   it depends on the number of samples, the scale of the outcome, and how
#   strongly the predictors correlate with it.
#
#   This code previously used five hardcoded values (0.01, 0.1, 1, 10, 100).
#   That had two failure modes:
#
#     1. WASTED GRID POINTS. glmnet can compute lambda.max, the smallest lambda
#        at which EVERY coefficient is zero. Any lambda above it fits the same
#        intercept-only model, so those grid points are not merely uninformative
#        - they are duplicates of each other.
#
#     2. INVISIBLE BOUNDARY OPTIMA. If the best value was the smallest one
#        offered, that is the tuner saying "I wanted less regularisation than
#        you allowed". With a fixed grid and no check, that message was silently
#        discarded and the boundary was reported as though it were an optimum.
#
# HOW THE PATH IS BUILT
#   lambda.max = max|x'y| / (n * alpha) on standardised predictors - the
#   standard glmnet derivation. We then lay out n_lambda points log-spaced from
#   lambda.max down to lambda_min_ratio * lambda.max. Log spacing is used
#   because lambda acts multiplicatively: the step from 0.01 to 0.1 matters as
#   much as the step from 1 to 10.
#
# LEAKAGE
#   x and y here are TRAINING rows for the current fold. The path is therefore a
#   learned quantity like any other, computed inside the same boundary. It never
#   sees the held-out cancer.
#
# Dividing by alpha means a small alpha (near ridge) yields a larger lambda.max,
# which is correct: ridge penalties need larger lambda for equivalent shrinkage.
# alpha is floored at 0.01 to avoid dividing by zero if someone passes alpha=0.
lambda_path <- function(x,y,alpha,n_lambda=8L,lambda_min_ratio=0.001) {
  n <- length(y)
  yc <- y-mean(y)
  # crossprod gives x'yc; the largest absolute value sets the point at which the
  # first coefficient would become non-zero.
  lam_max <- max(abs(crossprod(x,yc)))/(n*max(alpha,0.01))
  # Degenerate guard: if x carries no signal at all, fall back to a token path
  # rather than returning NaN and taking the whole run down.
  if(!is.finite(lam_max)||lam_max<=0) return(c(1,0.1,0.01))
  exp(seq(log(lam_max),log(lam_max*lambda_min_ratio),length.out=n_lambda))
}


# -----------------------------------------------------------------------------
# check_lambda_boundary(): warn when the tuner picks an endpoint.
# -----------------------------------------------------------------------------
# A selected lambda at the TOP of the path means the model wanted even more
# shrinkage - often a sign there is little signal. At the BOTTOM it wanted less
# regularisation, which risks overfitting and suggests the path was too narrow.
#
# Either way the search has hit a wall we imposed rather than finding an
# interior optimum, and that fact belongs in the log and the audit trail. This
# warns rather than stops: a boundary optimum is informative, not fatal.
check_lambda_boundary <- function(selected_lambda,path,alpha,label="") {
  if(!length(path)||!is.finite(selected_lambda)) return(invisible(NULL))
  at_top <- isTRUE(all.equal(selected_lambda,max(path)))
  at_bot <- isTRUE(all.equal(selected_lambda,min(path)))
  if(at_top||at_bot) {
    warning(sprintf(
      "%sSelected lambda %.5g (alpha %.2g) is at the %s of its path [%.5g, %.5g]. %s",
      if(nzchar(label)) paste0(label,": ") else "",
      selected_lambda, alpha, if(at_top) "TOP" else "BOTTOM",
      min(path), max(path),
      if(at_top) "The tuner wanted MORE shrinkage; signal may be weak."
      else "The tuner wanted LESS regularisation; consider a wider path."),
      call.=FALSE)
  }
  invisible(list(at_top=at_top,at_bottom=at_bot))
}


# -----------------------------------------------------------------------------
# fit_en(): tune and fit the elastic net. The heart of the modelling code.
# -----------------------------------------------------------------------------
# Two phases:
#   Phase 1 - inner CV over an (alpha, lambda) grid to pick hyperparameters.
#   Phase 2 - refit preprocessing and model on ALL supplied rows using the
#             winning hyperparameters.
#
# Everything passed in here is training data for the current outer fold. The
# outer held-out cancer is never visible to this function.
fit_en <- function(x,y,patient,cancer,max_features=5000L,seed=260910L,
                   feature_rank=c("pooled","within_tissue","lineage_penalized"),
                   block_size=20000L,lambda_pen_grid=c(0,0.25,0.5,1,2)) {
  feature_rank <- match.arg(feature_rank)
  if (!requireNamespace("glmnet",quietly=TRUE)) stop("Install glmnet via scripts/setup.R")
  stopifnot(length(y)==nrow(x),all(is.finite(y)),length(patient)==length(y))
  # A constant target makes R2 undefined and the fit meaningless.
  if (sd(y)==0) stop("Constant training target")
  is_v4 <- feature_rank == "lineage_penalized"
  if (is_v4) {
    stopifnot(is.numeric(lambda_pen_grid), length(lambda_pen_grid) >= 1L,
              all(is.finite(lambda_pen_grid)), !anyDuplicated(lambda_pen_grid))
    lambda_pen_grid <- sort(lambda_pen_grid)
  }
  # On the two unsupervised paths lambda_pen is not a hyperparameter at all, so
  # the grid collapses to a single inert value and the tuning table keeps its
  # pre-V4 shape exactly.
  pens <- if (is_v4) lambda_pen_grid else NA_real_

  folds <- inner_folds(patient,cancer,seed)

  # ---------------------------------------------------------------------------
  # THE HYPERPARAMETER GRID
  # ---------------------------------------------------------------------------
  # alpha = the elastic-net mixing parameter.
  #   1.0 is pure lasso: aggressive, drives most coefficients to exactly zero,
  #       and when two probes are correlated it tends to keep one arbitrarily.
  #   0.1 is nearly ridge: keeps many small coefficients and shares weight
  #       across correlated probes.
  #   Methylation probes are correlated in blocks (neighbouring CpGs in the same
  #   island move together), so the low-alpha end genuinely matters here and is
  #   not just filler.
  #
  # lambda = regularisation strength. See the long note in lambda_path() below
  #   for why this is now DERIVED FROM THE DATA rather than hardcoded.
  # ---------------------------------------------------------------------------
  alphas <- c(0.1,0.5,1)

  # Build one lambda path per alpha. The path is computed from the TRAINING rows
  # of this fold only, using the same preprocessing the model will see, so it
  # respects the leakage boundary exactly as every other learned quantity does.
  #
  # Why per alpha: lambda.max depends on alpha (it is max|x'y|/(n*alpha)), so a
  # single shared path would be correctly scaled for at most one of them.
  #
  # EFFICIENCY: this preprocessing fit is deliberately hoisted out of Phase 2
  # below and reused there. fit_preprocess() is the most expensive operation in
  # the file (~73 s at full width), so computing it once and using it for both
  # the lambda path and the final refit costs nothing extra - the original code
  # already paid for one call here.
  # V3-abs: `cancer` here is the OUTER-TRAINING tissue vector. The held-out
  # cancer is absent from x and from cancer by construction (the caller subsets
  # both with the same mask), so the within-tissue ranking cannot see it.
  #
  # V4: lambda_pen is a HYPERPARAMETER, so the outer-training block yields one
  # transform PER candidate lambda_pen. They all share a single preprocess_base()
  # and a single lineage_meta_stats() sweep over these same outer-training rows,
  # because neither depends on lambda_pen - only the final ranking does. That is
  # the whole reason a 5-point penalty grid costs barely more than one point.
  # `y` here is the OUTER-TRAINING target, subset with the same mask as x.
  base_out <- preprocess_base(x)
  stats_out <- if (is_v4) lineage_meta_stats(base_out$z,y,cancer,med=base_out$med,
                                             block_size=block_size) else NULL
  pen_key <- function(p) sprintf("%.10g",p)
  pp_by_pen <- lapply(pens, function(p)
    preprocess_finalize(base_out,max_features,feature_rank=feature_rank,
                        cancer=cancer,y=y,
                        lambda_pen=if (is_v4) p else 0,
                        block_size=block_size,stats=stats_out))
  names(pp_by_pen) <- pen_key(pens)
  rm(base_out); invisible(gc(verbose=FALSE))
  z_by_pen <- lapply(pp_by_pen, function(q) apply_preprocess(x,q,FALSE))

  # One lambda path per (lambda_pen, alpha). The path is derived from the
  # standardised training matrix, which differs between penalties because the
  # selected probes differ, so a single shared path would be mis-scaled for all
  # but one penalty.
  lam_by_pen <- lapply(pen_key(pens), function(pk) {
    l <- lapply(alphas, function(a) lambda_path(z_by_pen[[pk]],y,a))
    names(l) <- as.character(alphas); l
  })
  names(lam_by_pen) <- pen_key(pens)
  # Back-compatible name for the single-penalty (pooled / within_tissue) paths.
  lam_by_alpha <- lam_by_pen[[pen_key(pens[1])]]

  grid <- do.call(rbind,lapply(pens,function(p) do.call(rbind,lapply(alphas,function(a) {
    d <- data.frame(alpha=a,lambda=lam_by_pen[[pen_key(p)]][[as.character(a)]])
    # The lambda_pen column exists only on the V4 path, so the tuning table of a
    # pooled or within_tissue fit keeps exactly its pre-V4 shape.
    if (is_v4) d$lambda_pen <- p
    d
  }))))
  losses <- matrix(NA_real_,nrow(grid),length(unique(folds)))

  # ---- Phase 1: inner cross-validation -------------------------------------
  for (f in sort(unique(folds))) {
    tr <- folds!=f;va <- !tr

    # CRITICAL: preprocessing is re-learned from scratch on THIS inner fold's
    # training rows. It is not hoisted out of the loop. Hoisting it would let
    # every inner validation fold see statistics computed from itself.
    #
    # (The hoisted `pp` above is a DIFFERENT object: it is fitted on all rows
    # passed to fit_en, which are all training rows for this outer fold, and is
    # used only for the lambda path and the Phase 2 refit - never to evaluate an
    # inner validation fold.)
    # V3-abs: the ranking statistic is recomputed here too, from x[tr,] and
    # cancer[tr] only. There is no hoisted or global ranking - each inner fold
    # ranks on its own rows, exactly as it re-learns medians and centre/scale.
    # V4: the meta-association sweep is likewise recomputed HERE, from x[tr,],
    # y[tr] and cancer[tr] only. stats_in is local to this inner fold and is
    # never reused by another fold; it is shared only across the lambda_pen grid
    # WITHIN this fold, which is legitimate because every penalty is ranking the
    # same inner-training rows.
    base_in <- preprocess_base(x[tr,,drop=FALSE])
    stats_in <- if (is_v4) lineage_meta_stats(base_in$z,y[tr],cancer[tr],
                                              med=base_in$med,block_size=block_size) else NULL

    for (p in pens) {
      pk <- pen_key(p)
      pp_in <- preprocess_finalize(base_in,max_features,feature_rank=feature_rank,
                                   cancer=cancer[tr],y=y[tr],
                                   lambda_pen=if (is_v4) p else 0,
                                   block_size=block_size,stats=stats_in)
      ztr <- apply_preprocess(x[tr,,drop=FALSE],pp_in,FALSE)
      zva <- apply_preprocess(x[va,,drop=FALSE],pp_in,FALSE)

      for (a in alphas) {
        # Each alpha gets ITS OWN lambda path, sorted decreasing because that is
        # glmnet's expected convention and lets it use warm starts along the path.
        # Fitting the whole path in one call is far cheaper than one call per
        # lambda.
        # standardize=FALSE: see the header note. ztr is already standardised
        # using training-only statistics and glmnet must not redo it.
        lam_a <- sort(lam_by_pen[[pk]][[as.character(a)]],decreasing=TRUE)
        mod <- glmnet::glmnet(ztr,y[tr],alpha=a,lambda=lam_a,standardize=FALSE)
        ids <- if (is_v4) which(grid$alpha==a & grid$lambda_pen==p) else which(grid$alpha==a)
        for (g in ids) {
          # s= selects one lambda from the fitted path.
          pred <- as.numeric(predict(mod,zva,s=grid$lambda[g]))
          # Equal weight to each inner-validation cancer, independent of size.
          # NOTE: on the primary path each inner fold IS a single cancer, so this
          # tapply collapses to a plain mean within the fold; the macro-averaging
          # actually comes from rowMeans(losses) below. The tapply only does real
          # work on the patient-fold fallback path.
          losses[g,match(f,sort(unique(folds)))] <- mean(tapply(abs(pred-y[va]),cancer[va],mean))
        }
      }
      rm(ztr,zva,pp_in)
    }
    rm(base_in,stats_in); invisible(gc(verbose=FALSE))
  }

  # Average each grid point's loss across inner folds and take the winner.
  # Because each fold is a cancer, this is a macro average over tissues: a
  # hyperparameter setting that is excellent on BRCA and terrible everywhere
  # else will lose to one that is uniformly decent.
  #
  # V4: lambda_pen is selected HERE, by exactly this macro-averaged inner loss,
  # jointly with alpha and lambda. It is never chosen by hand and never chosen
  # by looking at the held-out tissue - the held-out tissue is not in x at all.
  grid$inner_macro_mae <- rowMeans(losses)
  best <- which.min(grid$inner_macro_mae)
  best_alpha <- grid$alpha[best]; best_lambda <- grid$lambda[best]
  best_pen <- if (is_v4) grid$lambda_pen[best] else NA_real_
  best_pk <- pen_key(if (is_v4) best_pen else pens[1])

  # Was the winner at an endpoint of its own path? If so the tuner hit a wall we
  # imposed rather than finding an interior optimum. Warn, and record the fact in
  # the bundle so it survives into the audit trail rather than living only in a
  # console log that may be lost.
  best_path <- lam_by_pen[[best_pk]][[as.character(best_alpha)]]
  boundary <- check_lambda_boundary(best_lambda,best_path,best_alpha,
                                    label=sprintf("fit_en (n=%d)",length(y)))

  # ---- Phase 2: refit on all supplied training rows -------------------------
  # `pp` and `z` were fitted above on the FULL x - every row passed to fit_en,
  # all of which are training rows for this outer fold. The held-out cancer is
  # still excluded, so reusing them here is correct AND saves one full
  # fit_preprocess call (the single most expensive operation in the file).
  # V4: we pick the transform belonging to the winning lambda_pen. It too was
  # fitted on these same outer-training rows.
  pp <- pp_by_pen[[best_pk]]; z <- z_by_pen[[best_pk]]
  rm(z_by_pen); invisible(gc(verbose=FALSE))
  mod <- glmnet::glmnet(z,y,alpha=best_alpha,
                        lambda=sort(best_path,decreasing=TRUE),standardize=FALSE)

  # LEARNED QUANTITY 6: a crude out-of-distribution cutoff.
  # dist is the root-mean-square standardised distance of each training sample
  # from the training centroid. Samples far outside the training cloud get
  # flagged at prediction time.
  # Heuristic score only: not a calibrated domain-shift test.
  dist <- sqrt(rowMeans(z^2));ood_cut <- as.numeric(quantile(dist,0.99,names=FALSE))

  # The returned bundle is self-contained: it carries the transform, the model,
  # the winning hyperparameters, the full tuning table (so the grid search is
  # auditable after the fact), the fold assignments, the OOD cutoff, the seed,
  # the training mean that serves as the null comparator downstream, and the
  # lambda paths plus boundary status so the tuning can be reviewed later.
  out <- list(preprocess=pp,model=mod,alpha=best_alpha,lambda=best_lambda,
       tuning=grid,inner_folds=data.frame(patient_id=patient,cancer_type=cancer,fold=folds),
       lambda_paths=lam_by_alpha,
       lambda_at_boundary=if(is.null(boundary)) NA else (boundary$at_top||boundary$at_bottom),
       lambda_boundary_side=if(is.null(boundary)) NA_character_
                            else if(boundary$at_top) "top" else if(boundary$at_bottom) "bottom" else "interior",
       ood_cut=ood_cut,seed=seed,training_n=length(y),training_mean=mean(y),
       feature_rank=feature_rank,target="reference_HRDsum")
  # Again, V4-only fields. Adding them unconditionally would change the shape of
  # every bundle this project has ever produced.
  if (is_v4) {
    out$lambda_pen <- best_pen
    out$lambda_pen_grid <- lambda_pen_grid
    out$lambda_pen_paths <- lam_by_pen
    out$lambda_pen_at_boundary <- best_pen %in% range(lambda_pen_grid)
  }
  out
}


# -----------------------------------------------------------------------------
# predict_en(): apply a fitted bundle to new samples. Learns nothing.
# -----------------------------------------------------------------------------
# Returns one row per input sample with the prediction plus QC flags. Note it
# returns the raw prediction even when QC fails; it is the CALLER's job to
# decide whether to display it (predict_frozen.R blanks the display column but
# deliberately retains the raw value for analysis).
predict_en <- function(bundle,x) {
  # FALSE here means "flag, do not abort" - one bad sample in a batch should not
  # kill predictions for the rest.
  z <- apply_preprocess(x,bundle$preprocess,FALSE)
  raw <- as.numeric(predict(bundle$model,z,s=bundle$lambda))

  # Same distance metric as used to derive ood_cut at fit time, applied with the
  # stored cutoff rather than a freshly computed one.
  score <- sqrt(rowMeans(z^2));missing <- attr(z,"missing_fraction")
  fail <- missing>bundle$preprocess$max_missing
  ood <- score>bundle$ood_cut

  # reportable is the conjunction: a sample must both pass QC and sit inside the
  # training distribution before its prediction should be shown to anyone.
  data.frame(sample_id=rownames(x),predicted_reference_HRDsum=raw,
             missing_fraction=missing,ood_score=score,ood=ood,qc_fail=fail,
             reportable=!(fail|ood),stringsAsFactors=FALSE)
}


# -----------------------------------------------------------------------------
# conformal_q(): split-conformal prediction interval half-width.
# -----------------------------------------------------------------------------
# Given labels and predictions on a CALIBRATION set that the model was NOT fit
# on, return the residual quantile q such that [pred - q, pred + q] has the
# requested marginal coverage.
#
# k = ceiling((n+1) * coverage) is the textbook finite-sample-valid form (the
# +1 accounts for the new test point). If k > n there are too few calibration
# points to certify the requested coverage, and returning Inf is the honest
# answer: "no usable interval" rather than a falsely narrow one.
#
# IMPORTANT: this guarantee is marginal and assumes the test point is
# exchangeable with the calibration set. Applying a TCGA-calibrated interval to
# pediatric samples breaks exchangeability, so coverage there is not guaranteed.
conformal_q <- function(y,pred,coverage=0.95) {
  good <- is.finite(y)&is.finite(pred);res <- sort(abs(y[good]-pred[good]));n <- length(res)
  k <- ceiling((n+1)*coverage)
  if (!n || k>n) return(Inf)
  res[k]
}


# -----------------------------------------------------------------------------
# metrics(): evaluation panel for one set of predictions.
# -----------------------------------------------------------------------------
# Non-finite pairs are dropped, but n is reported so any such drop is visible.
metrics <- function(y,p) {
  ok <- is.finite(y)&is.finite(p);y<-y[ok];p<-p[ok];n<-length(y)
  if (!n) return(data.frame(n=0,MAE=NA,RMSE=NA,R2=NA,Pearson=NA,Spearman=NA,calibration_intercept=NA,calibration_slope=NA,bias=NA))

  # Correlations and the calibration regression need genuine variation in both
  # vectors; guard so a constant predictor yields NA rather than an error.
  variable <- n>2 && sd(y)>0 && sd(p)>0

  # Calibration: regress OBSERVED on PREDICTED. A perfectly calibrated model has
  # intercept 0 and slope 1. Slope < 1 means predictions are over-dispersed
  # (too extreme); slope > 1 means they are shrunk toward the mean.
  cal <- if (variable) coef(lm(y~p)) else c(NA,NA)

  data.frame(n=n,MAE=mean(abs(y-p)),RMSE=sqrt(mean((y-p)^2)),
             # R2 IS COMPUTED AGAINST THE MEAN OF y AS PASSED IN. Callers pass
             # one held-out cancer at a time, so mean(y) is that cancer's OWN
             # mean. This is the strict, conservative baseline: scoring above
             # zero requires beating an oracle that already knows the held-out
             # cancer's average HRDsum. It is NOT the training mean, and it is
             # therefore a harder bar than the null_MAE reported by
             # train_baseline.R - be careful not to compare the two directly.
             R2=if(n>1 && sum((y-mean(y))^2)>0) 1-sum((y-p)^2)/sum((y-mean(y))^2) else NA,
             # Pearson measures linear agreement; Spearman measures rank
             # agreement and is the more robust of the two if the relationship
             # turns out to be monotone but not linear.
             Pearson=if(variable) cor(y,p) else NA,Spearman=if(variable) cor(y,p,method="spearman") else NA,
             # bias is signed mean error: positive means systematic over-prediction.
             calibration_intercept=unname(cal[1]),calibration_slope=unname(cal[2]),bias=mean(p-y))
}


# =============================================================================
# TISSUE AND PERMUTATION CONTROLS
# =============================================================================
# Methylation is among the strongest tissue-of-origin signals in genomics, and
# HRDsum varies substantially BY cancer type. A model handed ~336k CpGs and a
# tissue-correlated outcome can therefore score well by learning lineage and
# nothing about HRD biology. Leave-one-cancer-out does not prevent this: it
# stops the model memorising a held-out type's mean, but not from reasoning
# "this looks squamous, squamous tumours score around X".
#
# The controls below separate two distinct estimands (see docs/22):
#
#   E1 between-tissue : can the model rank cancer TYPES by typical scar burden?
#   E2 within-tissue  : given two tumours OF THE SAME TYPE, can it tell which
#                       has the higher burden?
#
# E2 is the clinically meaningful question and the one an N-of-1 pediatric
# application actually requires. A pooled metric silently averages the two.
#
# SCOPE WARNING. Within a single LOCO fold every sample shares one cancer type,
# so tissue_mean_null() collapses to that fold's own mean and within-tissue
# centering collapses to plain centering. Per-fold output from these functions
# is still meaningful, but the between-tissue confound only becomes visible
# when predictions from ALL folds are POOLED and then compared against the
# tissue-mean null. Run them both ways.
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# tissue_mean_null(): the decisive null comparator.
# -----------------------------------------------------------------------------
# Returns, for each sample, the mean outcome of its OWN cancer type. This is an
# oracle - it uses labels unavailable at deployment - which is exactly the point.
# If the model cannot beat an opponent that knows nothing except the tissue's
# average, then the model has learned tissue lineage, not HRD biology.
#
# Note this is the same baseline the R2 in metrics() already uses per fold; the
# function makes it explicit, available on MAE, and usable on pooled data.
tissue_mean_null <- function(y,cancer) {
  stopifnot(length(y)==length(cancer))
  cancer <- as.character(cancer)
  mu <- tapply(y,cancer,function(v) mean(v[is.finite(v)]))
  as.numeric(mu[cancer])
}


# -----------------------------------------------------------------------------
# permute_within_tissue(): shuffle the outcome WITHIN each cancer type.
# -----------------------------------------------------------------------------
# Permuting globally would also destroy tissue structure, producing a null so
# weak that any lineage-aware model beats it trivially. Permuting within tissue
# holds the between-tissue signal fixed and asks only whether the within-tissue
# ordering carries information - the correct null for E2.
#
# The length(i) > 1L guard is deliberate: sample() on a length-1 vector would
# reinterpret it as sample(1:x, ...) and return an arbitrary index.
permute_within_tissue <- function(y,cancer,seed=NULL) {
  if(!is.null(seed)) set.seed(seed)
  cancer <- as.character(cancer);out <- y
  for(g in unique(cancer)) {
    i <- which(cancer==g)
    if(length(i)>1L) out[i] <- y[i][sample.int(length(i))]
  }
  out
}


# -----------------------------------------------------------------------------
# permutation_test_within_tissue(): is within-tissue agreement above chance?
# -----------------------------------------------------------------------------
# Holds the fitted predictions fixed and permutes the observed outcome within
# tissue n_perm times, building a null distribution of MAE.
#
# WHAT THIS DOES NOT TEST. Because predictions are not recomputed, this cannot
# detect a pipeline that manufactures signal from noise (that requires refitting
# on permuted labels, which costs a full LOCO run per permutation). Treat a
# small p-value as "these predictions carry within-tissue information", not as
# "the pipeline is leakage-free".
#
# p-value uses the standard (1 + #{as extreme}) / (n_perm + 1) form, which is
# never zero and stays valid for small n_perm.
permutation_test_within_tissue <- function(y,pred,cancer,n_perm=1000L,seed=260910L) {
  ok <- is.finite(y)&is.finite(pred)
  y<-y[ok];pred<-pred[ok];cancer<-as.character(cancer)[ok];n<-length(y)
  empty <- data.frame(n=n,observed_MAE=NA_real_,null_MAE_median=NA_real_,
                      null_MAE_q05=NA_real_,perm_p_value=NA_real_,n_perm=0L)
  if(n<3L) return(empty)

  # A permutation is only informative where some tissue has >1 sample.
  if(!any(table(cancer)>1L)) return(empty)

  obs <- mean(abs(y-pred));set.seed(seed)
  null_mae <- vapply(seq_len(n_perm),
                     function(k) mean(abs(permute_within_tissue(y,cancer)-pred)),
                     numeric(1))
  # Lower MAE is better, so "at least as extreme" means at least as SMALL.
  data.frame(n=n,observed_MAE=obs,null_MAE_median=median(null_mae),
             null_MAE_q05=as.numeric(quantile(null_mae,0.05,names=FALSE)),
             perm_p_value=(1+sum(null_mae<=obs))/(n_perm+1),
             n_perm=as.integer(n_perm))
}


# -----------------------------------------------------------------------------
# null_panel(): model against both nulls, side by side.
# -----------------------------------------------------------------------------
# skill = 1 - MAE_model / MAE_null. Positive means the model beats that null;
# zero means it matches it; negative means it is worse than predicting the mean.
#
# skill_vs_tissue_mean is the number to lead with. skill_vs_training_mean is the
# weaker, more flattering comparison and should never be quoted alone.
null_panel <- function(y,pred,cancer,training_mean=NA_real_) {
  ok <- is.finite(y)&is.finite(pred)
  y<-y[ok];pred<-pred[ok];cancer<-as.character(cancer)[ok];n<-length(y)
  if(!n) return(data.frame(n=0L,MAE_model=NA_real_,MAE_tissue_mean_null=NA_real_,
                           MAE_training_mean_null=NA_real_,
                           skill_vs_tissue_mean=NA_real_,skill_vs_training_mean=NA_real_))
  mae_model  <- mean(abs(y-pred))
  mae_tissue <- mean(abs(y-tissue_mean_null(y,cancer)))
  mae_train  <- if(is.finite(training_mean)) mean(abs(y-training_mean)) else NA_real_
  data.frame(n=n,MAE_model=mae_model,MAE_tissue_mean_null=mae_tissue,
             MAE_training_mean_null=mae_train,
             skill_vs_tissue_mean=if(mae_tissue>0) 1-mae_model/mae_tissue else NA_real_,
             skill_vs_training_mean=if(is.finite(mae_train)&&mae_train>0) 1-mae_model/mae_train else NA_real_)
}


# -----------------------------------------------------------------------------
# within_tissue_metrics(): the E2 estimand in isolation.
# -----------------------------------------------------------------------------
# Centers BOTH the outcome and the prediction by cancer type, removing each
# tissue's offset, then measures agreement on what remains. Centering the
# prediction as well as the outcome is what makes this a within-tissue question:
# a model whose entire skill is tissue-level offsets scores ~0 here.
#
# Spearman is the more robust summary if the relationship is monotone but not
# linear, which is plausible for a bounded scar score.
within_tissue_metrics <- function(y,pred,cancer) {
  ok <- is.finite(y)&is.finite(pred)
  y<-y[ok];pred<-pred[ok];cancer<-as.character(cancer)[ok];n<-length(y)
  if(n<3L) return(data.frame(n=n,within_Pearson=NA_real_,within_Spearman=NA_real_,within_MAE=NA_real_))
  yc <- y-tissue_mean_null(y,cancer);pc <- pred-tissue_mean_null(pred,cancer)
  v <- sd(yc)>0 && sd(pc)>0
  data.frame(n=n,
             within_Pearson=if(v) cor(yc,pc) else NA_real_,
             within_Spearman=if(v) cor(yc,pc,method="spearman") else NA_real_,
             within_MAE=mean(abs(yc-pc)))
}


# =============================================================================
# PURITY CONFOUNDING
# =============================================================================
# Tumour purity confounds this design through THREE paths at once, which is why
# it needs its own treatment rather than being folded into the tissue controls:
#
#   1. Purity is directly readable from methylation. An observed beta is a
#      mixture, beta_obs ~= pi*beta_tumour + (1-pi)*beta_normal, so the feature
#      matrix carries purity information whether or not anyone wants it to.
#   2. Purity shapes the LABEL. Reference HRDsum comes from SNP6 + ABSOLUTE
#      segmentation, and low-purity samples yield attenuated, noisier scar
#      calls.
#   3. Purity plausibly correlates with scar burden biologically.
#
# Paths 1 and 2 together are sufficient for a model to score well by inferring
# purity from methylation and exploiting purity's correlation with the label.
# That is a real statistical association and a scientifically empty one.
#
# train_baseline.R deliberately excludes purity from the FEATURES. That does not
# adjust for it - it only makes it unmeasured. These functions measure it.
#
# WHY NOT JUST REGRESS PURITY OUT. Purity influences both the features and the
# label, so naive residualisation can induce bias rather than remove it (see
# docs/22). Stratification and matching are the primary tools here;
# residualisation is a sensitivity analysis at most.
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# purity_confound_legs(): does the confounding pathway actually exist?
# -----------------------------------------------------------------------------
# Both legs must hold for purity confounding to operate, so test both before
# worrying about it:
#
#   leg 1  predicted HRDsum ~ purity   (does the model's output track purity?)
#   leg 2  observed  HRDsum ~ purity   (does the label track purity?)
#
# Correlations are computed WITHIN cancer type and macro-averaged, because a
# pooled correlation here would itself be confounded by tissue.
purity_confound_legs <- function(y,pred,purity,cancer) {
  ok <- is.finite(y)&is.finite(pred)&is.finite(purity)
  y<-y[ok];pred<-pred[ok];purity<-purity[ok];cancer<-as.character(cancer)[ok];n<-length(y)
  empty <- data.frame(n=n,n_types=0L,
                      cor_pred_purity_within=NA_real_,cor_label_purity_within=NA_real_,
                      cor_pred_purity_pooled=NA_real_,cor_label_purity_pooled=NA_real_)
  if(n<3L) return(empty)

  # Per-type Spearman, then macro-average over types with usable variation.
  per_type <- function(a,b) {
    vals <- vapply(split(seq_along(a),cancer),function(i) {
      if(length(i)<3L||sd(a[i])==0||sd(b[i])==0) return(NA_real_)
      suppressWarnings(cor(a[i],b[i],method="spearman"))
    },numeric(1))
    vals[is.finite(vals)]
  }
  cpp <- per_type(pred,purity); clp <- per_type(y,purity)
  pooled_cor <- function(a,b) if(sd(a)>0&&sd(b)>0) suppressWarnings(cor(a,b,method="spearman")) else NA_real_

  data.frame(n=n,n_types=length(unique(cancer)),
             cor_pred_purity_within=if(length(cpp)) mean(cpp) else NA_real_,
             cor_label_purity_within=if(length(clp)) mean(clp) else NA_real_,
             cor_pred_purity_pooled=pooled_cor(pred,purity),
             cor_label_purity_pooled=pooled_cor(y,purity))
}


# -----------------------------------------------------------------------------
# purity_stratified_metrics(): is the model only a purity detector?
# -----------------------------------------------------------------------------
# Splits samples into purity tertiles (or n_strata quantile bins) and reports
# the full null panel plus within-tissue agreement inside each one.
#
# The pattern to watch for is performance that holds in the high-purity stratum
# and collapses in the low-purity stratum. That is the signature of a model
# reading tumour content rather than scar burden. Roughly equal skill across
# strata is the reassuring result.
#
# Strata are cut on quantiles of the observed purity distribution, so bins are
# balanced by construction rather than by arbitrary cutpoints.
purity_stratified_metrics <- function(y,pred,purity,cancer,training_mean=NA_real_,n_strata=3L) {
  ok <- is.finite(y)&is.finite(pred)&is.finite(purity)
  y<-y[ok];pred<-pred[ok];purity<-purity[ok];cancer<-as.character(cancer)[ok];n<-length(y)
  if(n<3L*n_strata) return(NULL)

  qs <- unique(quantile(purity,probs=seq(0,1,length.out=n_strata+1L),names=FALSE))
  if(length(qs)<3L) return(NULL)   # too little variation to stratify
  bin <- cut(purity,breaks=qs,include.lowest=TRUE,labels=FALSE)

  do.call(rbind,lapply(sort(unique(bin)),function(b) {
    i <- which(bin==b)
    np <- null_panel(y[i],pred[i],cancer[i],training_mean)
    wt <- within_tissue_metrics(y[i],pred[i],cancer[i])
    data.frame(purity_stratum=b,
               purity_min=min(purity[i]),purity_median=median(purity[i]),purity_max=max(purity[i]),
               np,within_Pearson=wt$within_Pearson,within_Spearman=wt$within_Spearman,
               row.names=NULL)
  }))
}


# -----------------------------------------------------------------------------
# purity_matched_subset(): sensitivity analysis on a narrow purity band.
# -----------------------------------------------------------------------------
# Restricting to a narrow band removes most purity variation, so if the signal
# survives it is unlikely to be purity-driven. Power drops - that is expected
# and is why this is a sensitivity analysis, not the primary result. Judge the
# DIRECTION of the estimate, not its significance.
purity_matched_subset <- function(y,pred,purity,cancer,lower=0.5,upper=0.8,training_mean=NA_real_) {
  i <- which(is.finite(purity)&purity>=lower&purity<=upper&is.finite(y)&is.finite(pred))
  if(length(i)<10L) return(NULL)
  cancer <- as.character(cancer)
  cbind(data.frame(purity_band=paste0("[",lower,",",upper,"]"),n_retained=length(i)),
        null_panel(y[i],pred[i],cancer[i],training_mean),
        within_tissue_metrics(y[i],pred[i],cancer[i])[,c("within_Pearson","within_Spearman")])
}


# =============================================================================
# THE LOCK: one source of truth for which samples are locked
# =============================================================================
# The locked CNS partition is spelled out in SIX places (train_baseline.R,
# loco_one_fold.R, c1_rank_one_fold.R, loco_merge.R, R/rank_model.R and
# config/analysis_protocol.json) as the hardcoded list c("GBM","LGG"), while
# data/processed/master_samples.tsv carries an independent `partition` column
# written by scripts/build_master.py. Two sources of truth that are never
# compared will eventually disagree, and the failure is silent in the worst
# possible direction: a locked sample entering development, or a development
# sample being withheld from a fold that claims to have used it.
#
# CNS_LOCKED_TYPES is the canonical constant (R/rank_model.R's CNS_LOCKED is now
# an alias of it), and assert_partition_matches_cns() is the executable
# comparison. It does NOT change which samples are locked - it refuses to
# proceed when the two definitions disagree.
CNS_LOCKED_TYPES <- c("GBM", "LGG")
CNS_LOCKED_LABEL <- "locked_CNS"
DEVELOPMENT_LABEL <- "development"

# -----------------------------------------------------------------------------
# assert_partition_matches_cns(): the hardcoded list and the column must AGREE.
# -----------------------------------------------------------------------------
# Errors unless the agreement is exact and two-way:
#   * every row whose cancer_type is in `locked_types` has partition == locked
#   * every row whose partition is locked has a cancer_type in `locked_types`
# Also refuses NA/blank partitions and labels outside {locked, development},
# both of which make the column unusable as an authority.
#
# `require_partition = TRUE` (the default) makes a MISSING column fatal too: a
# master table that cannot state its own partition is exactly the drift this
# assertion exists to catch. Callers holding a legacy table can pass FALSE, and
# get a warning instead.
assert_partition_matches_cns <- function(cancer_type, partition,
                                         locked_types = CNS_LOCKED_TYPES,
                                         locked_label = CNS_LOCKED_LABEL,
                                         development_label = DEVELOPMENT_LABEL,
                                         require_partition = TRUE,
                                         context = "master table") {
  if (is.null(partition)) {
    msg <- sprintf(paste0("No `partition` column in %s, so the hardcoded locked list (%s) ",
                          "cannot be checked against it. Rebuild with scripts/build_master.py."),
                   context, paste(locked_types, collapse = "/"))
    if (require_partition) stop(msg, call. = FALSE)
    warning(msg, call. = FALSE)
    return(invisible(FALSE))
  }
  cancer_type <- as.character(cancer_type)
  partition   <- as.character(partition)
  if (length(cancer_type) != length(partition)) {
    stop(sprintf("cancer_type (%d) and partition (%d) differ in length in %s",
                 length(cancer_type), length(partition), context), call. = FALSE)
  }
  bad_na <- is.na(partition) | !nzchar(trimws(partition))
  if (any(bad_na)) {
    stop(sprintf("%d row(s) in %s have a missing/blank partition; the lock is unverifiable.",
                 sum(bad_na), context), call. = FALSE)
  }
  unknown <- setdiff(unique(partition), c(locked_label, development_label))
  if (length(unknown)) {
    stop(sprintf("Unrecognised partition label(s) in %s: %s. Expected only '%s' or '%s'.",
                 context, paste(unknown, collapse = ", "), locked_label, development_label),
         call. = FALSE)
  }
  by_type      <- cancer_type %in% locked_types
  by_partition <- partition == locked_label
  if (!identical(by_type, by_partition)) {
    locked_type_not_partition <- sort(unique(cancer_type[by_type & !by_partition]))
    locked_partition_not_type <- sort(unique(cancer_type[!by_type & by_partition]))
    stop(sprintf(paste0("LOCK MISMATCH in %s: the hardcoded CNS list and the `partition` column disagree.\n",
                        "  %d row(s) have cancer_type in {%s} but partition != '%s' (types: %s)\n",
                        "  %d row(s) have partition == '%s' but a cancer_type outside that list (types: %s)\n",
                        "  Refusing to proceed: which samples are locked must not be ambiguous."),
                 context,
                 sum(by_type & !by_partition), paste(locked_types, collapse = ", "), locked_label,
                 if (length(locked_type_not_partition)) paste(locked_type_not_partition, collapse = ", ") else "-",
                 sum(!by_type & by_partition), locked_label,
                 if (length(locked_partition_not_type)) paste(locked_partition_not_type, collapse = ", ") else "-"),
         call. = FALSE)
  }
  invisible(TRUE)
}


# -----------------------------------------------------------------------------
# tissue_identity_r2(): how much of a quantity is explained by tissue alone.
# -----------------------------------------------------------------------------
# R^2 of values ~ factor(cancer_type). For the PREDICTION this is the headline
# lineage-confounding number (56.2% in run 01); for the OBSERVED label it is the
# reference point that says how much of that is a property of HRD itself rather
# than of the model. Quoting either alone is misleading, so loco_merge.R now
# ships both on every run instead of the number existing only inside
# scripts/c1_zeroshot_percentile.R.
#
# Computed as 1 - SS_within/SS_total directly (identical to
# summary(lm(y ~ factor(cancer)))$r.squared but without fitting a design matrix
# with one column per tissue). Returns NA rather than erroring on degenerate
# input: zero total variance, fewer than two tissues, or fewer than three rows.
tissue_identity_r2 <- function(values, cancer) {
  ok <- is.finite(values)
  values <- values[ok]; cancer <- as.character(cancer)[ok]
  n <- length(values)
  if (n < 3L || length(unique(cancer)) < 2L) return(NA_real_)
  ss_tot <- sum((values - mean(values))^2)
  if (!(ss_tot > 0)) return(NA_real_)
  ss_within <- sum(vapply(split(values, cancer),
                          function(v) sum((v - mean(v))^2), numeric(1)))
  1 - ss_within/ss_tot
}

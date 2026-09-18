#!/usr/bin/env Rscript
# =============================================================================
# scripts/benchmark_v4_rank.R - time ONE V4 ranking pass at full probe count.
# =============================================================================
#
# PURPOSE
#   docs/27 section 6 pre-registers the V4 ranker but says nothing about its
#   cost, and a statistic that cannot finish a fold in the sentinel's time box
#   is not usable however correct it is. This measures the thing that actually
#   matters: ONE lineage_meta_stats() sweep over the real 336k-probe matrix on a
#   realistic outer-training block (29 tissues, all rows except one cancer),
#   plus the marginal cost of turning that sweep into a ranked top-5,000 at each
#   point of the lambda_pen grid.
#
#   No model is fitted and no result is produced. Timing harness only.
#
# USAGE
#   Rscript scripts/benchmark_v4_rank.R <beta.tsv> <master_samples.tsv> <out_dir>
#   (through bsub - loading the matrix needs far more memory than a login node.)
# =============================================================================

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) stop("Usage: benchmark_v4_rank.R <beta.tsv> <master_samples.tsv> <out_dir>")
beta_path <- args[1]; meta_path <- args[2]; out_dir <- args[3]
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

cat("host        :", Sys.info()[["nodename"]], "\n")
cat("started     :", format(Sys.time()), "\n\n")
source("R/model.R")

cat("reading matrix ...\n")
t_read <- system.time({
  beta <- data.table::fread(beta_path, data.table=FALSE, check.names=FALSE)
})[["elapsed"]]
x <- t(as.matrix(beta[,-1,drop=FALSE])); storage.mode(x) <- "double"; colnames(x) <- beta[[1]]
rm(beta); invisible(gc(verbose=FALSE))
cat(sprintf("  fread+transpose: %.1f s -> %d samples x %d probes\n", t_read, nrow(x), ncol(x)))

meta <- read.delim(meta_path, check.names=FALSE, stringsAsFactors=FALSE)
meta <- meta[match(rownames(x), meta$sample_id),]
# GBM/LGG stay locked: they enter neither the benchmark nor any statistic.
cns <- meta$cancer_type %in% c("GBM","LGG")
types <- sort(unique(meta$cancer_type[!cns]))
# Hold out BRCA (fold 3), the largest tissue and the worst case for this sweep.
tr <- !cns & meta$cancer_type != types[3]
cat(sprintf("held out %s; training n=%d over %d tissues\n\n",
            types[3], sum(tr), length(unique(meta$cancer_type[tr]))))

xt <- x[tr,,drop=FALSE]; yt <- meta$HRDsum[tr]; ct <- meta$cancer_type[tr]
rm(x); invisible(gc(verbose=FALSE))

cat("preprocess_base() (missingness filter + medians + imputation) ...\n")
t_base <- system.time({ base <- preprocess_base(xt) })[["elapsed"]]
cat(sprintf("  %.1f s -> %d probes survive\n", t_base, ncol(base$z)))

cat("lineage_meta_stats() ONE FULL RANKING PASS ...\n")
t_stat <- system.time({
  st <- lineage_meta_stats(base$z, yt, ct, med=base$med, block_size=20000L)
})[["elapsed"]]
cat(sprintf("  %.1f s (%.2f min)\n", t_stat, t_stat/60))

pens <- c(0,0.25,0.5,1,2)
t_fin <- vapply(pens, function(p) system.time({
  pp <- preprocess_finalize(base, 5000L, feature_rank="lineage_penalized",
                            cancer=ct, y=yt, lambda_pen=p, stats=st)
})[["elapsed"]], numeric(1))
cat(sprintf("  preprocess_finalize() per lambda_pen: %s s\n",
            paste(sprintf("%.1f", t_fin), collapse=", ")))

tot <- t_base + t_stat + sum(t_fin)
cat(sprintf("\nONE RANKING PASS (base + sweep + 5-point penalty grid): %.1f s = %.2f min\n",
            tot, tot/60))
# A fold pays this once for the outer refit and once per inner fold (29 of them).
cat(sprintf("projected per-fold ranking cost (1 outer + %d inner): %.1f min\n",
            length(unique(ct)), tot * (1 + length(unique(ct))) / 60))

write.table(data.frame(step=c("preprocess_base","lineage_meta_stats",
                              paste0("finalize_pen_", pens), "one_ranking_pass_total"),
                       seconds=c(t_base, t_stat, t_fin, tot)),
            file.path(out_dir, "v4_rank_benchmark.tsv"), sep="\t",
            row.names=FALSE, quote=FALSE)

# ---------------------------------------------------------------------------
# The number that actually sets the grid size: what ONE INNER FOLD costs, split
# into the part SHARED across the lambda_pen grid (base + meta sweep, paid once)
# and the part paid PER PENALTY (finalize + apply_preprocess + three glmnet
# paths). A 5-point grid is only affordable if the per-penalty part is small.
# ---------------------------------------------------------------------------
cat("\n--- inner-fold cost decomposition ---\n")
itr <- ct != sort(unique(ct))[1]           # drop one training tissue = one inner fold
xi <- xt[itr,,drop=FALSE]; yi <- yt[itr]; ci <- ct[itr]
xv <- xt[!itr,,drop=FALSE]
cat(sprintf("inner-train n=%d  inner-valid n=%d\n", nrow(xi), nrow(xv)))
i1 <- system.time({ bi <- preprocess_base(xi) })[["elapsed"]]
i2 <- system.time({ si <- lineage_meta_stats(bi$z, yi, ci, med=bi$med) })[["elapsed"]]
i3 <- system.time({ ppi <- preprocess_finalize(bi, 5000L, feature_rank="lineage_penalized",
                                               cancer=ci, y=yi, lambda_pen=0.5, stats=si) })[["elapsed"]]
i4 <- system.time({ zi <- apply_preprocess(xi, ppi, FALSE)
                    zv <- apply_preprocess(xv, ppi, FALSE) })[["elapsed"]]
alphas <- c(0.1,0.5,1)
lam <- lapply(alphas, function(a) lambda_path(zi, yi, a))
i5 <- system.time({ for (k in seq_along(alphas))
  glmnet::glmnet(zi, yi, alpha=alphas[k], lambda=sort(lam[[k]], decreasing=TRUE),
                 standardize=FALSE) })[["elapsed"]]
cat(sprintf("  base %.1f  sweep %.1f  finalize %.1f  apply %.1f  glmnet(3 alphas) %.1f\n",
            i1, i2, i3, i4, i5))
shared <- i1 + i2; per_pen <- i3 + i4 + i5
cat(sprintf("  SHARED per inner fold: %.1f s   PER PENALTY: %.1f s\n", shared, per_pen))
nin <- length(unique(ct))
for (G in c(1,2,3,5)) cat(sprintf("  grid of %d -> %.1f s/inner fold -> %.1f min over %d inner folds (+ outer refit)\n",
                                  G, shared + G*per_pen, nin*(shared + G*per_pen)/60, nin))
cat("\nfinished    :", format(Sys.time()), "\n")

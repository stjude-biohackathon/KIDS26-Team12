# 26 — C1: Unseen-Tissue Calibration Offset and the n-of-1 Workaround

**Status:** investigation complete for Phases A–B (real data), Phase C–D (synthetic
validation), and the **four sentinel folds, which RAN and FAILED their gates**
(LSF array `323195423`, `priority` queue, 4/4 complete 2026-09-17 ~21:00, zero
errors). See §7 for the scorecard and §8 for the revised verdict.

**Headline:** the full 30-fold rank array is **NOT authorised**. Three of seven
pre-registered gates failed, including both load-bearing ones. The relative-target
model is *more* tissue-confounded than the absolute model it was meant to replace.


**Scope.** This document addresses blocker C1: the model carries a tissue-specific
additive offset (mean |offset| = 3.46 HRD units) that cannot be estimated for a tumour
type the model has never seen, which breaks true n-of-1 inference. It asks two separate
questions:

1. If a small labelled calibration panel *were* available in a new tumour type, could a
   hierarchical (shrunk) estimator make k = 3–5 labels safe, where the current unshrunk
   mean is harmful?
2. With *zero* target labels, can the model emit a meaningful **relative** (within-tissue
   percentile) score instead of an absolute HRDsum?

**A percentile or relative score is not absolute HRDsum. It is not a clinical HRD
determination and is not a treatment-selection result.** Nothing in this document should
be read as supporting PARP-inhibitor decision-making.

---

## 1. Leakage boundaries and data provenance

| Constraint | How it was enforced |
|---|---|
| GBM/LGG partition stays locked | `stopifnot(!any(cancer_type %in% c("GBM","LGG")))` executes at the top of `scripts/c1_shrinkage_calibration.R` and again inside `scripts/c1_rank_one_fold.R` *before* any target transform, feature selection, fit, calibration, or summary. `R/rank_model.R::assert_cns_locked()` is the shared guard. |
| No PBTP / PBTA / pediatric / protected data | No such path was opened. Phases A–B read exactly one file. Phases C–D used synthetic matrices only. |
| Target-tissue labels never inform fitting | Every prior parameter, every mapping, and every preprocessing statistic is estimated with the target tissue fully excluded. Test 2 in `tests/test_c1_rank_model.R` proves that replacing the held-out tissue's labels with garbage leaves preprocess statistics, glmnet coefficients, and calibration knots byte-identical. |
| No interference with the concurrent log1p run | No `bsub`, `bkill`, or `bmod` was issued. Nothing under `results/loco_run01/` or `results/c2_log1p/` was written. `docs/21_BLOCKER_RESOLUTION_PLAN.md` and `docs/23_PIPELINE_SCHEMATIC.md` were not modified. |

### Inputs

- `results/loco_run01/loco_predictions.tsv` — 7,065 LOCO predictions across 30 non-CNS
  tumour types. **Sole real-data input for Phases A and B.**
- Phases C–D: synthetic matrices generated in-script (300 × 400, 6 tissues).

### Outputs

```
results/c1_shrinkage/shrinkage_priors.tsv
results/c1_shrinkage/shrinkage_by_tissue_k.tsv
results/c1_shrinkage/shrinkage_macro_by_k.tsv
results/c1_shrinkage/shrinkage_draw_distribution.tsv
results/c1_shrinkage/shrinkage_by_offset_stratum.tsv
results/c1_rank_probe/zeroshot_percentile_by_tissue.tsv
results/c1_rank_probe/zeroshot_percentile_macro.tsv
results/c1_rank_probe/zeroshot_reliability_bins.tsv
results/c1_rank_probe/zeroshot_tissue_r2.tsv
```

### Code

```
scripts/c1_shrinkage_calibration.R   Phase A
scripts/c1_zeroshot_percentile.R     Phase B
R/rank_model.R                       Phase C + D library
scripts/c1_rank_one_fold.R           Phase C LOCO fold driver
tests/test_c1_rank_model.R           12-case leakage / invariance suite
```

---

## 2. Phase A — empirical-Bayes shrinkage of the few-shot offset

### Method

For each target tissue *t*, tissue *t* is excluded from all prior estimation. From the
remaining 29 tissues, using per-sample residual `r = predicted − actual`:

- `mu0` = equal-weight mean of per-tissue mean residuals (macro, not sample-weighted)
- `sigma2` = equal-weight mean of within-tissue residual variances
- `tau2 = max(0, Var_s(rbar_s) − mean_s(s2_s / n_s))` — between-tissue offset variance,
  corrected for finite-sample offset noise (uncorrected `tau2_raw` retained as a
  sensitivity column)

With k labelled target samples, `w = tau2 / (tau2 + sigma2/k)` and
`b_hat = w * rbar_t + (1 − w) * mu0`.

1,000 draws per (tissue, k); k ∈ {1, 3, 5, 10, 20}; eligibility `n_t ≥ k + 10`; scoring
only on rows not used for calibration. OV (n = 10) is ineligible at every k.

### Results (macro over 29 tissues, equal weight per tissue)

| k | m1 uncorr | m2 unshrunk | m3 **shrunk** | m4 oracle | m5 flat null | m6 **hier. null** | w̄ | helped | beat m5 | beat m6 | oracle gain (shrunk) | (unshrunk) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 8.434 | 10.548 | **8.217** | 7.444 | 12.028 | 11.213 | 0.167 | 16/29 | 27 | 25 | 0.219 | −2.137 |
| 3 | 8.434 | 8.649 | **8.040** | 7.444 | 9.983 | 10.111 | 0.375 | 14/29 | 22 | 26 | 0.398 | −0.218 |
| 5 | 8.433 | 8.214 | **7.934** | 7.444 | 9.518 | 9.686 | 0.500 | 15/29 | 21 | 25 | 0.504 | 0.222 |
| 10 | 8.432 | 7.853 | **7.775** | 7.442 | 9.132 | 9.263 | 0.666 | 17/29 | 20 | 23 | 0.664 | 0.585 |
| 20 | 8.434 | 7.665 | **7.653** | 7.444 | 8.945 | 9.024 | 0.799 | 18/29 | 20 | 21 | 0.789 | 0.777 |

`m6` is the fair hierarchical null: a label-only predictor given the *same* source-tissue
prior information and the *same* k target labels, shrinking `mean(shot actual)` toward the
cross-tissue mean of tissue means.

Variance components (median across the 30 leave-one-tissue-out fits): **τ² = 22.3**
(τ²_raw = 23.9; the finite-sample correction removes ~7%), **σ² = 108.6**, τ²/σ² ≈ 0.205,
so `w ≈ k / (k + 4.87)`. `mu0` ranges only 0.63–1.45 — the pan-cancer offset is close to
zero, so the prior is effectively "apply almost no correction".

### Stratification by true offset magnitude

| stratum | tissues | k=3 m1 | m2 | m3 | helped |
|---|---|---|---|---|---|
| \|offset\| ≥ 5 | 6 | 11.10 | 9.01 | 9.37 | **6/6** |
| \|offset\| 2–5 | 12 | 8.03 | 8.54 | 7.79 | 8/12 |
| \|offset\| < 2 | 11 | 7.42 | 8.57 | 7.59 | **0/11** |

### Interpretation

Shrinkage converts the few-shot correction from harmful to beneficial at small k. At k=1
the unshrunk estimator is catastrophic (oracle gain −214%); shrunk it is mildly positive.
At k=3 unshrunk remains harmful (−22% of oracle gain) while shrunk recovers 40%. At k=5,
50% vs 22%.

Three honest qualifications:

- **Macro-mean safety is not per-tissue safety.** Only 14/29 tissues are helped at k=3 and
  15/29 at k=5 — a minority. Per draw, roughly 62% of k=3 draws improve on no correction.
  The macro gain is carried by large wins in a few high-offset tissues.
- **All the benefit lives in the high-offset stratum.** For the 11 tissues with |offset| < 2,
  shrinkage never helps at any k. It does, however, bound the damage (7.59 vs 7.42 for no
  correction) where the unshrunk estimator inflicts 8.57. In the |offset| ≥ 5 stratum the
  shrunk estimator *over-shrinks* a genuinely large offset and is beaten by the unshrunk
  mean (9.37 vs 9.01). A magnitude-adaptive rule would dominate uniform shrinkage.
- **The decisive comparison passes.** Against the fair hierarchical null m6 the shrunk model
  wins in 25–26 of 29 tissues at every k (8.04 vs 10.11 at k=3). The genomic model carries
  within-tissue signal that no amount of prior information plus k labels can reproduce.

**This does not enable zero-label inference.** Every method in this table requires k ≥ 1
labelled samples from the target tumour type. The purpose of Phase A is solely to
determine whether a small future pediatric/CNS calibration panel could be used *safely*.
The answer is: yes at k ≥ 3 with shrinkage, provided the panel is used with the
understanding that it helps most where the offset is large and is near-neutral where it is
small.

---

## 3. Phase B — zero-shot percentile triage

### Method

Leave-one-tissue-out over all 30 tissues. In each source tissue, the true label is
converted to a midrank percentile `q = (rank(actual, ties="average") − 0.5) / n_tissue`.
A frozen mapping `raw prediction → q` is fitted on source rows only, with equal total
weight per source tissue, in two variants: weighted linear, and weighted isotonic (a
hand-written weighted PAVA, unit-checked against `stats::isoreg` at equal weights). The
mapping is applied to the target tissue sample-by-sample. Target `q` is used for
evaluation only.

**Single-sample invariance was asserted, not assumed:** scoring the target tissue one row
at a time reproduces block scoring exactly, for all 30 tissues and all three mappings.

### Results

| metric | linear | monotone | null |
|---|---|---|---|
| macro MAE in q | **0.2263** | 0.2281 | 0.2452 (constant 0.5) |
| tissues beating null | 27/30 | 25/30 | — |
| macro Spearman | 0.4735 | 0.4560 | 0 |
| top-quartile AUROC (29 tissues) | 0.729 | 0.710 | 0.5 |
| Brier (frozen source-only logistic) | 0.1740 | — | 0.1835 (constant 0.25) |
| tissues beating Brier null | 21/30 | — | — |
| pooled calibration slope / intercept | 1.104 / −0.062 | 0.986 / −0.003 | — |
| macro per-tissue calibration slope | 2.61 | 2.22 | — |
| tissue-identity R² of output | **0.593** | 0.513 | q_true = 0.000 |
| tissue-identity R² of raw prediction | 0.562 | | |

### Interpretation

The mapping beats its null on MAE (7.7% reduction) and gives a usable top-quartile AUROC
of 0.73. But three findings undercut any stronger reading:

- **The mapping adds exactly zero ranking power.** Macro Spearman for the linear map equals
  the Spearman of the raw prediction against the actual label within tissue, to the last
  digit (max per-tissue gap = 0). A strictly increasing transform cannot reorder anything.
  The isotonic map scores *worse* (0.4560) purely because flat PAVA blocks create ties that
  destroy ranking information. The mapping only relocates the absolute placement.
- **Within-tissue calibration is poor.** Pooled slope is near 1 only because tissue-level
  errors cancel. The macro per-tissue slope is 2.61 and the intercept −0.82: within any one
  tissue, q̂ is badly under-dispersed. 2,666 of 7,065 samples fall in the single bin (0.4, 0.5].
- **Tissue identity is not removed.** A true within-tissue percentile must have tissue R² ≈ 0
  (q_true gives exactly 0.000). The linear map gives **0.593**, slightly *higher* than the raw
  prediction's 0.562 on the same pooled set. The "percentile" is still mostly reporting which
  tumour type the sample came from.

**Failures:** linear fails against the null in CHOL, KIRC, UVM; monotone additionally in OV
and UCS. CHOL is anti-correlated (ρ = −0.134, AUROC 0.463, n = 35). UVM (ρ = 0.154) and
PCPG (ρ = 0.064) carry essentially no within-tissue signal, consistent with run 01.

> **WARNING — stacked diagnostic.** The predictions consumed here come from 30 *different*
> LOCO outer-fold models. The mapping is therefore stacked on heterogeneous inputs, and no
> single frozen model produced this table. This is a feasibility screen, not a validated
> frozen-model result. The warning is printed to stdout by the script and repeated in its
> header.

---

## 4. Phase C — leakage-safe relative-target model (synthetic validation only)

`R/rank_model.R` implements a model separate from the production `R/model.R`, which is
unmodified. Two pointwise relative targets:

- **centered:** `y − mean(y within training tissue)`
- **normal_score:** `z = qnorm(q)` on the within-tissue midrank percentile, with defensive
  clipping at the half-rank limits and an explicit constant-tissue branch mapping to
  q = 0.5, z = 0 rather than NaN.

Tissue statistics are computed from training rows only and recomputed inside every inner
validation fold from that fold's inner-training partition alone. Validation-row statistics
are used for *scoring* validation rows only and are labelled `EVAL-ONLY` in code; nothing
learned from them enters the fit. Hyperparameter selection is macro-averaged over tissues
**in relative space** (macro MAE in z, with macro 1 − Spearman as a secondary diagnostic);
an explicit guard rejects an absolute-HRDsum tuning objective. The frozen bundle carries
`score_type = "relative_score"` and deliberately stores nothing permitting back-transform
to absolute HRDsum.

Variants: **V1** pooled features + relative target + sample weights; **V2** = V1 with equal
total glmnet weight per training tissue; **V3** = V2 with pooled *within-tissue* variance
feature ranking.

### V3 feasibility

Within-tissue variance is computed by the blocked sufficient-statistics identity
(`colSums(x_t)`, `colSums(x_t²)` per tissue; within-SS = `Σ_t (SS_t − S_t²/n_t)`), never
materialising a residualised copy. Measured at 7,000 × 2,000 × 28 tissues: **0.79 s**.
Extrapolated: **~8 s per 20,000-column block, ~2.2 min for all 336,480 probes, ~1.04 GB per
block (~2.5 GB peak)**. Negligible against the ~45 min/fold matrix load. There is no silent
fallback — the path errors rather than degrading to a pooled or leaky computation.

### Synthetic evidence (held-out tissue T6, never seen in fitting)

| scenario | V1 | V2 | V3 | centered |
|---|---|---|---|---|
| pure tissue intercept, no within-tissue signal (seed 808) | 0.000 (constant score) | 0.000 (constant) | — | — |
| shared within-tissue signal (seed 909) | **0.964** | **0.964** | **0.968** | 0.961 |
| negative control, noise features (seed 1010) | — | 0.000, permutation p = 1.00 | — | — |

End-to-end via the fold script on a 4-tissue TSV fixture: held-out Spearman 0.931 against a
permutation null of [−0.346, 0.318], p = 0.002.

The relative target transfers to an unseen tissue on synthetic data and manufactures nothing
when there is nothing to find. **This is a correctness check, not an effect-size forecast** —
synthetic linear signal with a β shared exactly across tissues is the easiest possible version
of this hypothesis.

---

## 5. Phase D — within-tissue pairwise ranker (contingency prototype)

Pairs are constructed within training tissue only, features `x_i − x_j`, outcome
`I(y_i > y_j)`, penalised logistic regression with per-tissue equal total weight, bounded
and balanced pair sampling, single sample scored as β'x (intercept asserted exactly zero).
Synthetic: signal ρ = 0.957, noise ρ = −0.157; pair counts bounded, outcome mean in
[0.4, 0.6] globally and per tissue.

It works, but gives no advantage over V2/V3 on synthetic data at O(pairs) additional cost.
**Recommendation: retain as a documented contingency. Escalate only if the pointwise relative
target underperforms on real data.** Per the sentinel policy, no full pairwise array should be
launched without pointwise evidence.

---

## 6. Tests

`Rscript tests/test_c1_rank_model.R` — **12/12 pass, 7.8 s.**

Coverage: GBM/LGG cannot enter fitting or summaries (including a `tryCatch` confirming the
assertion *errors* when CNS rows are forced in); held-out tissue labels cannot affect target
transforms or calibration parameters; permuting held-out labels leaves predictions identical;
midrank behaviour with ties and constant tissues; equal total weight per tissue; one-sample
scoring equals batch scoring; feature-order invariance; pure tissue intercept handled without
manufacturing performance; within-tissue signal transfers to an unseen tissue; negative control;
pairwise sampling bounded/balanced/within-tissue; tuning objective is a relative-space loss.

Two tests exposed real problems, both fixed in the code rather than by weakening the test:

1. On pure-noise data the tuner correctly shrank to an intercept, making `relative_score`
   constant and `cor()` return `NA` — which would have made a *correct refusal to order* look
   like a crash. `rank_null_panel()` now has an explicit constant-score branch (ρ = 0, p = 1)
   and the harness reports `constant_score` separately, so "no claim" is distinguishable from
   "a wrong claim".
2. Test 9's fixture, not the model: signal probes drawn i.i.d. uniform like the background never
   survived the *unsupervised* top-variance filter (ρ = 0.107). The filter behaved as designed;
   the fixture was unrealistic. Informative probes are now drawn bimodally (0.1/0.9 mixture), as
   real informative CpGs are.

---

## 7. Sentinel results — RAN 2026-09-17, GATES FAILED

LSF array `323195423`, `priority` queue, variant V2 (`normal_score` / `tissue` /
`pooled`). All four folds completed, zero errors, every stderr empty, every
selected lambda **interior** to its data-derived path (alpha = 1 in all four).
Fits took 64–86 min against run 01's 45–54 min; the difference is node
contention on the shared `priority` hosts (49–64 jobs, load 41–60), not extra
model cost. Peak RSS 129–174 GB against 240 GB reserved.

```bash
# Reproduce:
bsub < scripts/lsf_c1_rank_sentinel.bsub
```

### Per-fold results against run 01

| Fold | n | Rank ρ | Run 01 ρ | Δρ | Run 01 r | Perm p |
|---|---:|---:|---:|---:|---:|---:|
| UCEC | 403 | 0.550 | 0.703 | **−0.154** | 0.758 | 0.002 |
| BRCA | 743 | 0.261 | 0.666 | **−0.405** | 0.632 | 0.002 |
| KICH | 65 | 0.443 | 0.543 | −0.101 | 0.803 | 0.002 |
| THCA † | 464 | 0.125 | −0.029 | +0.154 | −0.032 | 0.004 |

† **THCA is the negative control, not a win.** Its within-tissue target is
near-constant, so the rank transform is degenerate there. Its apparent
"improvement" is exactly the artefact the gate was written to catch — counting it
would repeat the C2/log1p preview error (see B16).

Macro ρ **excluding THCA**: **0.418** (rank) vs **0.637** (run 01), Δ = **−0.219**.

### Gate scorecard

| # | Gate | Result | Verdict |
|---|---|---|---|
| 1 | Macro discrimination not materially worse than run 01 | ρ 0.418 vs 0.637, **−0.219** | ❌ **FAIL** |
| 2 | BRCA shows no large discrimination loss | ρ 0.261 vs 0.666, **−0.405** | ❌ **FAIL** |
| 3 | THCA not presented as successful | flagged degenerate, excluded from macro | ✅ honoured |
| 4 | Tissue identity explains substantially less than 56.2% | **57.6%** vs 52.7% for the absolute model on the same samples | ❌ **FAIL** |
| 5 | ≥3 of 4 folds show positive rank association | 4/4 positive | ✅ PASS |
| 6 | Beats an honest relative-target null | perm p ≤ 0.004 in all folds | ✅ PASS |
| 7 | All leakage / n-of-1 invariance tests pass | 12/12, plus 7/7 (B13) and 10/10 (B14) | ✅ PASS |

**Three gates failed, including both load-bearing ones. The full 30-fold rank
array is NOT authorised.**

### Tissue-identity R², apples to apples

Same four tissues, same samples, both quantities pooled across the four held-out
folds:

| Quantity | Tissue R² |
|---|---:|
| Rank model `relative_score` | **0.576** |
| Run 01 absolute prediction | 0.527 |
| **True relative target (the floor)** | **0.001** |

This is the decisive number. A genuine within-tissue score should approach the
0.001 floor. The rank model sits at 0.576 — *higher* than the absolute prediction
it was designed to improve on.

### What this means, mechanistically

The model learned something real: every fold beats its within-tissue permutation
null at p ≤ 0.004, and the negative controls in `tests/test_c1_rank_model.R`
confirm the machinery does not manufacture signal. But it learned **less** than
the absolute model, and it **did not remove the tissue confound it was built to
remove**.

The likely mechanism is feature selection. The 5,000-probe filter ranks by
**pooled** variance and runs *before* the target is ever consulted, so it
preferentially selects lineage-discriminating probes. Changing the target to a
within-tissue quantity cannot undo a feature set already chosen for between-tissue
variance. **V3 (`feature_rank = "within_tissue"`) is the untested variant that
addresses exactly this**, and its blocked implementation is already written and
benchmarked (~2.2 min for all 336,480 probes). It is the one remaining pointwise
option that has not been falsified.

This also confirms, with a properly trained model, what Phase B found with a
post-hoc mapping: **a monotone remap of a tissue-confounded score is still a
tissue-confounded score, and training on a relative target is not by itself
sufficient to fix it.**

---

## 7b. What C1 means for the pediatric transfer goal

This is the section that matters for the flash talk, because pediatric HGG
transfer is the project's stated purpose. **The failure closes both routes to an
n-of-1 pediatric answer, for the same underlying reason.**

### Route 1 — report an absolute HRDsum. Blocked.

The per-tissue offset is not estimable without target-domain labels. Across the
30 development tissues the offsets span **−15.44 (SARC) to +8.29 (PCPG)**, mean
|offset| 3.46, SD 4.83. Treating those 30 as the predictive distribution for the
next unseen tissue — which is what they are, since label-free regression on tissue
covariates failed at LOTO R² = −0.116 and the offset is uncorrelated with every
available covariate (|r| ≤ 0.267) — gives a

> **95% predictive interval for a new tumour type's offset: −9.7 to +8.1 HRD
> units, a width of 17.8 units — 74% as wide as the interquartile range of the
> label itself (4 to 28).**

A pediatric tumour type is **one draw** from that distribution, and nothing we can
measure narrows it. GBM and LGG remain locked, so their offsets are unknown; the
nearest available lineage proxies are not reassuring (SARC −15.4, UVM −2.3,
TGCT +2.0, THCA +8.0, PCPG +8.3).

**On the scale of the decision.** At the exploratory threshold of 42, the
population-wide flip rate looks mild at ~1.5% — but that is an artefact of where
the threshold sits, since only 10.6% of the cohort exceeds it. The operative
question is conditional:

> **Among the 1,171 patients within ±10 units of the threshold — precisely the
> patients for whom a test is supposed to add information — a tissue-level offset
> flips the HRD-high call for a median of 9.1% of them, and up to 79.2% in the
> worst tissue.**

Roughly one borderline patient in eleven would be reclassified by a constant
belonging to their tumour type rather than by their own biology.

### Route 2 — report a within-tissue rank instead. Now also blocked.

This was the designed escape hatch: decline the absolute claim and report ordering
within the tumour type, which is offset-invariant by construction. **Today's
sentinel closed it.** The relative score is 57.6% tissue identity against a 0.1%
floor, and BRCA discrimination fell by 0.405.

**This matters more for pediatrics than for the TCGA tissues.** 89% of the cohort
sits below the threshold, and pediatric HGG is expected to sit in that low-HRD
regime, where a binary "HRD-high" call is almost never the useful output and
**relative ordering is** — which is exactly the output now shown to be majority
tissue identity.

### What survives

**The shrinkage result (§2) is the one intact path**, and it is a real one: with
~10 labelled samples from the new tumour type, k=10 recovers 58% of achievable
gain, and EB shrinkage makes k=3 safe where the unshrunk mean is harmful. For a
pediatric cohort this is viable — but it **requires PBTP labels**, which we do not
have and did not touch, and it is **not n-of-1**.

### Honest statement for the talk

> We can rank tumours within a known type, and we can calibrate to a new type
> given ~10 labelled examples. We cannot yet place a single patient from an unseen
> tumour type on an absolute HRD scale, and our attempt to sidestep that with a
> relative score produced an output that is 58% tissue identity. The CNS partition
> is still sealed, because the calibration problem that would make a CNS number
> interpretable is not solved.

---

## 7c. Sentinel plan as originally pre-registered (retained for audit)

Recorded before the run, unchanged, so the scorecard above cannot be read as
post-hoc:

Four outer folds only — **BRCA, UCEC, THCA, KICH** — chosen a priori to span the
failure modes: BRCA as the largest tissue, KICH as the strongest run-01
within-tissue signal in a small cohort, UCEC as high-signal with the largest
positive offset, THCA as a **negative control** whose near-constant target makes
rank evaluation degenerate.

Gates: macro discrimination not materially worse than run 01; no large BRCA loss;
THCA not counted as success; tissue identity substantially below 56.2%; ≥3 of 4
folds positive; beats an honest relative-target null; all leakage tests pass.


---

## 8. Verdict

**Which is supported: (A) hierarchical few-shot calibration, (B) a zero-shot relative score,
(C) both, or (D) neither?**

**A only. B is now REFUTED, not merely unproven.**

- **A — supported with caveats.** Empirical-Bayes shrinkage rescues k = 1, 3, and 5 in the macro
  mean, converting a harmful correction into a beneficial one, and beats the fair hierarchical
  null in 25–26 of 29 tissues. But it helps only a minority of individual tissues at k = 3–5, and
  the entire benefit is concentrated in tissues whose true offset is large. It requires labels and
  therefore does **not** solve n-of-1 inference. Its value is that a future small pediatric or CNS
  calibration panel (k ≥ 3, ideally k ≥ 10) could be used without making predictions worse.
- **B — REFUTED on the evidence now available.** The earlier draft of this document recorded B as
  "unproven, not refuted", pending the sentinel folds. **The sentinel folds have run.** Both
  implementations of a zero-shot relative score have now failed on the same axis:
  - the post-hoc percentile remap (§3) adds *zero* ranking information and raises tissue R² from
    0.562 to 0.593;
  - the properly trained relative-target model (§7) loses 0.219 macro ρ against run 01, loses
    0.405 on BRCA, and raises tissue R² to **0.576** against a 0.001 floor.

  Two independent methods, the same failure. Training on a within-tissue target is **not
  sufficient** to remove lineage when the feature set was selected for between-tissue variance.

**One pointwise option remains untested and is not covered by this refutation: V3**
(`feature_rank = "within_tissue"`), which replaces pooled-variance probe ranking with
training-only *within-tissue* variance ranking. That addresses the diagnosed mechanism directly
rather than working around it, the blocked implementation is written and benchmarked (~2.2 min
across 336,480 probes), and it has not been run on real data. It is the single highest-value next
experiment for C1c. The Phase D pairwise ranker remains a contingency behind it.

### What remains unidentified without target-domain labels

The additive tissue offset is not identifiable from source-domain data alone. Run 01 already showed
label-free offset regression failing (best LOTO R² = −0.116), and Phase A explains why: the
cross-tissue offset spread (τ² = 22.3) is small relative to within-tissue residual noise
(σ² = 108.6), so there is little between-tissue structure for a label-free predictor to exploit,
and the pan-cancer mean offset is near zero. Phase B shows that a monotone remap cannot recover it
either, and §7 shows that a trained relative target does not either.

**Any absolute HRDsum reported for an unseen tumour type carries an unquantified additive bias with
a 95% predictive interval of −9.7 to +8.1 HRD units** (§7b) — 74% as wide as the label's own
interquartile range. A relative/percentile output sidesteps the offset by refusing to make an
absolute claim, but the sentinel shows our current relative output is majority tissue identity, so
it does not yet deliver a trustworthy within-tissue ordering either.


### Standing caveats

1. Held-out performance is measured against a relative truth that requires the held-out tissue's own
   labels. That quantity is not computable at deployment, so reported ρ is an **upper bound** on what
   an n-of-1 user learns.
2. Inner-fold hyperparameter selection is graded against validation-row tissue statistics. Nothing
   flows into the fit (test 2), but tuning is selected under a slightly more favourable regime than
   deployment.
3. Percentile calibration assumes exchangeability between source tissues and the new tissue —
   untestable from source data and exactly the assumption that a pediatric tumour type is likely to
   break. `rule=2` clamping makes extreme scores report as "at or beyond the most extreme source
   score" rather than extrapolating.
4. Batch and purity confounding are absent from the synthetic fixtures, so a real sentinel run could
   still ride purity rather than scar biology. C2's purity exoneration applies to the production
   model, not automatically to the relative-target model.
5. Nothing upstream of `beta.tsv` is exercised by these tests.

### Reiteration

A percentile or relative score is **not** absolute HRDsum, **not** a clinical HRD determination, and
**not** a treatment-selection result. The GBM/LGG partition remains locked and unopened. No PBTP,
PBTA, or pediatric data was accessed at any point in this work.

# 30 — Final analysis, 2026-09-18

Written 2026-09-18, revised 2026-09-22. This is the terminal document for the
hackathon phase. It supersedes `docs/24` as the description of the *shipped*
model, without replacing it as the record of run 01.

---

## 1. Executive summary

A frozen DNA-methylation elastic net, trained on 7,065 adult TCGA specimens
across 30 non-CNS cancer types and never exposed to CNS tissue, was scored
exactly once against 642 locked GBM/LGG specimens.

**Within-tissue rank transferred, weakly. Absolute calibration did not.**

GBM Pearson 0.320 (95% CI [0.173, 0.461]), LGG 0.329 ([0.234, 0.421]) — both
exclude zero, so there is genuine HRD-associated methylation signal in a lineage
the model had never seen. But both sit at the 13th–17th percentile of the 30
source-domain held-out tissues, and **the oracle tissue-mean null beats the model
in every group** (skill −1.009 GBM, −0.048 LGG, −0.233 pooled). Knowing only
"this is GBM, use the GBM average" outperforms the model's absolute predictions.

This is outcome **A** of the three pre-declared in `docs/28` §7. The deployment
gate remains shut. Nothing here is a validated clinical HRD assay.

The single most informative result is a negative one: a feature-selection change
that removed **42%** of the model's excess lineage structure *within* the
training distribution produced **no** improvement in calibration *outside* it.

---

## 2. Scientific question

Can a frozen methylation model predict independently measured reference HRDsum in
a cancer type it has never seen — ultimately pediatric high-grade glioma?

Canonical HRDsum = HRD-LOH + LST + TAI, measured from **allele-aware** genomic
data. Methylation arrays observe neither LOH nor telomeric allelic imbalance. The
model is therefore a *predictor of an independently measured continuous reference
value*, never a re-derivation of it. This framing is load-bearing and appears in
every claim below.

---

## 3. Data

| | |
|---|---|
| Matrix | 336,480 probes × 7,707 samples, 28 GB, probes in rows (sha256 `e3642d30…`) |
| Development | 7,065 specimens, 30 non-CNS cancer types |
| Locked CNS | 642 specimens — **135 GBM, 507 LGG** |
| Label | HRDsum; median 14, IQR 4–28; 10.6% ≥ 42; 14.3% exactly 0 |
| Probe bridge | 384,640-probe HM450/EPIC-v1 allowlist (sha256 `f359e43b…`) |

There is **no separate CNS matrix**; all 7,707 samples share one file, so scoring
CNS requires a full matrix load.

**Limitations that predate this analysis.** OV is n=10 — TCGA ovarian methylation
is mostly 27k and excluded by the 450k bridge, so the canonical HRD cancer is
effectively absent. HRD labels join at 15-character sample-type level, not vial.
Per-specimen array QC has never gated this cohort.

---

## 4. Modelling strategy

Elastic net on unsupervised-filtered CpGs → reference HRDsum. Nested
leave-one-cancer-out: the outer loop holds out an entire cancer type, and **inner
folds are one per training cancer type**, so hyperparameters are selected for
cross-tissue generalisation rather than within-tissue fit.

`alpha ∈ {0.1, 0.5, 1}`, data-derived lambda path anchored to each fold's own
`lambda.max`, `standardize=FALSE` (inputs pre-standardised), 5,000 features,
seed 260910. Adopted C2 rule: predictions clipped at 0, since negative HRDsum is
physically meaningless.

---

## 5. Leakage controls

Every preprocessing statistic — missingness filter, imputation medians, feature
ranking, centring, scaling — is relearned inside each training fold, outer and
inner, on that fold's rows only. `tests/test_v3_feature_rank.R` and
`tests/test_v4_lineage_penalized.R` assert this empirically: replacing held-out
feature rows *and* labels with garbage leaves preprocessing constants,
coefficients, alpha and lambda bit-identical.

The lock was enforced in code (`assert_partition_matches_cns()` requires two-way
agreement between the hardcoded `c("GBM","LGG")` and the `partition` column) and
by an append-only ledger on the inference path.

**One honest imperfection**, found by audit and recorded rather than quietly
fixed: the lambda grid is derived from all outer-training rows including the
inner-validation rows, so inner-CV model selection is not strictly nested. This
is a within-fold issue, not a held-out-cancer leak, and it affects both
candidates identically.

---

## 6. Run 01 — the original result

LOCO array `323078995`, 30/30 folds. Within-tissue Pearson 0.612 / Spearman
0.601 by the pooled-within estimator; permutation p = 0.001; positive in 29/30
tissues. Skill over the tissue-mean null only **+0.085**.

Per-tissue within-tissue r ranged from KICH 0.803 and UCEC 0.758 down to
**THCA −0.032, PCPG 0.110, CHOL 0.161** — several tissues carry no usable signal.

**The defect that defined the rest of the project:** tissue identity explained
**55.8%** of the *prediction* but only **34.1%** of the *truth*. The model encoded
substantially more lineage than the quantity it was predicting, and per-tissue
calibration offsets ranged SARC −15.44 to PCPG +8.29 (mean |offset| 3.46).

---

## 7. The main failure mode (C1)

A per-tissue additive offset. Removing it costs more than the model's entire
margin over the null. The diagnosed mechanism (`docs/26` §8): the 5,000-probe
filter ranked CpGs by **pooled** variance *before the target was consulted*, and
between-tissue methylation differences dwarf within-tissue ones — so the filter
selected lineage markers by construction.

---

## 8. Blocker-resolution experiments — five attacks, four failures

Reported in full. They are the scientific content, not an appendix.

| # | Attack | Outcome |
|---|---|---|
| 1 | Label-free tissue-offset regression | **FAILED** — best LOTO R² = −0.116, worse than predicting the mean |
| 2 | Zero-shot percentile remap | **FAILED** — adds exactly zero ranking information (a monotone map cannot reorder) and *raised* tissue R² 0.562 → 0.593 |
| 3 | Relative-target V1/V2 | **FAILED** 3 of 7 pre-registered gates (macro ρ −0.219, BRCA −0.405) |
| 4 | V4 supervised lineage-penalized selection | **FAILED** the mechanism gate G3 — tissue R² went **UP**, 0.5219 → 0.5361 |
| 5 | **V3-abs** within-tissue variance ranking | **PASSED** 7/7 sentinel gates and 7/7 selection criteria |

**V4 deserves emphasis.** It was the one pre-registered supervised alternative,
designed explicitly to penalize lineage, and it produced *more* lineage
imprinting than the unsupervised baseline. Its `lambda_pen` was selected by inner
folds and came out **0 for BRCA** — on the largest fold, the tuner preferred no
penalty at all. Plausible mechanism, offered as hypothesis: ranking by
|z_meta| × sign-consistency rewards probes correlating with HRDsum in the *same
direction* across tissues, but between-tissue HRD differences are themselves
lineage-structured, so consistency and lineage are not independent.

**What does work, but needs labels:** empirical-Bayes shrinkage few-shot
calibration (τ² = 22.3, σ² = 108.6). k=3 recovers +40% of oracle gain, k=10
recovers 58%. This is the *labelled* path — it is not n-of-1 inference.

---

## 9. V3-abs — definition and source-domain result

**One change from run 01:** the 5,000-probe filter ranks by pooled **within-tissue**
variance, computed on training-fold rows only, instead of pooled **total** variance.

$$\mathrm{Var}_{\text{within}}(j) = \frac{1}{N-T}\sum_{t}\left(Q_{t,j} - S_{t,j}^2/n_t\right)$$

**This is not `docs/26`'s "V3"**, which stacked that filter on V2's already-refuted
relative target. The distinction was pre-registered in `docs/27` §3 before any
result existed, precisely so the two would not be confused.

### Full 30-fold head-to-head (`results/tables/tableA2_A_vs_V3_full30.md`)

| Metric | A (run 01) | V3-abs | Δ |
|---|---|---|---|
| Macro within-tissue Pearson | 0.5197 | 0.5036 | −0.0161 |
| Macro within-tissue Spearman | 0.4747 | 0.4472 | −0.0275 |
| Pooled MAE | 8.9717 | 8.4456 | −0.526 |
| Mean absolute tissue bias | 3.4628 | 2.9737 | −0.489 |
| **Tissue R² of predictions** | **0.5584** | **0.4672** | **−0.0913** |
| Tissue R² of truth | 0.3414 | 0.3414 | cancels exactly |
| Macro skill vs oracle null | 0.419 | 0.630 | +0.211 |
| Tissues with positive r | 29/30 | 29/30 | — |

**Headline.** ΔR²_excess = **0.0913**, 95% CI [0.0824, 0.0996], p < 0.0001 —
**42.0%** (CI 37.9–46.5%) of excess lineage structure removed. Since both models
predict the same labels on the same cohort, R²_truth is identical and cancels, so
this quantity equals R²_pred,A − R²_pred,B.

**The caveat that must travel with it.** The paired per-tissue ranking change is
**not** statistically distinguishable from zero: Pearson mean −0.0161
(CI [−0.0410, 0.0051], Wilcoxon p = 0.428), Spearman −0.0275 (p = 0.477); 12/30
and 14/30 tissues improved. The evidence is asymmetric — a large, precisely
estimated lineage reduction against a *null result* on ranking, which is absence
of evidence for harm, not evidence of absence.

**V3-abs is a targeted reduction in lineage imprinting at no measurable ranking
cost. It is not a general accuracy improvement.**

### ROC, both framings

| Framing | A | V3-abs |
|---|---|---|
| Pooled, exploratory HRDsum ≥ 42 | 0.862 [0.850, 0.875] | 0.866 [0.855, 0.878] |
| **Tissue identity ONLY (control)** | **0.743** | **0.723** |
| Within-tissue top quartile | 0.779 [0.767, 0.791] | 0.777 [0.764, 0.789] |

The control curve is the point: replacing each sample with its own tissue's mean
prediction — no within-patient information whatsoever — still reaches ~0.74
against the pooled cutoff. **Most of the pooled AUC is lineage recognition.** The
within-tissue figure is the defensible claim, and the control falling
0.743 → 0.723 while within-tissue holds is the ROC view of the same lineage
reduction.

---

## 10. HRD component decomposition (secondary, exploratory)

Label-side only; no model fitting. HRDsum = LOH + LST + TAI holds exactly for all
7,065 development rows.

The prediction tracks **TAI** best (macro within-tissue r 0.543, vs LOH 0.427,
LST 0.421; paired over 30 tissues TAI−LOH = +0.116 [0.089, 0.145], 29/30 tissues
favour TAI).

The C1 offset correlates with each component's tissue **mean** (r −0.47 to −0.56)
but with **no** component **share** (all CIs span zero; adding shares to a
level-only model gives F-test p = 0.61). Tissue scar *mix* is also nearly constant
(share_LOH IQR 0.263–0.301). So C1 is regression toward the grand mean on tissue
**level** plus model-side lineage imprinting — **not** tissues carrying different
*kinds* of scar. With n=30 tissues this is "no evidence for the mix hypothesis",
not proof of absence.

---

## 11. Candidate selection

Pre-registered in `docs/28` §4 (commit `28bb913`) while `results/v3_loco_full/`
was verifiably empty. V3-abs passed all seven criteria S1–S7 and was designated
**primary**; run 01 became the pre-declared sensitivity analysis.

The counter-argument was recorded and declined on the record: Candidate A is
independently replicated in `pipeline_glmnet/` and the ranking comparison was a
statistical tie. The rule was followed rather than overridden.

---

## 12. Frozen artifacts

| | Primary — V3-abs | Sensitivity — run 01 |
|---|---|---|
| Path | `results/frozen_v3_2026-09-18/frozen_nonCNS.rds` | `results/frozen_2026-09-18/frozen_nonCNS.rds` |
| sha256 | `df9f7e82…` | `0e975cd6…` |
| feature_rank | `within_tissue` | `pooled` |
| alpha / lambda | 0.1 / 1.3598 | 0.1 / 1.4201 |
| Non-zero coefficients | 905 | 774 |
| Training / calibration | 5,664 / 1,401 | 5,664 / 1,401 |
| Conformal q (95%) | **21.10** | **21.73** |

Both were frozen from commit `28bb913` before the primary designation was made,
and the evaluation job verified both checksums with `sha256sum -c` before
scoring.

**The conformal interval is ±21 HRD units — wider than the label's own IQR of
4–28.** Single-sample absolute precision is poor, and any presentation of this
model must say so.

---

## 13. CNS external test

LSF job `323264626`, 2026-09-18. Ledger `results/LOCKED_EVALUATION_LEDGER.tsv`:
one session, two rows (one per candidate, permitted by `docs/28` §5). No
refitting, no recalibration, no preprocessing change, no post-hoc clipping change.

### Primary — V3-abs

| | n | Pearson | Spearman | MAE | bias | slope |
|---|---|---|---|---|---|---|
| GBM | 135 | 0.320 [0.173, 0.461] | 0.279 | 8.86 | **+7.99** | 0.358 |
| LGG | 507 | 0.329 [0.234, 0.421] | 0.297 | 5.13 | +1.62 | 0.578 |
| pooled | 642 | 0.261 [0.183, 0.335] | 0.242 | 5.92 | +2.96 | 0.358 |

### Sensitivity — Candidate A

GBM r = 0.293 [0.142, 0.446], MAE 9.79, bias +9.11; LGG r = 0.371 [0.278, 0.460],
MAE 5.95, bias +3.69. **CIs overlap heavily; neither candidate is clearly better
on CNS.** Reported because `docs/28` §5 clause 3 forbids reporting only the
better one.

### Null comparison

The oracle tissue-mean null — which uses CNS labels and is therefore unavailable
at deployment — **beats the model in every group**: skill −1.009 (GBM), −0.048
(LGG), −0.233 (pooled).

### Placement among the 30 source tissues

The primary interpretive frame: *would CNS look unusual as the 31st held-out
tissue?*

| Statistic | GBM | LGG |
|---|---|---|
| Pearson | 13.3rd pct (28th of 32) | 16.7th pct (26th of 32) |
| Spearman | 16.7th pct (27th of 32) | 16.7th pct (26th of 32) |
| MAE | 43.3rd pct (19th of 32) | 83.3rd pct (6th of 32) |
| bias | 10.0th pct (29th of 32) | 66.7th pct (11th of 32) |
| Calibration slope | 16.7th pct (27th of 32) | 36.7th pct (20th of 32) |

**Yes, CNS would have looked unusual.** It is a poor relation, not a typical
held-out tissue. LGG's favourable MAE placement is best read as a consequence of
its low, tight HRDsum distribution rather than good calibration — a judgement,
not a computed decomposition.

### Two defects the CNS run exposed

- **OOD flagged only 8 of 642** samples (all LGG). The frozen reportability rule
  did **not** warn that an entire unseen lineage was out of distribution. This is
  a deployment-blocking defect (**B17**), and it cannot be fixed by tuning against
  the now-consumed CNS labels.
- **The HRDsum ≥ 42 cutoff is vacuous in CNS**: GBM has **0** samples above it and
  LGG exactly **1**. Any AUC at that threshold is degenerate and no binary claim is
  possible in this lineage (**B18**).

---

## 14. Interpretation

**Rank transfer is real but weak.** Both lineages show correlations excluding
zero. Methylation carries HRD-associated signal that survives a lineage shift.
But at the 13th–17th percentile of source tissues, "survives" is the right verb,
not "generalises".

**Calibration transfer failed outright.** A slope of 0.358 in GBM against an ideal
of 1.0, a bias of +8 HRD units, and a model that loses to its own tissue's mean.
C1 is not a residual nuisance; out of distribution it is the dominant term.

**The finding that should shape the next model:** V3-abs measurably reduced
lineage imprinting *within* the source distribution — 42% of the excess, tightly
estimated — and that produced **no** calibration benefit on an unseen lineage.
Reducing lineage structure among tissues you have is not sufficient for
transferring to a tissue you do not have. Whatever mechanism carries the offset
into a new lineage is not the one V3-abs removed.

---

## 15. Limitations

- Development metrics are a *development* result: clipping, the log1p rejection
  and the entire C1 investigation were all developed against those same folds.
- The macro estimator used for A-vs-V3 comparison is **not** interchangeable with
  the pooled-within estimator behind run 01's 0.612 headline.
- Clipping is not rank-neutral *within* tissue even though pooled
  Spearman(raw, clipped) = 1.000; in KICH (31/65 predictions negative) it moves
  within-tissue Pearson 0.803 → 0.888.
- n=30 tissues gives limited power for paired comparisons; |r| > ~0.36 is needed
  for p < 0.05.
- OV n=10 — the canonical HRD cancer is effectively missing.
- The lambda grid is not strictly nested (§5).
- CNS is adult GBM/LGG. **Pediatric high-grade glioma is a different disease**,
  and nothing here validates pediatric transfer.
- The collaborator's conflicting "V3" report was never reconciled (`docs/28` §2).

---

## 16. Claims supported

- The model contains transferable, HRD-associated methylation signal in adult
  TCGA source cancers: within-tissue ranking positive in 29/30 tissues,
  permutation p = 0.001.
- That signal is **not purely lineage recognition**: within-tissue top-quartile
  AUC 0.777–0.779 with the between-tissue axis removed by construction, and the
  signal survives purity adjustment (partial correlation *rose* 0.612 → 0.621).
- Ranking by within-tissue variance removes **42%** (CI 37.9–46.5%) of excess
  lineage structure at no measurable ranking cost.
- Weak but non-zero rank transfer into an unseen lineage: GBM r = 0.320, LGG
  r = 0.329, both CIs excluding zero.
- Clipping at zero is defensible and resolves physically invalid negatives.
- Frozen inference works technically, end to end, with verified provenance.

## 17. Claims explicitly NOT supported

- **Not** an absolutely calibrated zero-shot HRDsum predictor in an unseen tissue.
  GBM bias +7.99, slope 0.358, and the model loses to the tissue-mean null.
- **Not** a validated clinical HRD assay, and not a treatment-selection or PARP
  inhibitor response result.
- **Not** an HRD probability. The output is a continuous regression estimate.
- **No** supported binary HRD-high claim in CNS — the ≥42 cutoff is vacuous there.
- **Not** pediatric-validated. PBTP remains untouched.
- **Not** demonstrated that reducing in-distribution lineage imprinting improves
  out-of-distribution calibration — this analysis is evidence *against* it.
- The pooled AUC of 0.862–0.866 does **not** represent clinical discrimination;
  most of it is lineage recognition (tissue-identity control 0.723–0.743).

---

## 18. Next steps

1. **Do not re-score CNS.** It is consumed. Any successor model is validated on
   PBTP or on nothing.
2. **Fix B17 first.** An OOD rule that fails to flag an entire unseen lineage is
   deployment-blocking, and fixing it needs no new labels.
3. **Attack the real question:** why did in-distribution lineage reduction not
   yield out-of-distribution calibration transfer? Candidates — hierarchical /
   multitask elastic net, anchor regression, group DRO, tissue-subspace removal,
   adversarial tissue-invariant representations, hurdle models for the
   zero-inflated low-HRD regime.
4. **Budget for labels.** The EB shrinkage few-shot path is the only approach with
   demonstrated traction (k=10 → 58% of oracle gain). For pediatric deployment,
   plan for ~10 labelled PBTP cases rather than hoping for zero-shot transfer.
5. **Protect PBTP.** It is the last genuinely independent test set this project
   has, and its value is destroyed by exactly one premature look.

# 24 — Results: LOCO run 01 (2026-09-17)

> **Status note added 2026-09-18.** This document is the historical record of
> LOCO run 01 and is unchanged below. Run 01 has since been **superseded as the
> shipped model by V3-abs**, and the locked CNS partition described as "still
> unopened" below **was opened on 2026-09-18**. See the UPDATE section at the
> end of this file before quoting anything here.

**Run.** LSF array `323078995`, 30/30 folds DONE, merged 06:20 CDT by
`scripts/watch_loco_array.sh`. Development cohort n = 7,065 across 30 non-CNS
cancer types. Locked CNS partition (GBM+LGG, 642 samples) **still unopened.**

**Verdict.** The tissue confound is **real but not disqualifying.** Modelling
HRDsum continues, with a changed claim: this is a **within-tissue relative
ranker**, not an absolute HRD calculator. See §6 for what must change.

---

## 1. Headline numbers

| Quantity | Value | File |
|---|---|---|
| Pooled MAE (model) | **9.049** | `pooled_null_panel.tsv` |
| Pooled MAE (tissue-mean null) | 9.890 | `pooled_null_panel.tsv` |
| **Skill vs tissue-mean null** | **+0.085** | `pooled_null_panel.tsv` |
| Skill vs training-mean null | +0.328 | `pooled_null_panel.tsv` |
| Macro MAE (unweighted over tissues) | 8.566 | `macro_metrics.txt` |
| **Within-tissue Pearson** | **0.612** | `pooled_within_tissue_metrics.tsv` |
| Within-tissue Spearman | 0.601 | `pooled_within_tissue_metrics.tsv` |
| Within-tissue MAE | 7.881 | `pooled_within_tissue_metrics.tsv` |
| Within-tissue permutation p | **0.001** (1000 perms) | `pooled_within_tissue_permutation.tsv` |
| Lambda at path boundary | **0 of 30 folds** | `loco_metrics.tsv` |
| Reportable rate | 99.2% | `loco_predictions.tsv` |

---

## 2. The apparent contradiction, and its resolution

Two numbers look incompatible:

- Skill vs tissue-mean null = **+0.085** — barely better than guessing each
  cancer type's average.
- Within-tissue correlation = **0.612**, permutation p = 0.001 — strongly
  better than chance at ranking patients inside a tissue.

Both are correct. They disagree because the model gets the **ranking** right and
the **absolute level** wrong.

Per-tissue calibration offset (`bias = mean(predicted − actual)`):

| Statistic | Value |
|---|---|
| Mean absolute offset | **3.46 HRD units** |
| SD of offset across tissues | 4.83 |
| Range | −15.44 (KIRC-like low) to +8.29 (PCPG) |
| Prediction shrinkage (mean sd_pred/sd_actual) | 0.772 |

Remove that per-tissue offset and pooled MAE falls **9.049 → 7.881**. The offset
alone costs ~1.17 MAE units — more than the entire margin over the tissue null.
So the model *has* learned HRD-relevant signal; it is squandered by systematic
mis-levelling of each tissue.

`results/loco_run01/figures/within_vs_absolute.png` shows this directly.

---

## 3. How much of the model is just tissue identity?

| Test | Result | Reading |
|---|---|---|
| Variance of **predictions** explained by tissue identity alone | R² = **0.562** | 56% of what the model says is recoverable from the label "this is a breast tumour" |
| Variance **not** explained by tissue | **0.438** | 44% is something else |
| Variance of **true HRDsum** explained by tissue identity | R² = **0.341** | The confound is real in the biology, not invented by the model |
| Within-tissue permutation test | p = 0.001 | Residual signal survives after tissue is removed |

The model leans on tissue harder than the truth warrants (0.562 vs 0.341) — it
**over-weights lineage**. But 44% of its behaviour is not tissue, and that
portion is what produces r = 0.61 within tissue. A pure tissue-lookup model
would show within-tissue r ≈ 0 and permutation p ≈ 0.5. It does not.

---

## 4. Per-tissue breakdown

29 of 30 tissues have within-tissue r > 0; 26 of 30 exceed r > 0.3; median
r = **0.551**. Full table: `per_tissue_within_correlation.tsv`.

**Best:** KICH 0.803, UCEC 0.758, STAD 0.738, LUAD 0.721, BLCA 0.675, PRAD 0.663.

**Failures:** THCA **−0.032** (bias +8.05), PCPG 0.110 (bias +8.29),
CHOL 0.161 (n=35), TGCT 0.162.

The failure pattern is coherent: THCA and PCPG are genomically quiet, low-HRD
tumours, and the model **over-predicts both by ~8 units**. Trained mostly on
scarred tumours, it cannot express "this genome is calm." Combined with the
zero-floor problem (§5) this is the single most actionable defect.

---

## 5. Limitations that constrain interpretation

**Ovarian is effectively missing — n = 10.** OV is the canonical HRD cancer and
the primary clinical use case for HRD testing. TCGA ovarian methylation is
predominantly 27k-array, so it was excluded by the 450k probe bridge. Any claim
about HRD-high tumours rests on UCEC/BRCA/STAD instead. **This must be stated in
the presentation.** Similarly thin: CHOL 35, DLBC 47, UCS 56.

**HRDsum is bounded at zero and 14.3% of samples sit exactly at 0.** The elastic
net is unbounded and emits negative predictions (visible in the figure). MAE is
inflated by predictions that are impossible a priori. Clipping at 0, or modelling
`log1p(HRDsum)`, is a free improvement not yet applied.

**Purity is a live, unresolved concern.**

| Quantity | Value |
|---|---|
| cor(prediction, purity) within tissue | **0.165** |
| cor(true HRDsum, purity) within tissue | **0.019** |

Predictions track tumour purity roughly **8× more strongly than the truth does**.
Some of what the model reads is "how much tumour is in this sample," not HRD.
Skill degrades monotonically as purity rises:

| Purity stratum | n | MAE model | MAE tissue null | Skill |
|---|---|---|---|---|
| Low (0.08–0.51) | 2351 | 8.599 | 10.611 | **+0.190** |
| Mid (0.52–0.73) | 2344 | 9.184 | 9.726 | +0.056 |
| High (0.74–1.00) | 2202 | 9.484 | 9.112 | **−0.041** |

Note the direction: model MAE *worsens* with purity (8.60→9.48) while the tissue
null *improves* (10.61→9.11). In the cleanest, highest-purity samples the model
**loses to the tissue mean.** That is the opposite of what a genuine biological
signal should do and is not yet explained. Purity-matched subset
(0.5–0.8, n=3312) retains skill +0.045, so the effect is attenuated but not
abolished.

**No CNS validation yet.** The locked partition is untouched. Everything above
is development-cohort performance and does not establish transfer.

---

## 6. Decision and required changes

**Continue modelling HRDsum.** Justification: within-tissue r = 0.61 with
permutation p = 0.001 across 29/30 tissues is not a tissue-lookup artefact, and
0/30 folds hit a lambda boundary so the fit is not grid-limited.

Required before any result is presented or the CNS lock is opened:

1. **Reframe the claim.** Report as a within-tissue relative ranker. Lead with
   within-tissue r = 0.61 and permutation p = 0.001, not with skill = +0.085.
   Do **not** describe the output as an HRD score comparable across tissues.
2. **Fix per-tissue calibration.** Worth ~1.17 MAE units. Must be fitted inside
   the LOCO loop on training tissues only — a per-tissue offset fitted on the
   held-out tissue is leakage and invalidates the fold.
3. **Enforce the zero floor.** Clip at 0 or model `log1p(HRDsum)`.
4. **Resolve the purity inversion** (§5). Until explained, the high-purity
   result blocks any strong biological claim.
5. **State the OV limitation** in every presentation of these numbers.
6. **Keep the CNS lock closed** until 1–4 are done. It opens once.

---

## 7. Provenance

All numbers from `results/loco_run01/`, array `323078995`, commit `696e05e`
plus this run's outputs. Matrix `data/processed/beta.tsv`
(sha256 `e3642d30…`, `engineering_only: false`), 336,480 probes × 7,707 samples.
Per-fold artefacts in `results/loco_run01/folds/` (30 × 5 files).

---

# UPDATE 2026-09-18 — run 01 is SUPERSEDED as the shipped model

**Everything above §7 is unchanged and remains the correct record of LOCO run
01.** Nothing in it has been retracted. This section records only what happened
afterwards, so a reader arriving at this file does not mistake run 01 for the
project's final model.

## S1. Run 01 is now Candidate A, the sensitivity analysis

The shipped model is **V3-abs**: this same absolute-target LOCO model with
**exactly one change** — the 5,000-probe unsupervised filter ranks probes by
pooled **within-tissue** variance computed on training-fold rows only, rather
than pooled total variance. Pre-registered in
`docs/27_V3_PREREGISTRATION.md` §3 before any V3 output existed; it passed 7 of
7 sentinel gates (job `323220006`) and 7 of 7 full-30-fold selection criteria
S1–S7 (array `323242800`, scored under the rule frozen in
`docs/28_CANDIDATE_SELECTION_AND_CNS_PROTOCOL.md` §4).

That rule was written against §3 of this document. §3 asked how much of the
model is just tissue identity and answered R² = 0.562 of the prediction against
0.341 of the truth. V3-abs attacks precisely that excess.

## S2. What changed, and what did not

Head-to-head on all 30 development folds: **`results/tables/tableA2_A_vs_V3_full30.md`**
(macro estimator under adopted C2 clipping, seed 5813, 10,000 bootstrap
replicates). Note that macro estimator is *not* the pooled-within estimator that
produced the headline 0.612 in §1; the two must not be quoted side by side.

| Quantity | Run 01 (A) | V3-abs | Reading |
|---|---|---|---|
| Tissue R² of **predictions** | 0.5584 | **0.4672** | against 0.3414 for the truth |
| Excess tissue R² (pred − truth) | 0.2171 | **0.1258** | `ΔR²_excess = 0.0913` [0.0824, 0.0996], p < 0.0001 |
| Fraction of excess lineage structure removed | — | **42.0%** | CI 37.9–46.5% |
| Mean absolute tissue bias (C1) | 3.4628 | **2.9737** | §2's 3.46 units, reduced not resolved |
| Pooled MAE | 8.9717 | **8.4456** | — |
| Macro skill vs oracle tissue-mean null | 0.4192 | **0.6304** | — |
| Macro within-tissue Pearson | 0.5197 | 0.5036 | paired mean −0.0161, CI [−0.0410, 0.0051], **p = 0.428** |
| Macro within-tissue Spearman | 0.4747 | 0.4472 | paired mean −0.0275, **p = 0.477** |
| Tissues with positive within-tissue r | 29 of 30 | 29 of 30 | — |

Only 12 of 30 tissues improved on Pearson and 14 of 30 on Spearman, and the
paired difference is not significant either way. **V3-abs is a targeted lineage
reduction at no measurable ranking cost — it is not a general accuracy
improvement**, and §6's required change 1 (report as a within-tissue relative
ranker) still stands unaltered.

Shipped artifact: `results/frozen_v3_2026-09-18/frozen_nonCNS.rds`, sha256
`df9f7e82b0b371b81ecca6a1d99a1a50128a95e78daa63307d5f3b0633a4aa72`, alpha 0.1,
lambda 1.3598, 905 non-zero coefficients, 95% conformal q = 21.10. Run 01 was
frozen alongside it as the pre-declared sensitivity candidate,
`results/frozen_2026-09-18/frozen_nonCNS.rds`, sha256 `0e975cd6…`, alpha 0.1,
lambda 1.4201, 774 non-zero, conformal q = 21.73.

## S3. The §6 gate logic was corrected, not satisfied

§6 item 6 said "keep the CNS lock closed until 1–4 are done." `docs/27` §2
records that this was **circular**: it required zero-shot calibration to be
proven before running the only experiment able to measure zero-shot calibration
in an untouched lineage. The gate was split into a DEPLOYMENT gate (still shut)
and an EXTERNAL-EVALUATION gate (satisfied 2026-09-18). C1 was reclassified from
a precondition into a finding to be tested on CNS. C3, the purity inversion of
§5, remains unresolved.

## S4. The CNS result

The locked partition was opened **exactly once**, 2026-09-18, LSF `323264626`,
642 samples (135 GBM, 507 LGG), ledger row in
`results/LOCKED_EVALUATION_LEDGER.tsv`. Primary (V3-abs): GBM r = 0.320
[0.173, 0.461], MAE 8.86, bias +7.99; LGG r = 0.329 [0.234, 0.421], MAE 5.13,
bias +1.62. The **oracle** within-lineage tissue-mean null beats the model in
both (skill −1.009 GBM, −0.048 LGG). Rank transfers weakly; absolute
cross-lineage calibration does not. **The 42% lineage reduction measured on the
source tissues did not translate into better CNS calibration.**

Full result: `results/cns_eval_2026-09-18/analysis_v3/cns_analysis_summary.md`
(primary) and `analysis_A/` (sensitivity); protocol `docs/28` §6–7; pre-unlock
audit `docs/29_PRE_UNLOCK_RECORD.md`.

## S5. Also recorded: V4 failed

A fourth variant, V4 lineage-penalized, was sentinelled (job `323241956`) and
**failed the mechanism gate G3**: tissue R² of the predictions rose
0.5219 → 0.5361 against a ceiling of 0.4919, while passing the other six gates.
The branch was terminated. Scorecard: `results/v4_sentinel/gate_scorecard.tsv`.


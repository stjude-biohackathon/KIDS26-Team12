# 28 — Candidate selection rule and CNS evaluation protocol

**Written 2026-09-18 11:58 CDT, BEFORE the full 30-fold V3-abs array (LSF
`323242800`, submitted 11:55) returned a single fold.** `results/v3_loco_full/`
was empty at the time of writing. Every threshold and every metric below is
fixed now, while the inputs are still unknown.

This document governs (a) which model is frozen as the final source-only
candidate, and (b) exactly what happens when the locked CNS cohort is opened.

---

## 1. Why this is written now

Choosing between run 01 and V3-abs after seeing both 30-fold results, with no
rule written down first, is model selection by inspection. It would not
invalidate the CNS test — CNS is untouched either way — but it would make the
source-side claim unfalsifiable, and this project has an explicit standard
against gates that cannot fail (`docs/21` B5).

So the rule goes first.

---

## 2. State of the evidence at time of writing

### V3-abs sentinel result (scored 11:53, gates frozen at 10:27 in `docs/27` §5)

All 7 gates passed. Scorecard in `results/v3_sentinel/gate_scorecard.tsv`:

| Gate | Threshold | Baseline | Observed | Pass |
|---|---|---|---|---|
| G1 macro within-r ex-THCA | ≥ 0.7096 | 0.7596 | **0.7688** | ✓ |
| G2 BRCA within-r | ≥ 0.5823 | 0.6323 | **0.6729** | ✓ |
| G3 tissue R² of predictions | ≤ 0.4919 | 0.5219 | **0.2991** | ✓ |
| G4 pooled MAE | ≤ 12.143 | 11.039 | **8.6137** | ✓ |
| G5 mean abs tissue bias | ≤ 5.674 | 5.1581 | **3.2029** | ✓ |
| G6 non-THCA tissues not degraded | ≥ 2 of 3 | 3 of 3 | **3 of 3** | ✓ |
| G7 provenance recorded | all TRUE | — | TRUE | ✓ |

The mechanism gate G3 is the one that matters. Tissue identity explained 52.2%
of the run-01 prediction on these four tissues and 29.9% of the V3-abs
prediction, against 29.4% for the observed HRDsum. **The excess lineage
imprinting that `docs/26` §8 diagnosed as the cause of C1 is largely gone on the
sentinel tissues, and the per-tissue offsets shrank with it** (BRCA +4.96 →
−2.27, UCEC +6.49 → +3.23, mean |bias| 5.16 → 3.20) while within-tissue ranking
improved rather than degraded.

This is the first intervention in the project to move C1 at all.

### A conflicting external report, and why the two are not comparable

A collaborator reports having run "the V3 sentinel" overnight with no
improvement over the existing models. That is most likely **docs/26's V3**,
defined in `docs/26` §4 as "**V2** with pooled within-tissue variance feature
ranking" — i.e. the within-tissue filter stacked on top of V2's relative
target, which independently failed 3 of 7 gates in sentinel `323195423`.

**V3-abs is a different experiment**: run 01's absolute-target model with the
feature filter changed and *nothing else*. `docs/27` §3 records this distinction
in writing at 10:26, before either result existed, precisely because stacking a
new factor on a refuted one is uninterpretable.

The two reports are therefore not in conflict; they are different models. The
collaborator's result, if it is docs/26-V3, is additional evidence that V2's
relative target is the problem. **This must be checked before the write-up
asserts anything about it** — if the collaborator did run the single-factor
variant, the discrepancy is real and matters, and my full 30-fold run is the
tiebreak. Recorded here as an open question, not a resolved one.

### What is NOT yet known

The sentinel is four tissues. `docs/21` B16 exists because a 4-fold log1p
preview inverted the verdict of the full 30-fold run. **A sentinel may stop a
run; it may not authorise a conclusion.** The full 30-fold V3-abs array is the
evidence that counts, and it is still running.

---

## 3. Run-01 reference values, all 30 development tissues

Computed 11:58 CDT under the adopted C2 clipping, from
`results/loco_run01/loco_predictions.tsv`:

| Statistic | Run 01 |
|---|---|
| Macro within-tissue Pearson | 0.5197 |
| Macro within-tissue Spearman | 0.4747 |
| Pooled MAE | 8.9717 |
| Mean absolute tissue bias (C1) | 3.4628 |
| **Tissue R² of predictions** | **0.5584** |
| Tissue R² of observed HRDsum | 0.3414 |
| Tissues with positive within-tissue r | 29 of 30 |

Note these are computed per-tissue then macro-averaged under clipping, which is
not identical to the pooled-within estimator that produced the headline 0.612 in
`docs/24`. Both are legitimate; they are different estimators and must not be
quoted side by side as if comparable. Comparisons below use the macro estimator
consistently for both candidates.

---

## 4. THE SELECTION RULE (frozen)

Two candidates, both to be frozen as inference artifacts before CNS opens:

- **Candidate A — run 01 / pooled.** `feature_rank="pooled"`, freeze job
  `323234897` → `results/frozen_2026-09-18/`.
- **Candidate B — V3-abs.** `feature_rank="within_tissue"`, freeze job
  `323242801` → `results/frozen_v3_2026-09-18/`.

**Candidate B is selected as the primary frozen candidate if and only if ALL of
the following hold on the complete 30-fold LOCO run:**

| # | Criterion | Threshold |
|---|---|---|
| **S1** | Macro within-tissue Pearson | ≥ 0.4697 (A − 0.05) |
| **S2** | Tissue R² of predictions | ≤ 0.5084 (A − 0.05), **and** strictly below A |
| **S3** | Mean absolute tissue bias | ≤ 3.4628 (no worse than A) |
| **S4** | Pooled MAE | ≤ 9.8689 (A × 1.10) |
| **S5** | Tissues with positive within-tissue r | ≥ 27 of 30 |
| **S6** | Breadth, not one-tissue artefact | B's within-r ≥ A's − 0.05 in ≥ 20 of 30 tissues |
| **S7** | All 30 folds completed; every `metrics_*.tsv` records `feature_rank=within_tissue`; all 8 test scripts pass | — |

**If any of S1–S7 fails, Candidate A (run 01) is frozen as primary** and V3-abs
is reported as a source-only mechanism result that did not survive the full run.

Rationale for the asymmetry: A is the incumbent with a complete, already-merged
30-fold run and an independent replication of its architecture
(`pipeline_glmnet/`). B must earn the displacement. S2 is the reason B exists at
all — if B does not reduce lineage imprinting across all 30 tissues, it has not
reproduced its own sentinel and should not ship, no matter how good its MAE is.

S6 guards against the failure mode `docs/21` warns about repeatedly: a macro
average improved by two or three tissues while most get worse.

### Tie-breaking and simplicity

If B passes every criterion but the margins on S1/S2 are within ±0.01 of A —
i.e. the two are performing equivalently — **A is retained**, on the grounds
that an unsupervised pooled-variance filter is simpler, already replicated, and
already documented. Displacement requires a real improvement, not noise.

---

## 5. Multiplicity: how many candidates may be scored on CNS

**The locked CNS cohort will be scored with the primary candidate. That single
result is the project's headline external test.**

The secondary candidate may *also* be scored, in the same session, but only
under conditions that are fixed here, before any CNS label is visible:

1. Both artifacts are frozen and checksummed **before** the first CNS label is
   read. No artifact may be built or modified after unlock.
2. The primary candidate is designated **before** unlock, by §4, and recorded in
   the pre-unlock audit file.
3. Both results are reported. **Reporting only the better one is forbidden**,
   and is the specific abuse this clause exists to prevent.
4. The secondary result is labelled a pre-declared sensitivity analysis. It may
   not be promoted to headline after the fact, whatever it shows.

Scoring both is legitimate here only because the choice between them is made on
source data by §4 and written down first. What would *not* be legitimate — and
is prohibited — is scoring one, disliking the number, and swapping. If the team
prefers maximum conservatism, scoring the primary alone is always the safe
option and costs nothing scientifically.

**Under no circumstances may a third model be built after unlock.**

---

## 6. CNS evaluation protocol (frozen)

Predeclared, to be executed exactly once by `scripts/predict_frozen.R`.

**Cohort.** 642 locked samples, `partition == "locked_CNS"`, GBM and LGG.
Reported **separately first**, pooled CNS second.

**No refitting, no recalibration, no preprocessing changes, no CpG changes, no
alpha/lambda changes, no post-hoc clipping changes.** Raw predictions are
preserved and hashed before any plotting.

**Metrics, fixed now:**

- *Ranking*: Pearson, Spearman, with bootstrap CIs (10,000 resamples, seed
  recorded). This is the primary endpoint.
- *Absolute*: MAE, RMSE, median absolute error.
- *Calibration*: mean signed error (bias), calibration intercept and slope from
  `lm(observed ~ predicted)`, observed vs predicted mean and SD, range
  compression.
- *Null*: the tissue-mean null, computed **within** GBM and within LGG. Note
  honestly that this null uses CNS labels and so is an *oracle* baseline — it is
  not available at deployment for a genuinely new tissue. It answers "does the
  model beat knowing the tissue's mean?", which is the question that matters,
  but it flatters nothing and must be labelled as an oracle.
- *OOD*: the frozen `ood_score` / `reportable` flags. Report n reportable, n
  flagged, score distribution, and whether accuracy differs between groups.
- *Purity*: secondary residual diagnostic if `purity` is populated for CNS.

**Threshold analysis is exploratory only.** HRDsum ≥ 42 may be shown, labelled
"ranking performance at an exploratory cutoff". It is not a primary endpoint and
no continuous prediction will be presented as a probability.

**The primary interpretive frame** is §15 of the handoff: place GBM and LGG
inside the empirical distribution of the 30 source LOCO tissues, and ask *"would
CNS look unusual if it had simply been the 31st held-out tissue?"* That is more
informative than any pass/fail line, and no pass/fail line is being drawn.

**The CNS result cannot change the model.** If transfer is poor, that is the
finding, and the next-generation hypothesis gets documented rather than
implemented.

### Few-shot calibration on CNS

Permitted **only** as an explicitly separate, clearly-labelled secondary
analysis, using the already-developed EB shrinkage (τ²=22.3, σ²=108.6) at
k ∈ {1,3,5,10}, with calibration subsets drawn by a seeded random rule fixed
here: `set.seed(4127)`, 200 random draws per k, disjoint from the evaluation of
the remaining samples. It must never be described as zero-shot, and the
zero-shot numbers must be reported alongside it.

---

## 7. What is claimed either way

Fixed in advance so the conclusion is not written to fit the number:

- **If rank transfers and calibration does not** — the expected outcome given
  C1 — the claim is: *methylation carries HRD-associated signal that transfers
  in rank to an unseen lineage, but absolute genomic-scar calibration remains
  lineage-dependent.* Deployment gate stays shut.
- **If both transfer**, the claim is: *evidence of ranking and quantitative
  transfer into adult CNS tumours, warranting independent pediatric testing.*
  Still not a clinical claim, still not validated, still needs PBTP.
- **If neither transfers**, the claim is: *the pan-cancer methylation–HRD
  association did not generalise to the held-out CNS lineage under a strict
  frozen zero-shot test.* This is a real and publishable result.

In all three cases the deployment gate (`docs/27` §2A) remains shut. Nothing
today produces a validated clinical HRD assay.

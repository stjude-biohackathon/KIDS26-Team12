# 27 — V3-abs pre-registration, and the gate-logic correction

**Written 2026-09-18 10:26 CDT. LSF job `323220006` submitted 10:20 CDT and was
RUNning at the time of writing; `results/v3_sentinel/folds/` was verified EMPTY
(`ls -la`, 10:26 CDT) immediately before this file was committed.**

Nothing in this document was written with knowledge of a V3 result. That is the
entire point of the file, and it is the reason it is committed as its own commit
before the sentinel finishes. If you are reading this in the git history, check
that its commit timestamp precedes the first file in `results/v3_sentinel/`.

---

## 1. Why this document exists

`docs/21` B16 records the rule this project operates under:

> A sentinel may **stop** a full run; it may not on its own **authorise** one
> unless the pre-registered gates in `docs/26_C1_N_OF_1_WORKAROUND.md` §7 are met.

B16 exists because a 4-fold log1p preview inverted the verdict of the full
30-fold run (`docs/21` C2). Sentinels are not representative by construction.
So the gates have to be numeric, they have to be written down first, and they
have to be scored honestly afterwards even when the answer is inconvenient.

---

## 2. The gate-logic correction (supersedes the C1–C4 gate in docs/21, 23, 24)

Two different questions have been conflated throughout this project under the
single phrase "the CNS lock". They are now separated.

### Gate A — the DEPLOYMENT gate

> Can this model be claimed to produce trustworthy **absolute** zero-shot HRDsum
> estimates in a tissue it has never seen?

**Status: NO, and nothing today is expected to change that.** C1 — a per-tissue
additive offset of mean magnitude 3.46 HRD units, range SARC −15.44 to PCPG
+8.29 — is unresolved after three independent failed attacks (`docs/26` §8). The
95% predictive interval for an unseen tissue's offset is −9.7 to +8.1 units,
74% as wide as the label's own IQR. This gate stays shut. No presentation,
figure, caption or README line may imply otherwise.

### Gate B — the EXTERNAL-EVALUATION gate

> Has the architecture, preprocessing, model-selection rule, evaluation plan and
> interpretation been frozen firmly enough that the locked CNS cohort can be
> opened **exactly once** as an external domain test?

**This is the gate being satisfied today.**

The correction: docs/21, docs/23 and docs/24 all previously gated the CNS unlock
on *resolving* C1. That is circular. It requires zero-shot calibration to be
proven before running the only experiment capable of measuring zero-shot
calibration in a genuinely untouched lineage. C1 is a **finding to be tested on
CNS**, not a precondition for testing.

The CNS evaluation is therefore explicitly permitted to **fail**. A result
showing that rank transfers but calibration does not is a complete, publishable,
honest answer to the project's scientific question. It is not a project failure.

What Gate B actually requires — all of which must be true *before* any CNS label
is read:

1. A single final source-only candidate selected under a decision rule written
   before its inputs were known (this document).
2. That candidate frozen as an inference artifact with a recorded checksum.
3. The frozen inference path validated on non-CNS data.
4. The CNS metrics, and the script that computes them, predeclared.
5. An auditable pre-unlock record: git commit, candidate, rationale, artifact
   checksum, timestamp, and an explicit statement that no CNS outcome label was
   inspected during candidate selection.

**The CNS result may not be used to choose among candidate models.** Once the
labels are open, the candidate is fixed forever, whatever the numbers say.

---

## 3. What V3-abs is, and why it is not the V3 in docs/26

`docs/26` §4 defines V3 as "**V2** with pooled *within-tissue* variance feature
ranking", where V2 is the relative-target, tissue-weighted ranker. That
definition bundles **two** changes on top of the absolute model: a changed
prediction target *and* a changed feature filter. V2's target change already
failed its own gates (sentinel `323195423`, 3 of 7 gates failed, `docs/26` §7).
Running docs/26-V3 would therefore stack a new factor on a refuted one and
produce an uninterpretable result.

**V3-abs** is the single-factor experiment instead:

> The run-01 absolute-target LOCO model, entirely unchanged, except that the
> 5,000-probe unsupervised filter ranks probes by pooled **within-tissue**
> variance computed on training-fold rows only, instead of pooled **total**
> variance.

Everything else is held fixed: identity target transform, clipping applied at
evaluation, `alpha ∈ {0.1, 0.5, 1}`, data-derived lambda path, inner folds = one
per training cancer, `standardize=FALSE`, seed 260910, 5,000 features.

### The hypothesis being tested

`docs/26` §8 diagnosed the mechanism behind C1: the filter ranks by **pooled**
variance *before* the target is consulted. Across a pan-cancer matrix, the
probes with the largest pooled variance are overwhelmingly those separating
lineages, because between-tissue methylation differences dwarf within-tissue
ones. The model is therefore handed a lineage-discriminating basis and can do
little else but encode tissue identity — which is exactly what run 01 shows
(tissue identity explains 56.2% of the prediction, versus far less of the truth).

Within-tissue variance ranking removes the between-tissue component from the
ranking statistic:

$$ \mathrm{Var}_{\text{within}}(j) = \frac{1}{N-T}\sum_{t}\left(Q_{t,j} - S_{t,j}^2/n_t\right) $$

selecting probes that vary *inside* tumour types. If C1's mechanism is correctly
diagnosed, this should reduce lineage imprinting while preserving — possibly
improving — within-tissue HRD ranking.

**Leakage requirement.** The statistic is recomputed inside every training fold,
outer and inner, on that fold's rows only. Enforced at `R/model.R:275` with call
sites at `:516` (outer refit) and `:541` (inner fold), and tested in
`tests/test_v3_feature_rank.R` check 2, which replaces held-out feature rows and
labels with garbage and asserts the preprocessing constants, selected features,
coefficients, alpha and lambda are `identical()`.

---

## 4. Sentinel design

**Tissues: BRCA (3), KICH (10), THCA (26), UCEC (28).** These are *not* newly
chosen. They are the same four pre-registered in `docs/26` §7c for the C1c
sentinel, reused deliberately and without modification so the two experiments
are directly comparable and so the choice cannot be accused of being tuned to
V3. Their documented rationale, verbatim from `docs/26` §7c:

> BRCA as the largest tissue, KICH as the strongest run-01 within-tissue signal
> in a small cohort, UCEC as high-signal with the largest positive offset, THCA
> as a **negative control** whose near-constant target makes rank evaluation
> degenerate.

**THCA remains a negative control and can never constitute a win.** Its
within-tissue target is near-constant (run-01 r = −0.036 after clipping), so any
rank metric there is numerically unstable and improvements are meaningless. THCA
is scored only as a *harm detector*: it is reported, and it may cause a gate to
fail, but it is excluded from every macro statistic used to pass a gate.

### Baseline for comparison

Run 01 (`results/loco_run01/loco_predictions.tsv`, array `323078995`), with the
adopted C2 clipping at 0 applied, restricted to the same four tissues, computed
at 10:26 CDT 2026-09-18 before any V3 output existed:

| Tissue | n | within r | within ρ | MAE | bias |
|---|---|---|---|---|---|
| BRCA | 743 | 0.6323 | 0.6662 | 13.177 | +4.961 |
| KICH | 65 | 0.8881 | 0.5620 | 2.102 | +1.088 |
| THCA | 464 | −0.0361 | −0.0292 | 8.237 | +8.091 |
| UCEC | 403 | 0.7584 | 0.7031 | 11.766 | +6.492 |

| Macro statistic (4 sentinel tissues) | Baseline |
|---|---|
| Macro within-tissue Pearson, all 4 | 0.5607 |
| Macro within-tissue Spearman, all 4 | 0.4755 |
| **Macro within-tissue Pearson, ex-THCA** | **0.7596** |
| **Macro within-tissue Spearman, ex-THCA** | **0.6437** |
| Pooled MAE | 11.039 |
| Mean absolute tissue bias | 5.158 |
| **Tissue R² of predictions** | **0.5219** |
| Tissue R² of observed HRDsum | 0.2940 |

That last pair is the heart of the matter. The model encodes **0.522** of its
prediction variance as lineage while the truth only carries **0.294**. The model
is substantially more lineage-determined than the thing it is predicting. V3-abs
targets precisely that excess.

---

## 5. Pre-registered gates — NUMERIC, declared before results exist

Scored on the four sentinel folds only. All thresholds are fixed now.

| # | Gate | Threshold | Rationale |
|---|---|---|---|
| **G1** | Macro within-tissue Pearson, ex-THCA | **≥ 0.7096** (baseline 0.7596 − 0.050) | Ranking is the surviving claim. A 0.05 macro loss is the largest drop defensible as noise at these n. |
| **G2** | BRCA within-tissue Pearson | **≥ 0.5823** (baseline 0.6323 − 0.050) | Largest tissue; the C1c sentinel died here (−0.405). A single-tissue floor stops a macro average hiding one collapse. |
| **G3** | Tissue R² of predictions | **≤ 0.4919** (baseline 0.5219 − 0.030) | **The mechanism gate.** V3-abs exists only to cut lineage imprinting. If it does not measurably fall, the hypothesis is wrong regardless of other metrics. |
| **G4** | Pooled MAE | **≤ 12.143** (baseline 11.039 × 1.10) | Generalisation must not materially deteriorate. |
| **G5** | Mean absolute tissue bias | **≤ 5.674** (baseline 5.158 × 1.10) | C1 must not be made worse. |
| **G6** | Non-THCA tissues not individually degraded | **≥ 2 of 3** of BRCA/KICH/UCEC lose ≤ 0.05 within-tissue Pearson | Breadth: improvement must not be carried by one tissue. |
| **G7** | Leakage and provenance tests | All 7 test scripts pass; every fold's `metrics_*.tsv` records `feature_rank=within_tissue` | Integrity. |

### Decision rule

**ALL of G1–G7 must pass to authorise the full 30-fold V3-abs LOCO array.**

G3 is non-negotiable and cannot be traded against the others. V3-abs is a
mechanism test; a version that ranks better while remaining just as
lineage-determined has not tested the hypothesis and provides no reason to
displace run 01.

### If the sentinel fails

1. **STOP V3-abs.** No re-tuning against these same four tissues — that is how
   a sentinel becomes a training set. No second V3 variant.
2. Record the failure in `docs/27` §7 and in the blocker ledger, scored gate by
   gate, including any gate that passed.
3. Proceed to **exactly one** alternative source-only experiment: supervised
   within-tissue HRD meta-association feature selection with a lineage penalty
   (§6), using the same four tissues and the same gate structure — **or**, if
   the clock does not permit it, freeze run 01 as the final candidate and
   proceed to the CNS evaluation.
4. Under no circumstances does failure here authorise an uncontrolled
   architecture search on the day of the presentation.

### Time budget

Written now so it cannot be relaxed later. A fold costs 45–85 min. The sentinel
should return by ~11:50 CDT. A full 30-fold array, submitted unthrottled to
`priority`, needs ~90–120 min and must be submitted **by 13:15 CDT** to leave
time for freeze, CNS evaluation, figures and write-up before 16:00. If the
sentinel passes but lands after 13:15, the full array is not run and run 01
remains the frozen candidate; the V3-abs sentinel is then reported as a
4-tissue mechanism probe only, never as the shipped model. **A model that has
not completed a full 30-fold LOCO cannot be frozen as the final candidate.**

---

## 6. The one alternative experiment, if V3-abs fails

Declared now so it cannot be invented to fit a disappointing result.

Within each training fold, for each probe $j$ and each training cancer $t$:
estimate the within-tissue association between methylation and HRDsum after
standardising both within tissue; Fisher-z transform; aggregate across training
cancers with a fixed-effect meta-statistic; quantify sign consistency and
between-tissue heterogeneity; and independently quantify lineage discriminability
(between/within variance ratio, the same decomposition as §3). Rank by

$$ \text{score}_j = |z_j^{\text{meta}}| \cdot c_j - \lambda \cdot \ell_j $$

for HRD consistency $c_j$ and lineage score $\ell_j$, with $\lambda$ selected in
the inner folds. Elastic net stays downstream. All statistics training-fold only;
no manual curation of probes from full-development coefficients. Same four
sentinel tissues, same gates G1–G7.

---

## 7. Scorecard

*To be completed after the sentinel returns. Section 5 above is frozen and must
not be edited when the numbers arrive.*

| Gate | Threshold | Observed | Verdict |
|---|---|---|---|
| G1 | ≥ 0.7096 | — | — |
| G2 | ≥ 0.5823 | — | — |
| G3 | ≤ 0.4919 | — | — |
| G4 | ≤ 12.143 | — | — |
| G5 | ≤ 5.674 | — | — |
| G6 | ≥ 2 of 3 | — | — |
| G7 | all pass | — | — |

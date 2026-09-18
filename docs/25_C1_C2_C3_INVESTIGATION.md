# 25 — C1/C2/C3 Investigation Results (2026-09-17)

Follow-up to `docs/24_RESULTS_LOCO_RUN01.md`. Covers the three defects that
gate opening the CNS lock. **Two of the three produced negative or surprising
results**, which changes the project's direction.

---

## C3 — The purity inversion: NOT what it looked like

### The original alarm

Run 01 reported skill falling monotonically with tumour purity
(+0.190 → +0.056 → −0.041) and `cor(pred, purity)` of 0.165 against
`cor(label, purity)` of 0.019. The natural reading was that the model had
learned to detect tumour content rather than HRD biology.

### The decisive test

If the model's signal *is* purity, then controlling for purity should collapse
the within-tissue correlation. It does not — **it slightly increases it**:

| Quantity | Value |
|---|---|
| Within-tissue `cor(pred, actual)` | 0.6124 |
| **Partial** `cor(pred, actual │ purity)` | **0.6210** |
| Change from removing purity | **+0.0087** |

A model whose apparent skill came from purity would show the opposite. Purity is
a **nuisance variable the model partially encodes, not the source of its
signal.** The 0.165 correlation is real but orthogonal to the HRD prediction.

### So what caused the inversion?

The per-tissue calibration offset (C1) interacting with a shrinking null.
Decomposed by purity tertile:

| Tertile | n | skill (raw, as run 01) | skill (offset removed) | mean abs offset |
|---|---|---|---|---|
| Low | 2361 | 0.114 | **0.242** | 3.50 |
| Mid | 2314 | 0.082 | **0.207** | 3.55 |
| High | 2222 | 0.029 | **0.159** | 3.49 |

Two things are visible. First, **the offset is essentially constant** across
tertiles (3.50 / 3.55 / 3.49) — it is not a purity effect. Second, once the
offset is removed, **skill is positive in all three tertiles**, including the
high-purity third that previously appeared negative.

`cor(per-tissue offset, per-tissue mean purity)` across 30 tissues is **−0.124**
— weak, and the wrong sign to explain an inversion.

A residual gradient remains after correction (0.242 → 0.159), so purity is not
entirely innocent. But the headline "skill goes negative" was **an artefact of a
fixed ~3.5-unit offset consuming a margin that shrinks as the tissue-mean null
gets easier to beat**, not evidence of purity capture.

### Status

**C3 DOWNGRADED** from "most serious open defect" to "secondary effect, largely
explained by C1." It is no longer the blocking scientific question.

> **Scope note.** This is a coefficient- and prediction-level analysis, chosen
> over a full probe-level test (which needs a ~4 h cluster job loading the 28 GB
> matrix). It establishes that purity does not *drive* the prediction. It does
> **not** establish which individual probes are purity-sensitive. If a reviewer
> asks "are the selected CpGs purity-associated?", that test has not been run.
> Probe selection is stable enough to make it worthwhile: 4,080 of 5,000 probes
> are shared across all 5 folds examined.

---

## C1 — Cross-tissue calibration: the covariate approach FAILS

### The design problem

A per-tissue offset is normally estimated from labelled samples of that tissue.
In LOCO — and in the pediatric application — **the held-out tissue has no labels
by definition.** A lookup table keyed by cancer type cannot transfer, and
fitting one on the held-out tissue's own labels is leakage.

The approach tested: learn a regression from **label-free tissue-level
covariates** to the tissue's offset, using the 29 training tissues as
observations, then apply it to the unseen 30th. Covariates were chosen so that
every one is computable without any HRD measurement: `mean_pred`, `sd_pred`,
`mean_purity`, `sd_purity`, `mean_ood`, `n`.

### Result: it does not work

Evaluated by leave-one-tissue-out CV, which mirrors deployment exactly. The
baseline to beat is simply predicting every tissue's offset as the global mean
offset (naive MAE = 3.466).

| Covariate set | LOTO MAE | Naive MAE | LOTO R² | Beats naive? |
|---|---|---|---|---|
| all four | 3.916 | 3.466 | −0.217 | no |
| `mean_pred` | 3.706 | 3.466 | −0.124 | no |
| `mean_purity` | 3.679 | 3.466 | −0.116 | no |
| `mean_ood` | 3.683 | 3.466 | −0.127 | no |
| `mean_pred + mean_purity` | 3.811 | 3.466 | −0.175 | no |
| `mean_pred + sd_pred` | 3.953 | 3.466 | −0.226 | no |

**Every configuration is worse than doing nothing.** Negative R² means the
calibrator generalises worse than a constant.

The marginal correlations explain why — the offset is essentially uncorrelated
with every available covariate:

| Covariate | cor with offset |
|---|---|
| `n` | +0.267 |
| `mean_purity` | −0.120 |
| `sd_pred` | +0.086 |
| `mean_pred` | −0.040 |
| `mean_ood` | +0.000 |

### Is the negative result trustworthy?

Yes — the machinery was validated against a positive control.
`tests/test_calibration.R` case 8 constructs synthetic tissues whose offset
genuinely is a linear function of mean purity; the calibrator recovers it
(LOTO R² > 0.8). Case 9 confirms it does *not* manufacture signal from noise.
So "no relationship found" reflects the data, not a broken function.

> **Methodological note.** The first version of the positive control failed at
> LOTO R² = 0.52, which briefly looked like a calibrator bug. The cause was the
> fixture: drawing every tissue's purity from one common distribution and then
> averaging over 30 samples left a between-tissue spread of only 0.089, barely
> above the injected noise. Real tissues differ systematically in mean purity.
> The fixture now gives each tissue its own centre. Worth recording because it
> is the same averaging effect that makes the real 30-tissue regression hard.

### What this means

The per-tissue offset — worth 1.17 MAE units, more than the model's entire
margin over the tissue null — **is not recoverable from label-free tissue
features.** With 30 tissues as 30 observations, there is not enough information.

This pushes the project toward the alternative in `docs/24` §6: **report
within-tissue relative position rather than an absolute HRD value.** That is
also the better match to the clinical question ("is this tumour more scarred
than typical for its type?"), but it has a hard limitation for N-of-1 pediatric
use: a relative rank needs a reference cohort of the same tumour type, which a
single pediatric case does not have.

**This is an open scientific problem, not a coding task.** Options not yet
tested are listed in §4 below.

---

## C1b — Few-shot calibration: works, but only for the tissues that need it

Tested with `scripts/fewshot_calibration.R`. No refit required: an offset
correction is a post-hoc shift of predictions that already exist, so this is
computed exactly from run 01's held-out predictions. 200 random draws per
(tissue, k), scored only on the samples **not** used to estimate the offset.

### The comparator matters more than the result

It is tempting to compare a few-shot-corrected model against the original
tissue-mean null. That comparison is rigged: it gives the model *k* labels and
the null none. If a clinician has *k* labelled samples from a new tumour type,
they can ignore the model entirely and predict `mean(k labels)` for every future
patient. **That** is the honest competitor, and it also improves with *k*.

### Results (29 tissues, equal weight)

| k | model uncorrected | model few-shot | model oracle | **null few-shot** | offset error | tissues helped | **beats fair null** |
|---|---|---|---|---|---|---|---|
| 3 | 8.43 | 8.67 | 7.44 | 9.99 | 4.35 | 8/29 | 18/29 |
| 5 | 8.43 | 8.20 | 7.44 | 9.52 | 3.31 | 11/29 | 20/29 |
| 10 | 8.43 | **7.86** | 7.44 | 9.14 | 2.36 | 14/29 | 20/29 |
| 20 | 8.44 | **7.68** | 7.45 | 8.95 | 1.59 | 17/29 | 20/29 |

Fraction of the achievable (oracle) gain captured: **k=3 → −24%**, k=5 → 23%,
**k=10 → 58%**, k=20 → 77%.

### Reading

**k=3 actively hurts.** The offset estimated from 3 samples has a mean error of
4.35 units against a true offset SD of 4.83 — the estimate is as noisy as the
quantity being estimated, so correction injects more error than it removes.

**k=10 is the practical knee.** It captures 58% of the oracle gain and cuts MAE
8.43 → 7.86. Below that the estimate is too noisy; above it returns diminish.

**The model beats the fair null in 20 of 29 tissues at every k ≥ 5**, and that
count does not improve with more labels — which is the important point. The
model's advantage over "just average your k labels" is a *fixed* property of
whether it has within-tissue signal in that tissue, not something more labels
buy.

### The gain is entirely concentrated in mis-levelled tissues

At k=10, `cor(gain, |true tissue offset|) = **0.991**`. Few-shot calibration
helps exactly the tissues that were badly mis-levelled and mildly harms the
already-calibrated ones:

| Helped most | MAE before → after | | Harmed | MAE before → after |
|---|---|---|---|---|
| THCA | 8.28 → **3.45** | | KIRC | 4.78 → 5.06 |
| SARC | 17.95 → **14.18** | | PAAD | 7.57 → 7.78 |
| PCPG | 9.08 → **5.65** | | KICH | 3.97 → 4.16 |
| UCEC | 11.77 → **9.42** | | ESCA | 9.96 → 10.12 |

Only 14 of 29 tissues improve. The mean gain is real but comes from a minority
of tissues with large offsets.

Note THCA and PCPG again: few-shot fixes most of their offset, but both still
lose to their own tissue-mean null (THCA null MAE = 0.816 vs corrected 3.45).
For near-constant tissues nothing rescues the model.

### Status

**C1b VIABLE with caveats.** Roughly 10 labelled samples from a new tumour type
recover over half the achievable calibration gain. For the pediatric
application this is a concrete, modest ask — but it is still an ask, and it
assumes ~10 pediatric HRD measurements can be obtained. It does **not** solve
the N-of-1 case, where by definition there is one patient and no cohort.

---

## C2 — Zero floor

### Clipping: free, small, safe

Clipping needs **no refit** — it is a post-hoc transform of predictions already
on disk, so it was evaluated exactly on all 7,065 development samples.

| Quantity | Value |
|---|---|
| Negative predictions | 227 (**3.2%** of samples) |
| Most negative prediction | −11.53 |
| MAE raw | 9.0494 |
| **MAE clipped at 0** | **8.9717** |
| Gain | 0.0777 MAE units |
| Skill raw | 0.0850 |
| **Skill clipped** | **0.0928** |
| Spearman(raw, clipped) | **1.000** — ranking exactly preserved |

**Adopt clipping.** It is free, removes a priori impossible values, and cannot
change any ranking metric. It is not, however, a large gain — 3.2% of samples
are affected.

### log1p: large MAE gain, mixed correlation effect

Array `323169191`, four folds, completed 2026-09-17 (46–60 min each). Compared
sample-for-sample against the matched run 01 folds.

| Tissue | n | MAE base | MAE log1p | MAE null | skill base | skill log1p | r base | r log1p |
|---|---|---|---|---|---|---|---|---|
| BRCA | 743 | 13.18 | **12.06** | 15.47 | 0.148 | **0.220** | 0.632 | 0.538 |
| UCEC | 403 | 11.77 | **8.63** | 15.22 | 0.227 | **0.433** | 0.758 | 0.745 |
| PCPG | 160 | 9.10 | **4.57** | 4.25 | −1.142 | **−0.077** | 0.110 | 0.108 |
| THCA | 464 | 8.28 | **1.99** | 0.79 | −9.552 | **−1.540** | −0.032 | **0.171** |

Sample-weighted over the four folds: **MAE 11.204 → 7.964** (−29%), mean
within-fold r 0.440 → 0.450.

**Every tissue improves on MAE, several dramatically.** The largest gains are
exactly where predicted: THCA (8.28 → 1.99) and PCPG (9.10 → 4.57), the
genomically quiet tumours the model was over-predicting by ~8 units. On the log
scale the model can finally express "this genome is calm."

### The catch, stated plainly

**Correlation falls for the two well-behaved tissues** (BRCA 0.632 → 0.538,
UCEC 0.758 → 0.745) while rising for the two failures (THCA −0.032 → 0.171).
This is the predicted trade-off appearing in the data: compressing the upper
range improves calibration for low-HRD tumours and costs discrimination among
high-HRD ones. BRCA is where HRD-high cases actually matter, and it loses 0.09
of correlation.

Note also that THCA and PCPG remain **worse than their tissue-mean null** even
after the transform (skill −1.54 and −0.08). For near-constant tissues the null
is extremely hard to beat — THCA's null MAE is 0.785 — so log1p converts a
catastrophic failure into a modest one rather than into a success.

### Recommendation

**Adopt clipping unconditionally** (free, ranking-preserving, +0.008 skill).

**REJECT log1p for the primary model** — see the full 30-fold result below.

---

## C2 final — full 30-fold log1p run: REJECTED

Array `323176856` (26 folds) + `323169191` (4 folds) = 30/30, merged
2026-09-17. Zero failures, 0/30 lambda boundary hits.

| Metric | Baseline (run 01) | log1p | Verdict |
|---|---|---|---|
| Pooled MAE | 9.049 | **8.827** | log1p better |
| Skill vs tissue null | 0.0850 | **0.1075** | log1p better |
| **Within-tissue Pearson** | **0.6124** | 0.5221 | **baseline better by 0.090** |
| Within-tissue Spearman | **0.6013** | 0.5899 | baseline better |
| Within-tissue MAE | **7.881** | 8.234 | baseline better |
| Permutation p | 0.001 | 0.001 | tie |
| Negative predictions | 227 | **0** | log1p better |

### The four-fold preview was misleading

The 4-fold test suggested log1p was close to neutral on correlation
(0.440 → 0.450 pooled). **The full run reverses that**: within-tissue Pearson
drops **0.6124 → 0.5221**, a loss of 0.090. The four folds chosen for the
preview happened to include both tissues where log1p helps most (THCA, PCPG);
they were not representative.

Across all 30 tissues, **correlation falls in 20 of 30** while MAE improves in
only 16 of 30.

| Biggest correlation losses | r base → log1p | | Biggest gains | r base → log1p |
|---|---|---|---|---|
| ESCA | 0.508 → 0.350 | | THCA | −0.032 → **0.171** |
| BLCA | 0.675 → 0.529 | | KIRP | 0.500 → **0.672** |
| LUSC | 0.605 → 0.467 | | KIRC | 0.577 → **0.686** |
| STAD | 0.738 → 0.623 | | DLBC | 0.555 → 0.654 |

### Why this settles it against log1p

The MAE and skill gains are real but small (skill +0.022). The correlation loss
is large (−0.090) and lands precisely where it hurts: **BLCA, STAD, LUSC and
ESCA are high-HRD tissues where discriminating more- from less-scarred tumours
is the clinically useful task.**

Decisively, `docs/25` §C1 established that absolute calibration on an unseen
tissue is **not achievable** without labels. That forces the project toward
**within-tissue ranking**, and ranking quality is measured by correlation — the
exact metric log1p degrades. log1p optimises the metric we are being forced to
abandon (absolute MAE) at the cost of the one we must rely on.

### Status

**C2 RESOLVED.** Adopt clipping. Reject log1p for the primary model.

Retain log1p as a **documented option for the quiet-tumour regime** — it is the
only thing tested that gives THCA any within-tissue signal at all
(−0.032 → 0.171). A per-tissue choice of transform is *not* recommended without
a principled selection rule fitted inside the LOCO loop, which does not exist.

**Note on zero inflation:** `log1p` maps 0 → 0, so the 14.3% point mass at
exactly zero remains irreproducible by any continuous regressor. `log1p`
addresses right skew, **not** zero inflation. Genuine zero inflation needs a
hurdle / two-part model (classify zero vs non-zero, then regress the positives),
which is recorded as future work.

---

## 4. Revised priorities

| Was | Now |
|---|---|
| C3 purity — most serious | **Downgraded** — artefact of C1, not purity capture |
| C1 calibration — tractable fix | **Escalated** — covariate approach empirically fails |
| C2 zero floor — cheap win | **Confirmed cheap, confirmed small** (clip: +0.008 skill) |

Untested options for C1, in rough order of promise:

1. **Report within-tissue ranks only.** Honest; matches the clinical question;
   needs a same-type reference cohort at prediction time.
2. **Few-shot calibration.** A handful of labelled samples from the new tissue.
   Strongest correction, but assumes pediatric HRD labels exist — the very thing
   that is scarce.
3. **Richer tissue covariates.** The ones tested are summaries of the model's
   own output. Methylation-derived tissue descriptors (mean beta over selected
   probes, cell-composition estimates) were not tested and are not obviously
   more informative, but are cheap to add.
4. **Accept the offset and report it.** Publish the model with a stated
   per-tissue calibration error of ±3.5 units and let readers judge.

The CNS lock stays closed.

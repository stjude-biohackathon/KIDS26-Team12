# KIDS26 Team 12 — methylation prediction of HRD genomic-scar burden

St. Jude KIDS26 Biohackathon. Repo:
`/research/groups/shelagrp/home/esavage/KIDS26/KIDS26-Team12`
(GitHub `stjude-biohackathon/KIDS26-Team12`, branch `main`).

## The question

Can a frozen DNA-methylation model predict independently measured HRDsum in a
cancer type it has never seen — ultimately pediatric high-grade glioma?

**Framing that matters:** canonical HRDsum = HRD-LOH + LST + TAI and is measured
from *allele-aware* genomic data. Methylation arrays cannot observe LOH or
telomeric allelic imbalance. Methylation is therefore a **predictor of an
independently measured continuous reference HRDsum**, never a re-derivation of it.
Nothing here is a validated clinical HRD assay or a treatment-selection result.

## STATE AS OF 2026-09-22: the experiment is finished

**The CNS lock was opened on 2026-09-18 and is CONSUMED. It cannot be reused.**
Ledger: `results/LOCKED_EVALUATION_LEDGER.tsv` (two rows — one session, one per
candidate, permitted by `docs/28` §5). LSF job `323264626`.

The answer to the project's question, in one line: **within-tissue rank transfers
weakly to an unseen lineage; absolute cross-lineage calibration does not.**
This is outcome A of the three pre-declared in `docs/28` §7. The deployment gate
stays shut.

Do not re-score CNS. Do not tune against it. It is evaluation evidence now.
**PBTP remains untouched and is the only genuinely independent test left** —
protect it.

## Data

| Path | What |
|---|---|
| `data/processed/beta.tsv` | 28 GB, **336,480 probes × 7,707 samples**. Probes in ROWS, sample barcodes in COLUMNS, first col `probe_id`. sha256 `e3642d30…`. There is **no separate CNS matrix** — all 7,707 samples live in this one file |
| `data/processed/master_samples.tsv` | 7,707 rows, one per specimen. Join key `sample_id`. Has `cancer_type`, `HRDsum`, `HRD_LOH`, `LST`, `TAI`, `purity`, `ploidy`, `patient_id`, `partition`. **CRLF line endings** — `fread` handles it, naive `awk -F'\t'` field matching does not |
| `config/shared_autosomal_probes.txt` | 384,640-probe HM450/EPIC-v1 bridge (allowlist; 336,480 is its intersection with the matrix) |

Cohort: 7,065 development across **30 non-CNS cancer types**, plus **642 locked
CNS = 135 GBM + 507 LGG**. HRDsum median 14, IQR 4–28; only 10.6% ≥ 42; 14.3%
are exactly 0.

**Known caveats.** OV is n=10 (TCGA ovarian methylation is mostly 27k, excluded by
the 450k bridge) — the canonical HRD cancer is effectively missing. HRD labels join
at 15-character sample-type level, not vial. Per-specimen array QC has never gated
this cohort.

### Two loading gotchas

1. **Never `read.csv`/`fread` the matrix casually** — 28 GB on disk, ~173 GB peak
   RSS through the transpose chain.
2. **Never pass probe IDs to `fread(select=)`.** `select=` filters COLUMNS; probes
   are ROWS. It warns rather than errors and the transpose yields **0 × n_probes**.
   Shipped once (`4ec3493`), reverted. `tests/test_load_path.R` guards it.

## THE SHIPPED MODEL: V3-abs

`results/frozen_v3_2026-09-18/frozen_nonCNS.rds`, sha256 `df9f7e82…`,
alpha 0.1, lambda 1.3598, **905 non-zero coefficients**, 5,664 training +
1,401 calibration, conformal q = **21.10** at 95%.

**V3-abs = run-01's absolute model with exactly ONE change:** the 5,000-probe
unsupervised filter ranks probes by pooled **within-tissue** variance
(training-fold rows only) instead of pooled **total** variance.

**This is NOT `docs/26`'s "V3"**, which was V2's already-refuted relative target
*plus* that filter — two changes stacked, one of them known broken. The
distinction was pre-registered in `docs/27` §3 before either result existed, and
it matters: a collaborator reported a null "V3" result that was most likely
docs/26-V3. **That discrepancy was never resolved** and is recorded as an open
question in `docs/28` §2.

Candidate A (run 01, `feature_rank="pooled"`) is frozen at
`results/frozen_2026-09-18/frozen_nonCNS.rds`, sha256 `0e975cd6…`, 774 non-zero,
conformal q = 21.73.

**Both conformal intervals are ±21 HRD units — wider than the label's own IQR
(4–28).** Single-sample absolute precision is poor and saying so is mandatory.

## Results

### Source domain, full 30-fold LOCO (`results/tables/tableA2_A_vs_V3_full30.md`)

| Metric | A (run 01) | V3-abs | Δ |
|---|---|---|---|
| Macro within-tissue Pearson | 0.5197 | 0.5036 | −0.0161 |
| Macro within-tissue Spearman | 0.4747 | 0.4472 | −0.0275 |
| Pooled MAE | 8.9717 | 8.4456 | −0.526 |
| Mean \|tissue bias\| | 3.4628 | 2.9737 | −0.489 |
| **Tissue R² of predictions** | **0.5584** | **0.4672** | **−0.0913** |
| Tissue R² of truth | 0.3414 | 0.3414 | cancels |
| Macro skill vs oracle null | 0.419 | 0.630 | +0.211 |
| Tissues with positive r | 29/30 | 29/30 | — |

**Headline:** ΔR²_excess = **0.0913**, 95% CI [0.0824, 0.0996], p < 0.0001 —
**42.0%** (CI 37.9–46.5%) of the excess lineage structure removed. Because both
models predict the same labels on the same cohort, R²_truth cancels exactly and
this quantity *is* R²_pred,A − R²_pred,B.

**The honest caveat:** the paired per-tissue ranking change is **not**
significant — Pearson mean −0.0161 (CI [−0.0410, 0.0051], Wilcoxon p = 0.428),
Spearman −0.0275 (p = 0.477); only 12/30 and 14/30 tissues improved. V3-abs is a
**targeted reduction in lineage imprinting at no measurable ranking cost**, not a
general accuracy improvement. Do not oversell it.

### ROC, development LOCO (`results/figures/fig09_*`, `fig10_*`)

| Framing | A | V3-abs |
|---|---|---|
| Pooled, exploratory HRDsum ≥ 42 | 0.862 | 0.866 |
| **Tissue identity ONLY (control)** | **0.743** | **0.723** |
| Within-tissue top quartile | 0.779 | 0.777 |

The control is the point: replacing every sample with its own tissue's *mean*
prediction — zero within-patient information — still reaches 0.74. **Most of the
pooled AUC is lineage recognition.** The within-tissue number is the defensible
claim.

### CNS external test — the headline result

PRIMARY = V3-abs. Analysis in `results/cns_eval_2026-09-18/analysis_v3/`.

| | n | Pearson | MAE | bias | slope |
|---|---|---|---|---|---|
| GBM | 135 | 0.320 [0.173, 0.461] | 8.86 | **+7.99** | 0.358 |
| LGG | 507 | 0.329 [0.234, 0.421] | 5.13 | +1.62 | 0.578 |
| pooled | 642 | 0.261 [0.183, 0.335] | 5.92 | +2.96 | 0.358 |

SENSITIVITY = Candidate A (reported per `docs/28` §5 clause 3, which forbids
reporting only the better one): GBM r = 0.293, LGG r = 0.371. **CIs overlap
heavily; neither candidate is clearly better on CNS.**

**What failed, decisively:** the **oracle tissue-mean null beats the model
everywhere** — skill −1.009 (GBM), −0.048 (LGG), −0.233 (pooled). "This is GBM,
use the GBM average" outperforms the model's absolute predictions.

**Placement among the 30 source tissues** (the primary interpretive frame — "would
CNS look unusual as the 31st held-out tissue?"): GBM Pearson 13.3rd percentile
(28th of 32), LGG 16.7th (26th of 32), GBM bias 29th of 32. **CNS is a poor
relation, not a typical held-out tissue.**

**The most important scientific finding of the whole project:** V3-abs's 42%
lineage reduction, measured on source tissues, **did not translate into better
CNS calibration**. Reducing lineage imprinting *within* the training distribution
was not sufficient for transfer *outside* it. That is the result worth publishing
and the one that should shape the next model.

### Two CNS-specific defects

- **OOD flagged only 8 of 642** (all LGG). The frozen reportability rule did
  **not** warn that CNS was out of distribution. It should have. (Blocker B17.)
- **HRDsum ≥ 42 is vacuous in CNS**: GBM has **0** samples above it, LGG exactly
  **1**. No binary claim is possible in this lineage; the AUCs at that cutoff are
  degenerate. (Blocker B18.)

## What was tried against C1 and failed

Five pre-registered attacks; four failed. Record them — they are the scientific
content, and the project has an explicit rule against hiding them.

1. **Label-free tissue-offset regression** — best LOTO R² = **−0.116**, worse than
   predicting the mean.
2. **Zero-shot percentile remap** — adds *exactly zero* ranking information (a
   monotone map cannot reorder) and *raised* tissue R² 0.562 → 0.593.
3. **Relative-target V1/V2** — failed 3 of 7 pre-registered gates (macro ρ −0.219,
   BRCA −0.405).
4. **V4 supervised lineage-penalized selection** — failed the **mechanism gate
   G3**: tissue R² went **UP**, 0.5219 → 0.5361. The experiment designed to
   penalize lineage produced more of it. `lambda_pen` was inner-fold selected and
   came out **0 for BRCA** — on the largest fold the tuner preferred no penalty.
   Branch terminated per `docs/27` §6.
5. **V3-abs** — passed 7/7 sentinel gates and 7/7 selection criteria. Shipped.

**What works but needs labels:** EB shrinkage few-shot calibration (τ²=22.3,
σ²=108.6, w ≈ k/(k+4.87)); k=3 recovers +40% of oracle gain, k=10 recovers 58%.
This is the *labelled* path and is not n-of-1.

## Code layout

**Main pipeline.**
- `R/model.R` — `fit_preprocess` (now takes `feature_rank` ∈ {`pooled`,
  `within_tissue`, `lineage_penalized`} + `cancer` + `y`), `within_tissue_var`,
  `lineage_meta_stats`, `fit_en`, `predict_en`, `metrics`, `null_panel`,
  `within_tissue_metrics`, `tissue_identity_r2`, `assert_partition_matches_cns`.
- `R/rank_model.R` — relative-target (V1/V2) model, refuted but retained.
- `R/calibration.R`, `R/provenance.R`, `R/figures.R` (shared `theme_kids26()`).
- `scripts/loco_one_fold.R` — `<beta> <master> <out> <fold_index> [transform] [feature_rank]`.
  `fold_index` indexes `sort(unique(cancer_type[!cns]))` (3=BRCA, 10=KICH, 26=THCA, 28=UCEC).
- `scripts/freeze_final_model.R` — `<beta> <master> <out> [feature_rank]`. Phase-2 freeze only.
- `scripts/predict_frozen.R` — `<model> <beta> <out> [--allow-fixture] [--metadata=…] [--locked-evaluation]`.
- Scorers: `score_v3_sentinel.R`, `score_v4_sentinel.R`, `score_candidate_selection.R`,
  `compare_A_vs_V3.R`, `analyze_cns.R`, `component_decomposition.R`.
- Figures: `make_figures.R`, `make_roc_figure.R` (CLI-parameterised),
  `make_cns_figures.R`, `make_story_slide.R`.

**Independent replicate (Kayode).** `pipeline_glmnet/` — a *second implementation*
from `config/analysis_protocol.json`. **Do not refactor the two into one**; their
independence is the scientific point.

### Model design (do not silently change)

Nested LOCO. Inner folds are **one per training cancer type**, so hyperparameters
are selected for *cross-tissue* generalisation. Preprocessing relearned inside
every training fold; feature selection unsupervised by default; `standardize=FALSE`;
lambda path from the fold's own `lambda.max`. Seed 260910.

## Pre-unlock hardening (commit `adf420b`) — read before trusting old results

Four defects on the shipping path, all fixed 2026-09-18, plus one that was worse:

- **`predict_frozen.R` never loaded glmnet**, so `predict()` had no S3 method —
  **the shipping inference path could not score anything at all.** Pre-existing.
- It never applied the adopted C2 clipping, so it emitted negative HRDsum. Now
  clips via the bundle's own `target_transform`, retaining
  `predicted_reference_HRDsum_raw`.
- No locked-partition guard and no scoring record → `--locked-evaluation` opt-in
  plus the append-only ledger.
- Two sources of truth for the lock → `assert_partition_matches_cns()` requires
  two-way agreement between the hardcoded `c("GBM","LGG")` and `partition`.
- `loco_merge.R` completeness gate only globbed `predictions_*.tsv` → now checks
  metrics/nullpanel counts and that folds share one `target_transform` AND
  `feature_rank`; ships tissue-identity R² every run.

## Conventions

- **Tests are plain `stopifnot()` scripts, not testthat.** All 8 pass:
  `Rscript tests/{smoke_model,test_calibration,test_provenance_gate,test_c1_rank_model,test_load_path,test_app_governance,test_v3_feature_rank,test_v4_lineage_penalized}.R`
- Cluster is **LSF, not Slurm**. **`rusage[mem]` is PER SLOT** —
  `-n 4 -R "rusage[mem=60GB]" -M 240GB`. Getting this wrong leaves jobs PENDING forever.
- Use queue **`priority`** (rhel8_rome, 46×1 TB + 14×1.9 TB), not `biohackathon` (~2 nodes).
- A fold costs ~45–85 min and ~173 GB peak.
- `results/` is **gitignored** — presentation artifacts must be `git add -f`'d.
- Never commit without checking `git status`; teammates push concurrently.
- **Mermaid validation gotcha:** `quarto render` to HTML **passes on invalid
  mermaid** (it embeds source for client-side rendering). Real validation needs
  `mermaid-format: png`, which renders server-side.

## Working style that produced these results

**Pre-register gates before the run, in a committed file, while the output
directory is verifiably empty.** Every gate here was frozen in git before the
results it judged existed: `docs/27` §5 (V3/V4 G1–G7), `docs/28` §4 (S1–S7),
`docs/29` (pre-unlock designation). Score them mechanically with a script that
hard-codes the thresholds, so scoring is arithmetic rather than judgement applied
after seeing numbers.

**Sentinels can stop a run; they cannot authorise a conclusion.** B16 exists
because a 4-fold log1p preview inverted the full 30-fold verdict. It earned its
keep twice: V3-abs's sentinel tissue R² (0.2991) was far better than its full-run
value (0.4672).

**THCA is a negative control, never a win condition** (near-constant within-tissue
target, r ≈ −0.03).

**Pair every REJECT test with an ACCEPT test** — B5 was a QC check comparing a
value against the constant that produced it.

**Review teammate commits, don't just merge them.** Three shipped broken: B4's
`select=`, B10's allowlist, and the glmnet omission above.

## Blockers

**Open:** B10 (app governance — raw identifiers, per-session `results`), B11 (PBTP
EPIC generation unconfirmed), B12 (IDAT bridge, deferred), B15 (`renv.lock`
snapshot still owed), B16 (sentinel representativeness — process), **B17 (OOD rule
failed to flag CNS)**, **B18 (≥42 threshold vacuous in CNS)**, B20 (CNS bsub
script-contract mismatch: `predict_frozen.R` scores the whole matrix while
`analyze_cns.R` wants CNS-only input).

**Closed:** B1–B9, B13, B14, B19 (the `adf420b` hardening set). C2 resolved
(clipping adopted, log1p rejected). C3 downgraded (purity largely exonerated;
partial correlation *rose* 0.612 → 0.621). C1 **characterised and tested** —
mitigated ~42% in-distribution, still decisive out of it. C4 (OV n=10) open
disclosure.

## Key docs

`docs/30_FINAL_ANALYSIS_2026-09-18.md` (the full write-up — read this first),
`docs/27_V3_PREREGISTRATION.md` (gates + the deployment/external-evaluation gate
split), `docs/28_CANDIDATE_SELECTION_AND_CNS_PROTOCOL.md` (S1–S7 + CNS protocol),
`docs/29_PRE_UNLOCK_RECORD.md` (audit record), `docs/21` (blocker ledger),
`docs/23` (schematic), `docs/24` (run 01, superseded as shipped model),
`docs/25`/`docs/26` (C1 investigation), `pipeline_glmnet/RESULTS.md`.

## Next steps

Everything below is **post-hackathon** work. The hackathon experiment is done.

1. **Do not re-score CNS.** It is spent. Any new model is validated on PBTP or
   nothing.
2. **B17 first** — an OOD rule that fails to flag an entire unseen lineage is a
   deployment-blocking defect, and fixing it does not require new labels.
3. **The real research question the CNS result poses:** why did in-distribution
   lineage reduction not produce out-of-distribution calibration transfer?
   Candidate approaches — hierarchical/multitask elastic net, anchor regression,
   group DRO, tissue-subspace removal, explicit few-shot domain calibration.
4. **The labelled few-shot path is the only one with demonstrated traction.** If
   pediatric deployment is the goal, budget for k≈10 labelled PBTP cases rather
   than hoping for zero-shot transfer.
5. Protect PBTP. It is the last untouched test set.

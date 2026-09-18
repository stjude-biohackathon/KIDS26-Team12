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

## Current state (2026-09-18)

A real result exists and is reframed. The model works as a **within-tissue ranker**,
not an absolute predictor on unseen tissues. Two blockers gate everything:
**C1 (calibration offset) is unresolved and has resisted three independent
attacks.** The CNS lock is still closed.

Ask before: opening the CNS lock, running a full 30-fold array, or touching PBTP.

## Data

| Path | What |
|---|---|
| `data/processed/beta.tsv` | 28 GB, **336,480 probes × 7,707 samples**. Probes in ROWS, sample barcodes in COLUMNS, first col `probe_id`. sha256 `e3642d30…` |
| `data/processed/master_samples.tsv` | 7,707 rows, one per specimen. Join key `sample_id`. Has `cancer_type`, `HRDsum`, `HRD_LOH`, `LST`, `TAI`, `purity`, `ploidy`, `patient_id`, `partition` |
| `data/processed/beta.provenance.json` | `engineering_only:false`, `n_probes:336480`, `probe_allowlist_sha256:f359e43b…`. The provenance gate reads this |
| `config/shared_autosomal_probes.txt` | 384,640-probe HM450/EPIC-v1 bridge (the *allowlist*; 336,480 is its intersection with the matrix) |

Cohort: 7,065 development across **30 non-CNS cancer types**, plus **642 locked
GBM+LGG**. HRDsum median 14, IQR 4–28; only 10.6% ≥ 42 (the exploratory threshold
in `config/analysis_protocol.json`). 14.3% are exactly 0.

**Known data caveats.** OV is n=10 (TCGA ovarian methylation is mostly 27k, excluded
by the 450k bridge) — the canonical HRD cancer is effectively missing. HRD labels
join at 15-character sample-type level, not vial, so every join is a coarser match
than aliquot identity. Per-specimen array QC (detection p-values, sex checks) has
never gated this cohort.

### Two loading gotchas

1. **Never `read.csv`/`fread` the matrix casually** — 28 GB on disk, ~173 GB peak
   RSS through the transpose chain.
2. **Never pass probe IDs to `fread(select=)`.** `select=` filters COLUMNS; probes
   are ROWS. It warns rather than errors, returns only `probe_id`, and the
   transpose yields **0 × n_probes — zero samples with the probe count intact**.
   This shipped once (commit `4ec3493`, B4) and was reverted.
   `tests/test_load_path.R` guards it.

## Code layout

**Main pipeline (Evan).**
- `R/model.R` — `fit_preprocess`, `apply_preprocess`, `inner_folds`, `lambda_path`,
  `fit_en`, `predict_en`, `metrics`, null/purity controls.
- `R/rank_model.R` — relative-target (rank) model: `relative_target`, `fit_rank_en`,
  `predict_rank`, `score_to_percentile`, pairwise ranker prototype.
- `R/calibration.R`, `R/provenance.R` (shared provenance gate).
- `scripts/loco_one_fold.R` — `<beta> <master> <out_dir> <fold_index> [transform]`.
  `fold_index` indexes `sort(unique(cancer_type[!cns]))` (1=ACC … 3=BRCA, 10=KICH,
  26=THCA, 28=UCEC … 30=UVM).
- `scripts/c1_rank_one_fold.R` — same CLI plus `[target_type] [weight_mode] [feature_rank]`.
- `scripts/c1_shrinkage_calibration.R`, `scripts/c1_zeroshot_percentile.R`.
- `scripts/loco_merge.R`, `scripts/predict_frozen.R`, `scripts/train_baseline.R`, `app/app.R`.

**Independent replicate (Kayode).** `pipeline_glmnet/` — a *second implementation*
written from `config/analysis_protocol.json`, with its own `R/model.R`,
`01_prepare_features.R` … `04_domain_shift.R`, and `results_summary/`. Independently
reproduces C1 and adds a tissue-identity permutation control. **Do not refactor the
two into one**; their independence is the scientific point.

### Model design (do not silently change)

Nested leave-one-cancer-out. Inner folds are **one per training cancer type**, so
hyperparameters are selected for *cross-tissue* generalisation. Preprocessing is
relearned inside every training fold; feature selection is unsupervised (top 5,000
by **pooled** variance, tie-broken by probe name); `standardize=FALSE` because data
is pre-standardised; lambda path derived from the fold's own `lambda.max` (B8).
Hyperparameter loss is macro-averaged across tissues.

## Results

### Run 01 — the headline (array `323078995`, `results/loco_run01/`)

| Metric | Value |
|---|---|
| Within-tissue Pearson / Spearman | **0.612 / 0.601** |
| Pooled MAE | 9.049 (7.881 with oracle offset removed) |
| Skill vs tissue-mean null | +0.085 |
| Permutation p | 0.001; positive in 29/30 tissues |
| Tissue identity explains | **56.2%** of the prediction |

Kayode's independent replicate: macro MAE 8.68 vs null 13.16, pooled Spearman 0.738,
median within-cancer Spearman 0.476, pooled AUC(HRDsum≥42) **0.841**. His
tissue-identity permutation control collapses these to MAE 12.62, Spearman 0.374,
AUC 0.601 — i.e. **tissue identity alone explains much but not all of the signal.**

Per-tissue within-tissue r ranges from KICH 0.803 / UCEC 0.758 / BRCA 0.632 down to
**THCA −0.032, PCPG 0.110, CHOL 0.161** — several tissues carry no usable signal.

### C1 — per-tissue calibration offset. **UNRESOLVED, escalated.**

Mean |offset| 3.46 HRD units (SD 4.83, range SARC −15.44 to PCPG +8.29). Removing
it costs more than the model's entire margin over the null.

Three independent attacks have failed:
1. **Label-free offset regression** — best LOTO R² = **−0.116**; offset uncorrelated
   with every available covariate (|r| ≤ 0.267).
2. **Zero-shot percentile remap** — adds *exactly zero* ranking information (a
   monotone map cannot reorder) and *raises* tissue R² 0.562 → 0.593.
3. **Trained relative-target model (V1/V2)** — sentinel array `323195423`
   (BRCA/KICH/THCA/UCEC) **failed 3 of 7 pre-registered gates**: macro ρ −0.219,
   BRCA −0.405, and tissue R² **0.576 vs 0.527** for the absolute model (floor:
   0.001). Not manufacturing signal — all folds beat permutation at p ≤ 0.004 — it
   just learned less while staying more confounded.

**Diagnosed mechanism:** the 5,000-probe filter ranks by **pooled** variance *before*
the target is consulted, so it selects lineage-discriminating probes. Changing the
target cannot undo that.

**Still untested and the highest-value next experiment: V3**
(`feature_rank="within_tissue"`), which ranks probes by training-only within-tissue
variance. Implementation written and benchmarked (~2.2 min for 336,480 probes).

**What works: labelled calibration.** EB shrinkage toward a cross-tissue prior
(τ²=22.3, σ²=108.6, so w ≈ k/(k+4.87)) **inverts the k=3 verdict** — the unshrunk
mean is harmful (−22% of oracle gain), shrunk it recovers +40%. k=10 recovers 58%.
Beats a *fair hierarchical null* in 25–26 of 29 tissues. Caveats: only 14/29 tissues
helped at k=3, benefit confined to tissues with |offset| ≥ 5 (0/11 for |offset| < 2),
and uniform shrinkage over-shrinks genuinely large offsets.

### Pediatric consequence (docs/26 §7b)

Both n-of-1 routes are closed. **95% predictive interval for an unseen tissue's
offset: −9.7 to +8.1 HRD units** — 74% as wide as the label's own IQR. Among the
1,171 patients within ±10 units of threshold 42, a tissue-level offset flips the
HRD-high call for a median of **9.1%** (up to 79.2%). The rank escape hatch is
majority tissue identity, which hurts pediatrics *more* because 89% of the cohort
sits below threshold and pHGG is expected in that low-HRD regime where ordering —
not a binary call — is the useful output. Only the **labelled** path survives, and
it needs PBTP labels we do not have.

### C2 — resolved. C3 — downgraded. C4 — open disclosure.

- **C2:** clipping at 0 **ADOPTED** (MAE 9.049→8.972, Spearman with raw exactly
  1.000). **log1p REJECTED** on the full 30-fold run: better absolute metrics
  (MAE 8.827) but within-tissue r falls **0.612→0.522**, worse in 20/30 tissues,
  concentrated in high-HRD tissues (BLCA, STAD, LUSC, ESCA). Since C1 forces the
  project toward ranking, log1p optimises the metric we must abandon. A 4-fold
  preview had pointed the *opposite* way — it happened to contain both tissues where
  log1p helps most (see B16).
- **C3 purity:** largely exonerated. Partial correlation controlling for purity
  *rose* 0.612 → 0.621. The apparent inversion was an artefact of C1. Residual
  gradient survives, so downgraded not closed. Fine as a spoken caveat.
- **C4:** OV n=10 — state in every presentation.

## Conventions

- **Tests are plain `stopifnot()` scripts, not testthat.** Run explicitly:
  `Rscript tests/{smoke_model,test_calibration,test_provenance_gate,test_c1_rank_model,test_load_path,test_app_governance}.R`
  (Python: `python -m unittest discover -s tests`). All currently pass.
- Cluster is **LSF, not Slurm** (`bsub`/`bjobs`/`bqueues`). **`rusage[mem]` is PER
  SLOT** — `-n 4 -R "rusage[mem=60GB]" -M 240GB`. Getting this wrong leaves jobs
  PENDING forever, which looks like a busy queue.
- Queues: `biohackathon` has only ~2 usable nodes. **`priority` maps to `rhel8_rome`
  (46 nodes @ 1 TB, 14 @ 1.9 TB)** — use it; folds start immediately, though shared
  load stretches a 45–55 min fold to 65–85 min.
- A fold costs ~45–55 min and ~173 GB peak on a quiet node.
- Results land in `results/` which is **gitignored** — never assume tables are on
  GitHub.
- Never commit without checking `git status`; teammates push concurrently.
- Seeds: 260910 in the model; pick a random 1–10000 integer for new work.

## Working style for this project

Pre-register gates before sentinel runs and score them honestly afterwards
(`docs/26` §7c keeps the plan verbatim so scorecards can't be read as post-hoc).
**Sentinel previews are not representative by construction** — B16 exists because
a 4-fold preview inverted the full 30-fold log1p verdict. A sentinel may *stop* a
full run; it may not authorise one unless every gate clears.

**THCA is a negative control, never a win condition.** Its within-tissue target is
near-constant (r = −0.032), so rank metrics there are degenerate; it "improved" in
the C1c sentinel purely because of that.

Prefer a gate that can *fail*: B5 was a QC check comparing a value against the
constant that produced it. Pair every REJECT test with an ACCEPT test.

## Blockers (`docs/21_BLOCKER_RESOLUTION_PLAN.md`)

B1–B9 resolved. **Open:** B10 (app governance — allowlist/provenance done; raw
identifiers and per-session `results` still open), B11 (PBTP EPIC generation
unconfirmed), B12 (IDAT bridge, deferred), B15 (`renv.lock` snapshot still owed —
`matrixStats`/`digest` added to `setup.R` but not snapshotted, deliberately, since
this node isn't verified to match the cluster env), B16 (sentinel representativeness).
B13/B14 closed with `tests/test_load_path.R` (7/7) and `tests/test_app_governance.R` (10/10).

**Review teammate commits, don't just merge them.** Two shipped broken: B4's
`select=` (above) and B10's allowlist, which omitted the one path README documents,
required a sidecar `predict_frozen.R` never writes, and had unwired `renderUI`
blocks. Both were fixed in place with the reasoning recorded at the call site.

## Key docs

`docs/24_RESULTS_LOCO_RUN01.md` (headline), `docs/25_C1_C2_C3_INVESTIGATION.md`,
`docs/26_C1_N_OF_1_WORKAROUND.md` (C1 workaround + sentinel scorecard + pediatric
consequence), `docs/21_BLOCKER_RESOLUTION_PLAN.md` (ledger + mermaid),
`docs/23_PIPELINE_SCHEMATIC.md`, `pipeline_glmnet/RESULTS.md` (Kayode's replicate).
Mermaid diagrams are validated by rendering through Quarto after edits.

## Next steps as of 2026-09-18

Agreed deliverables for the flash talk, all lock-independent:
1. Training predicted-vs-genomic HRDsum scatter.
2. **ROC both ways side by side** — absolute at HRDsum ≥ 42 *and* within-tissue top
   quartile. The contrast makes the offset problem visible in one plot. Frame as
   "ranking performance at an exploratory cutoff", not diagnostic accuracy.
3. Probe genomic-location plot with **top 20 loci annotated** (needs
   `data/raw/annotation/HM450.hg19.manifest.202209.tsv.gz` for coordinates).

Decision taken: **keep the CNS lock closed** and present run 01's 30-fold LOCO as
held-out-by-tissue performance — while labelling it a *development* result, since
clipping, log1p and the C1 investigation were all developed against it. GBM/LGG is
the only genuinely untouched data left.

Optional: run the **V3 sentinel** (4 folds, same gates) — the last open C1 question.

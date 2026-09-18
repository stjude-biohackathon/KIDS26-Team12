# Methyl-HRD: independent LOCO replicate, tissue-identity control, and C1 quantification

Kayode Raheem (@kayoderaheem) — methylation preprocessing / model harmonization
Pipeline: `pipeline_glmnet/` (this folder). Runs completed 2026-09-17/18 on the team VM
(32 cores, 125 GB). R 4.6.1, glmnet 4.1-8, data.table. Seed 260910 throughout.

This is a **second, independent implementation** of the development model, written
from `config/analysis_protocol.json` after `R/model.R` was found to be missing from
the data share. It reproduces the C1 finding with different code, and adds a
tissue-identity permutation control that was not previously run.

---

## 1. Cohort QC (assigned task)

| Step | n |
|---|---:|
| Candidate matrix columns (`matching_audit.tsv`) | 9,664 |
| Removed: not primary tumour | 1,363 |
| Removed: missing/ambiguous HRD label | 319 |
| Removed: quality excluded | 235 |
| Removed: ambiguous multiple specimens per patient | 40 |
| **Clean cohort** | **7,707** |
| of which development | 7,065 (30 cancer types) |
| of which locked_CNS (GBM/LGG, untouched) | 642 |

`data/processed/beta.tsv` already contained exactly the 7,707 kept samples; the audit
exclusions were verified, not re-applied. Tables: `samples_kept.tsv`,
`samples_removed_by_audit.tsv` in `back_end/beta_subset/`.
All kept samples carry `quality_annotation = published_450K_no_exclusion`.
Note carried forward: every HRD label is matched at 15-character sample-type level,
not vial (`label_resolves_sample_type_not_vial`).

## 2. Method

- Probes: 336,480 scanned; label-free prefilter (≤5% missing, top 50,000 by variance);
  the protocol's 5,000-feature selection then runs **inside each training fold**.
- Model: elastic net (`cv.glmnet`), alpha ∈ {0.1, 0.25, 0.5, 0.75, 1}, full lambda path,
  `lambda.1se`. Inner folds = one per training cancer (nested leave-one-cancer-out), so
  hyperparameters are tuned for cross-tissue generalisation and are deterministic.
- Outer: leave-one-cancer-out over all 30 non-CNS development cancers. GBM/LGG are never
  read by any script in this folder.
- Imputation (training medians), feature selection and scaling are all learned on
  training rows only. Predictions clipped at 0 (consistent with C2).
- Runtime: feature prep 3.8 min; 30 folds + frozen model 2.2 h on 16 cores.

## 3. Development results (30 cancers, leave-one-cancer-out)

| Metric | Real labels | Tissue-identity control | Baseline |
|---|---|---|---|
| Median within-cancer Spearman | **0.476** | 0.134 | 0 |
| Pooled Spearman | 0.738 | 0.374 | 0 |
| Pooled AUC (HRDsum ≥ 42) | 0.841 | 0.601 | 0.5 |
| Macro MAE | 8.68 | 12.62 | 8.71 (each cancer's own mean) |
| Macro MAE after 20-label intercept | **7.67** | 8.86 | 8.71 |

The control trains on HRDsum shuffled **within** each cancer type, so cancer-level means
are preserved and only lineage is learnable. It is the direct test of the tissue confound.

**Per-cancer Spearman, real vs control** (full table: `loco_metrics.tsv`,
`domain_shift_per_cancer.tsv`):

| Signal clearly beyond lineage | Real | Control | | Margin thin — caveat | Real | Control |
|---|---|---|---|---|---|---|
| STAD | 0.75 | 0.37 | | PAAD | 0.57 | 0.51 |
| LUAD | 0.67 | 0.25 | | MESO | 0.62 | 0.45 |
| UCEC | 0.67 | 0.03 | | THYM | 0.50 | 0.43 |
| BRCA | 0.65 | 0.31 | | SARC | 0.56 | 0.35 |
| HNSC | 0.61 | 0.20 | | OV (n=10) | 0.53 | 0.62 |
| PRAD | 0.60 | 0.13 | | | | |
| ACC | 0.59 | 0.00 | | | | |
| BLCA | 0.58 | −0.12 | | | | |

Weak or absent in genomically quiet types: THCA −0.05, PCPG 0.06, TGCT 0.12, KIRC 0.23,
CHOL −0.07 (n=35). AUC(≥42) is strongest in PRAD 0.92, STAD 0.86, UCEC 0.85, LUAD 0.84,
BRCA 0.82.

## 4. C1 (per-tissue offset), quantified

Independent reproduction of the escalated C1 finding:

- Model vs tissue-mean baseline on the absolute scale: **8.68 vs 8.71** — effectively tied.
- Oracle offset (best possible per-tissue shift): **7.59**, so the whole offset problem is
  worth ~13% of MAE.
- Few-shot intercept calibration (labelled samples from the new tissue, intercept only,
  200 random draws per k, `calibration_curve.tsv`):

| k labelled samples | 0 | 3 | 5 | 10 | 20 |
|---|---|---|---|---|---|
| Mean MAE | 8.51 | 8.67 | 8.16 | 7.85 | 7.67 |

  Consistent with the team's "≈58% of achievable gain at k≈10" (this run, plain offset with
  no shrinkage, recovers ~70%).
- Biases are systematic shrinkage toward the centre, not noise: SARC −15.4, OV −11.2,
  SKCM +10.2, UCEC +8.3, THCA +8.0, ACC +7.4.

**Why label-free calibration cannot work, shown directly:** given 20 labels, the
lineage-only control lands exactly on the tissue-mean baseline (8.86 vs 8.71) — that is all
lineage ever knew. The real model given the same 20 labels goes below it (7.67). The offset
is a property of the tissue, not something present in the methylation signal.

## 5. What this supports

1. **Report within-tissue metrics as primary.** They are the part not confounded by lineage
   (0.476 vs 0.134 against the control).
2. **Always show pooled metrics next to the control.** With 30 training tissues, a
   lineage-only model reaches Spearman 0.374 / AUC 0.601 on held-out tissues, so pooled
   numbers overstate the HRD-specific signal.
3. **C1c (rank-only) is the defensible framing**; C1b (k≈10–20 labels) is the route to
   absolute values where labels exist. Zero-label absolute HRDsum on a new tissue is not
   achievable from methylation alone, and this run shows why.
4. **Per-cancer claims should be limited** to the left-hand column of §3.

## 6. Reproducing

```sh
SMOKE=1 bash run_pipeline.sh                                   # ~15 min sanity run
bash run_pipeline.sh                                           # full development LOCO
CONTROL=permute_within_cancer bash run_pipeline.sh             # tissue-identity control
Rscript 04_domain_shift.R --run_dir=results/full_none           # C1 tables
Rscript 02_train_loco.R --features=features_full \
  --out=results/full_within_z --target=within_z --cores=10      # C1c relative target
```
On LSF use `submit_hpc.sh` (memory is requested per slot, per `docs/23`).

## 7. Limitations

- The 50,000-probe prefilter uses variance across all development samples including the
  held-out cancer (label-free, but not fold-pure).
- Single permutation seed for the control; repeats would give a range.
- OV n=10 (C4) — no OV-specific claim is interpretable here.
- Not a clinical assay; predictions are of an independently measured reference HRDsum.

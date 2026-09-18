# KIDS26 Team 12 — glmnet HRDsum pipeline

Replaces `scripts/train_baseline.R` + missing `R/model.R` + `scripts/predict_frozen.R`.
Follows `config/analysis_protocol.json` (seed 260910, reference_HRDsum, elastic net,
<=5000 features, <=5% missing, GBM/LGG locked, 20% calibration, 95% conformal).

| File | Role |
|---|---|
| `R/model.R` | `fit_en()`, `predict_en()`, `metrics()`, `conformal_q()` |
| `01_prepare_features.R` | Streams beta.tsv, validates, keeps top 50k variable probes (label-free), dev only |
| `02_train_loco.R` | LOCO over non-CNS cancers + frozen model; parallel, resumable |
| `03_predict_frozen.R` | Streaming prediction from a frozen model |
| `run_pipeline.sh` | Plain server: `SMOKE=1`, `SUBSET=1`, or full |
| `submit_hpc.sh` | St. Jude LSF: `MODE=smoke`, `MODE=subset` (default), `MODE=full` |

## Prototype subset
`--cancers=BRCA,UCEC,BLCA,STAD,LUSC,HNSC` (default subset; override with `SUBSET_CANCERS=`).
GBM/LGG are refused. LOCO then runs over just those cancers (6 folds + frozen model).

## Blockers fixed
1. `R/model.R` missing -> written here.
2. Full 27 GB matrix loaded 3x (60+ GB RAM) -> streamed; model uses a 50k-probe matrix (~2.8 GB).
3. Training matrix copied per fold -> feature ranking from global sums minus held-out rows; only 5k columns copied.
4. Serial LOCO -> folds run in parallel (`--cores`), each saved on completion; rerun resumes.
5. QC annotation hard-coded exact match -> checked with clear error; override with `--expected_qc`.
6. Provenance path derived from input name -> explicit `--provenance`, no jsonlite needed.
7. predict_frozen loads whole matrix -> reads only model features.
8. Negative HRDsum predictions -> clipped at 0.

## Audit items from the original R/model.R header
- apply() over 336k columns (~42 h): stats come from global column sums; only <=5k columns copied per fold.
- 5-value hardcoded lambda grid, no boundary check: full data-derived glmnet path; `lambda_min_at_path_end` flagged in tuning_by_fold.tsv.
- No permuted-label or tissue-identity control: `--control=permute` and `--control=permute_within_cancer`.
- Inner folds: nested leave-one-cancer-out kept as the default (`--inner=cancer`), as in the original design.

## Lambda / alpha tuning (automated)
Per training fold: inner folds = one fold per training cancer (nested LOCO; `--inner=random` for stratified 5-fold), shared across
alphas {0.1,0.25,0.5,0.75,1}; `cv.glmnet` over the full lambda path; best alpha =
lowest CV MSE; lambda = `lambda.1se` (use `--lambda_rule=min` to compare).
Warnings are recorded if lambda.min hits the path end (`tuning_by_fold.tsv`).

## Outputs (`results/run01/`)
`loco_predictions.tsv`, `loco_metrics.tsv`, `macro_metrics.txt`, `tuning_by_fold.tsv`,
`frozen_nonCNS.rds`, `frozen_coefficients.tsv`, `final_partitions.tsv`, `sessionInfo.txt`.

## Known limits
- 50k prefilter uses variance across all dev samples incl. the held-out cancer (label-free).
- OOD flag (>25% of model features outside training 0.5–99.5% range) is a heuristic.
- Do not run 03 on GBM/LGG until development results are reviewed and locked.

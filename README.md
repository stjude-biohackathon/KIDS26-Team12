# KIDS26 - methylation prediction of genomic-scar burden

**Biohackathon overview:** [Project profile, team, mission, and three-day milestones](README_BIOHACKATHON.md)

**Team coordination:** [Project plan](project-management/project-plan.md) · [Team and roles](project-management/team.md) · [Team 12 Slack](https://stjudebiohackathon.slack.com/archives/C0BSC28M3U6)

Can a frozen methylation model predict independently measured HRDsum in a cancer type it has not seen, ultimately including pediatric high-grade glioma?

**Primary approach:** technically harmonized autosomal CpGs -> elastic-net regression -> predicted reference HRDsum. Canonical HRDsum=HRD-LOH+LST+TAI needs allele-aware genomic information; methylation total-CN profiles cannot directly identify LOH/TAI. CNV-only and fusion predictors are gated secondary comparisons. The prototype is a research scar-burden predictor, not a validated functional-HRD or therapy-selection test.

## Start today

Run from this cloned **KIDS26-Team12** directory using Python >=3.10:

```sh
python scripts/acquire_tcga.py published --tier labels --download
python scripts/acquire_tcga.py discover --kind beta --projects TCGA-BRCA TCGA-OV --limit 12 --download
python scripts/build_master.py
python scripts/prepare_beta.py --smoke-probes 1000
python -m unittest discover -s tests -p "test_*.py"
```

For the main historical cohort (41.5 GB; plan disk/RAM first):

```sh
python scripts/acquire_tcga.py published --tier beta --download
python scripts/index_publication.py
python scripts/build_master.py --metadata config/published_beta_metadata.tsv
python scripts/prepare_beta.py --published-matrix data/raw/pancanatlas/jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv --probes config/shared_autosomal_probes.txt
```

`shared_autosomal_probes.txt` now contains the deterministic 384,640-probe HM450/EPIC-v1 candidate bridge (SHA-256 `f359e43bc171c544e55e60000905bb0bba1c44ee0b37f3144c8ec733e052e5fc`). Confirm the PBTP EPIC generation before freezing it. The historical indexer reads only the header; subsequent beta extraction streams rows. Current GDC and historical betas have different preprocessing and must remain separate tracks until validated.

[Download guide](docs/DATA_ACQUISITION_NOW.md) covers tiers, verified resources, resume, checksums and manual fallbacks. All data/results are ignored by Git. Protected PBTP data belong in approved external storage with random project IDs and an internal crosswalk.

## Fit and infer after data/QC gates pass

```sh
Rscript scripts/setup.R
Rscript tests/smoke_model.R
Rscript -e "renv::snapshot(prompt=FALSE)"
Rscript scripts/train_baseline.R data/processed/beta.tsv data/processed/master_samples.tsv results/baseline
```

The R baseline rejects the small engineering fixture. It implements nested train-fold preprocessing, development LOCO excluding GBM/LGG, training-only mean nulls, and a final independent calibration reservation. Review development results and freeze the model/configuration before:

```sh
Rscript scripts/predict_frozen.R results/baseline/frozen_nonCNS.rds data/processed/locked_test_beta.tsv results/locked_predictions.tsv
Rscript -e "shiny::runApp('app')"
```

The demo defaults to prominently labeled synthetic fixtures. To use approved outputs, set `KIDS26_DEMO_RESULTS` to a precomputed TSV with the documented schema, including a provenance column. No PBTP upload or publication is automated.

The Shiny app now replaces the old methylation page with an **HRD Scores** page that loads `front_end_data/ddr_scars/` through `R/adapters/adapt_ddr_scores.R`. App prerequisites are installed by `Rscript scripts/setup.R`, including `VizModules`, `plotly`, `DT`, and `shinydashboard`.

## What exists now

- Repository/literature inventory and 70-row evidence matrix: 35 previously reviewed local PDFs, 25 newly reconciled PDFs (including one extra paper), 60 local PDFs total, and 10 unavailable entries with explicit gaps.
- Actionable public acquisition, publication-header indexing, conservative matching, beta preparation, total-CN proxy features and engineering simulator.
- R nested elastic-net/frozen inference implementation and Shiny scaffold; R runtime tests **now pass on the cluster** (2026-09-17).
- Four public label/metadata files and six beta files downloaded/verified; six uniquely matched BRCA specimens; 1000-probe engineering matrix.
- The team reports the full historical download complete outside the visible repository tree; confirm its 41,541,692,788-byte payload path and checksum before extraction. Only the header fixture was found under KIDS26 on 2026-09-13.
- 17 Python tests passed. Authoritative GDC project resolution and corrected primary-only duplicate logic retain 7,707 metadata-eligible historical columns across 32 cancers; the former 4,397 result is invalid. Actual interrupted download/resume and idempotent rerun were verified. Full HRD file has 10,647 unique sample IDs and no component-sum mismatches.
- **A trained biological model now exists (2026-09-17).** LOCO run 01 completed all 30 folds: within-tissue Pearson **0.612** (permutation p = 0.001, 29/30 tissues), but skill over the tissue-mean null is only **+0.085**. The model is a **within-tissue relative ranker**, not an absolute HRD calculator — it mis-levels each tissue by 3.46 units on average. Full results and limitations: [results](docs/24_RESULTS_LOCO_RUN01.md).
- Still absent: raw-IDAT preprocessing bridge, full CNV benchmark, calibrated HRD-high probability, any CNS or PBTP result. **The locked CNS partition (GBM+LGG, 642 samples) remains unopened.** See [verification status](docs/14_VERIFICATION_STATUS.md).

## Current status (2026-09-17)

| Item | State |
|---|---|
| LOCO run 01 (array 323078995) | **Complete**, 30/30 folds, merged 06:20 |
| Primary estimand E2 (within-tissue) | **Holds**, r = 0.612, p = 0.001 |
| Secondary estimand E1 (between-tissue) | Weak, skill +0.085 |
| Open defects gating the CNS unlock | C1 calibration, C2 zero floor, **C3 purity inversion**, C4 OV n=10 |
| CNS lock | **Closed** — opens once, after C1–C4 |

Read [the pipeline schematic](docs/23_PIPELINE_SCHEMATIC.md) for where each stage stands.

## Project map

| Area | Contents |
|---|---|
| docs/ | [Executive critique](docs/00_EXECUTIVE_PLAN.md), [literature](docs/01_LITERATURE_EVIDENCE.md), [new-paper reconciliation](docs/15_NEW_LITERATURE_RECONCILIATION.md), [method matrix](docs/02_METHOD_DECISION.md), [ground truth/CNV](docs/10_GROUND_TRUTH_AND_CNV.md) |
| docs/ | [Pre-event checklist](docs/04_PRE_HACKATHON_CHECKLIST.md), [runbook](docs/05_THREE_DAY_HACKATHON_RUNBOOK.md), [feature bridge](docs/17_450K_EPIC_FEATURE_BRIDGE.md), [PBTP readiness](docs/18_PBTP_READINESS.md), [validation](docs/06_VALIDATION_STRATEGY.md), [risk register](docs/11_RISK_REGISTER.md) |
| docs/ | [Synthetic data](docs/07_SYNTHETIC_DATA.md), [n-of-1](docs/08_N_OF_1_ROADMAP.md), [post-event](docs/09_POST_HACKATHON_PLAN.md), [environment](docs/12_ENVIRONMENT.md) |
| docs/ | [Blocker resolution plan](docs/21_BLOCKER_RESOLUTION_PLAN.md), [confounding and pan-cancer validity](docs/22_CONFOUNDING_AND_PANCANCER_PLAN.md) |
| config/ | Public manifests, curated evidence records, analysis protocol |
| scripts/, R/ | Acquisition/preparation/model/engineering utilities |
| data/raw, data/interim, data/processed | Ignored input, working and prepared data |
| results/ | Ignored tables/models/figures/run provenance |
| app/, tests/ | Approved-output demo and smoke tests |
| project-management/ | Original team/checklist files retained |

Prior root plans and papers remain unchanged outside this clone. The original template README is preserved in `docs/ORIGINAL_TEMPLATE_README.md`; existing contributor metadata/license/team organization remain in place. No Git commit or push was made.

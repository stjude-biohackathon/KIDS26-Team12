# KIDS26 - methylation prediction of genomic-scar burden

**Biohackathon overview:** [Project profile, team, mission, and three-day milestones](README_BIOHACKATHON.md)

**Team coordination:** [Project plan](project-management/project-plan.md) · [Team and roles](project-management/team.md) · [Team 12 Slack](https://stjudebiohackathon.slack.com/archives/C0BSC28M3U6)

Can a frozen methylation model predict independently measured HRDsum in a cancer type it has not seen, ultimately including pediatric high-grade glioma?

**Primary approach:** technically harmonized autosomal CpGs -> elastic-net regression -> predicted reference HRDsum. Canonical HRDsum=HRD-LOH+LST+TAI needs allele-aware genomic information; methylation total-CN profiles cannot directly identify LOH/TAI. CNV-only and fusion predictors are gated secondary comparisons. The prototype is a research scar-burden predictor, not a validated functional-HRD or therapy-selection test.

## The answer (2026-09-18): rank transfers weakly, absolute calibration does not

**The locked CNS partition was opened exactly once on 2026-09-18** (LSF `323264626`,
ledger row in `results/LOCKED_EVALUATION_LEDGER.tsv`). 642 adult TCGA samples,
135 GBM + 507 LGG, never touched during any fitting, tuning or model-selection
decision. This is the project's external domain test and it is now consumed.

Zero-shot result for the shipped model **V3-abs**:

| Group | n | Pearson [95% CI] | Spearman | MAE | bias | calib slope | skill vs ORACLE tissue-mean null |
|---|---|---|---|---|---|---|---|
| GBM | 135 | 0.320 [0.173, 0.461] | 0.279 | 8.86 | +7.99 | 0.358 | **−1.009** |
| LGG | 507 | 0.329 [0.234, 0.421] | 0.297 | 5.13 | +1.62 | 0.578 | **−0.048** |
| CNS pooled | 642 | 0.261 | 0.242 | 5.92 | +2.96 | 0.358 | **−0.233** |

Both lineages rank above chance; **neither beats an oracle that simply knows the
lineage's own mean HRDsum**. CNS is a poor relation rather than a typical
held-out tissue: against the 30 source LOCO tissues, GBM Pearson sits at the
13.3rd percentile (28th of 32) and LGG at the 16.7th (26th of 32), and GBM's
calibration offset is 29th of 32. The pre-declared sensitivity analysis
(Candidate A, run-01 pooled filter) gives GBM r = 0.293 / LGG r = 0.371 with
overlapping CIs — **neither candidate is clearly better in CNS**, and the
designation was fixed before unlock and cannot be revised (`docs/29`).

**The 42% reduction in excess lineage structure measured on the source tissues
did not translate into better CNS calibration.** That is the headline finding.
Verdict = outcome A of `docs/28` §7. **The deployment gate stays shut.**

Full record: [CNS protocol](docs/28_CANDIDATE_SELECTION_AND_CNS_PROTOCOL.md) ·
[pre-unlock audit](docs/29_PRE_UNLOCK_RECORD.md) ·
`results/cns_eval_2026-09-18/analysis_v3/cns_analysis_summary.md` (primary) and
`.../analysis_A/cns_analysis_summary.md` (sensitivity).

### Claims this project supports

- Methylation carries HRD-associated signal that **transfers in rank** to an
  unseen lineage: GBM and LGG Spearman are both positive with bootstrap 95% CI
  lower bounds above zero.
- Within the 30 development tissues the model is a **within-tissue relative
  ranker** of reference HRDsum, positive in 29 of 30 tissues.
- Ranking probes by pooled **within-tissue** variance measurably removes lineage
  imprinting from the predictor (42.0% of the excess, CI 37.9–46.5%).

### Claims this project does NOT support

- **No absolute zero-shot calibration.** In CNS an oracle tissue mean beats the
  model in both lineages (skill −1.009 GBM, −0.048 LGG).
- **No binary HRD call in CNS.** The exploratory HRDsum ≥ 42 cutoff is vacuous
  here: GBM has 0 samples above it and LGG exactly 1. No threshold claim is
  possible in this lineage.
- **No single-sample precision.** The frozen 95% conformal interval is ±21.10
  HRD units (Candidate A: ±21.73) — **wider than the label's own IQR of 4–28**.
  An individual prediction is not informative on its own.
- **No out-of-distribution warning.** The frozen reportability rule flagged only
  8 of 642 CNS samples (1.2%). It did **not** warn that an entire unseen lineage
  was out of distribution. This is a defect, recorded as such.
- No pediatric result, no PBTP data, no clinical validation, no HRD probability,
  no therapy-selection use.

### Don't be fooled by the pooled AUC

Both models were also scored as rankers at the exploratory HRDsum ≥ 42 cutoff on
**development** LOCO (`results/tables/roc_summary.tsv`, `roc_summary_v3.tsv`):

| Framing | Candidate A | V3-abs |
|---|---|---|
| Pooled AUC, HRDsum ≥ 42 | 0.862 [0.850, 0.875] | 0.866 [0.855, 0.878] |
| **Tissue-identity-ONLY control** | **0.743** | **0.723** |
| Within-tissue top quartile | 0.779 [0.767, 0.791] | 0.777 [0.764, 0.789] |

A pooled AUC near 0.86 looks strong until you notice that knowing nothing but
the cancer type already buys 0.72–0.74. The defensible number is the
tissue-centred one, ~0.78. The control falling 0.743 → 0.723 while the
within-tissue AUC holds (0.779 → 0.777) is the ROC view of the same lineage
reduction quantified in `results/tables/tableA2_A_vs_V3_full30.md`.

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

The R baseline rejects the small engineering fixture. It implements nested train-fold preprocessing, development LOCO excluding GBM/LGG, training-only mean nulls, and a final independent calibration reservation. Review development results and freeze the model/configuration before running inference:

```sh
Rscript scripts/predict_frozen.R <model.rds> <beta.tsv> <output.tsv> \
  [--allow-fixture] [--metadata=<master_samples.tsv>] [--locked-evaluation]

# non-locked scoring with the shipped V3-abs artifact:
Rscript scripts/predict_frozen.R results/frozen_v3_2026-09-18/frozen_nonCNS.rds \
  data/processed/beta.tsv results/predictions.tsv \
  --metadata=data/processed/master_samples.tsv

Rscript -e "shiny::runApp('app')"
```

`--metadata=` turns on the locked-partition guard and adds a `partition` column;
`--locked-evaluation` is the explicit opt-in required to score any `locked_CNS`
sample, and every such run appends an append-only row to
`results/LOCKED_EVALUATION_LEDGER.tsv`. **That path has already been used.** The
locked CNS cohort was scored once on 2026-09-18 with both frozen candidates (two
ledger rows, 20:05:08Z and 20:06:47Z). Re-running `--locked-evaluation` against
the same 642 samples does not produce a new held-out result; it produces a third
ledger row that documents the loss of held-out status. Do not do it.

The Shiny app now replaces the old methylation page with an **HRD Scores** page that loads `front_end_data/ddr_scars/` through `R/adapters/adapt_ddr_scores.R`. Launch it with `Rscript -e "shiny::runApp('app')"` after `Rscript scripts/setup.R` installs the app prerequisites (`VizModules`, `plotly`, `DT`, and `shinydashboard`). No PBTP upload or publication is automated.

## What exists now

- Repository/literature inventory and 70-row evidence matrix: 35 previously reviewed local PDFs, 25 newly reconciled PDFs (including one extra paper), 60 local PDFs total, and 10 unavailable entries with explicit gaps.
- Actionable public acquisition, publication-header indexing, conservative matching, beta preparation, total-CN proxy features and engineering simulator.
- R nested elastic-net/frozen inference implementation and Shiny scaffold; R runtime tests **now pass on the cluster** (2026-09-17).
- Four public label/metadata files and six beta files downloaded/verified; six uniquely matched BRCA specimens; 1000-probe engineering matrix.
- The team reports the full historical download complete outside the visible repository tree; confirm its 41,541,692,788-byte payload path and checksum before extraction. Only the header fixture was found under KIDS26 on 2026-09-13.
- 17 Python tests passed. Authoritative GDC project resolution and corrected primary-only duplicate logic retain 7,707 metadata-eligible historical columns across 32 cancers; the former 4,397 result is invalid. Actual interrupted download/resume and idempotent rerun were verified. Full HRD file has 10,647 unique sample IDs and no component-sum mismatches.
- **A trained biological model exists, and it has now been through an external domain test.**
  - **LOCO run 01 (2026-09-17, array `323078995`)** completed all 30 folds: within-tissue Pearson **0.612** (permutation p = 0.001, 29/30 tissues), skill over the tissue-mean null only **+0.085**, mean per-tissue mis-levelling 3.46 units. This is the historical record of the first run and remains correct as such: [results](docs/24_RESULTS_LOCO_RUN01.md).
  - **V3-abs (2026-09-18) is the SHIPPED model**, and supersedes run 01. It is the run-01 absolute model with **exactly one change**: the 5,000-probe unsupervised filter ranks probes by pooled **within-tissue** variance computed on training-fold rows only, instead of pooled total variance. It passed 7/7 pre-registered sentinel gates ([docs/27](docs/27_V3_PREREGISTRATION.md) §5) and 7/7 full-30-fold selection criteria S1–S7 ([docs/28](docs/28_CANDIDATE_SELECTION_AND_CNS_PROTOCOL.md) §4).
  - **What V3-abs improved** (full 30-fold, macro estimator, vs run 01 — `results/tables/tableA2_A_vs_V3_full30.md`): tissue R² of the predictions 0.5584 → 0.4672 against 0.3414 for the truth, i.e. `ΔR²_excess = 0.0913` [0.0824, 0.0996], p < 0.0001 = **42.0% of the excess lineage structure removed** (CI 37.9–46.5%); mean |tissue bias| 3.4628 → 2.9737; pooled MAE 8.9717 → 8.4456; macro skill vs the oracle tissue-mean null 0.419 → 0.630.
  - **What V3-abs did NOT improve:** per-tissue ranking. The paired change is indistinguishable from zero — Pearson mean −0.0161 (95% CI [−0.0410, 0.0051], p = 0.428), Spearman −0.0275 (p = 0.477), 12/30 and 14/30 tissues improved; macro within-tissue Pearson 0.5197 → 0.5036, Spearman 0.4747 → 0.4472; 29/30 tissues positive in both. **This is a targeted lineage reduction at no measurable ranking cost, not a general accuracy improvement.**
  - **Shipped artifact:** `results/frozen_v3_2026-09-18/frozen_nonCNS.rds`, sha256 `df9f7e82…3aaa4bf` (`df9f7e82b0b371b81ecca6a1d99a1a50128a95e78daa63307d5f3b0633a4aa72`), alpha 0.1, lambda 1.3598, 905 non-zero coefficients, 5,664 training + 1,401 calibration samples, 95% conformal q = **21.10**. Pre-declared sensitivity artifact (Candidate A, run-01 pooled filter): `results/frozen_2026-09-18/frozen_nonCNS.rds`, sha256 `0e975cd6…`, alpha 0.1, lambda 1.4201, 774 non-zero, conformal q = 21.73.
  - **V4 (lineage-penalized) FAILED and its branch was terminated.** It passed six of the seven sentinel gates but failed **G3, the mechanism gate**: tissue R² of the predictions went *up*, 0.5219 → 0.5361 against a ceiling of 0.4919. Penalising lineage explicitly made lineage imprinting worse. Recorded, not hidden (`results/v4_sentinel/gate_scorecard.tsv`).
  - **The locked CNS partition was OPENED and CONSUMED on 2026-09-18** (LSF `323264626`). See the section at the top of this README for the result.
- Still absent: raw-IDAT preprocessing bridge, full CNV benchmark, calibrated HRD-high probability, any pediatric or PBTP result, and any working out-of-distribution warning for an unseen lineage. See [verification status](docs/14_VERIFICATION_STATUS.md).

## Current status (2026-09-22; last experiment 2026-09-18)

| Item | State |
|---|---|
| Shipped model | **V3-abs**, frozen `results/frozen_v3_2026-09-18/frozen_nonCNS.rds` (sha256 `df9f7e82…`) |
| LOCO run 01 (array 323078995) | Complete 30/30; **superseded** as the shipped model, retained as record |
| V3-abs sentinel (323220006) / full 30-fold (323242800) | **7/7 gates**, **7/7 criteria S1–S7** |
| V4 lineage-penalized sentinel (323241956) | **FAILED G3** (tissue R² 0.5219 → 0.5361); branch terminated |
| Freezes (323234897 A, 323242801 B) | Both frozen from commit `28bb913` before unlock |
| Primary estimand E2 (within-tissue ranking) | **Holds** in development (29/30 tissues); transfers weakly to CNS (r ≈ 0.32–0.33) |
| Secondary estimand E1 (absolute/between-tissue) | **Fails** in CNS: oracle tissue-mean null beats the model in both lineages |
| **CNS lock** | **OPENED and CONSUMED 2026-09-18** (job 323264626, 642 samples); one-shot, ledger recorded |
| Deployment gate (`docs/27` §2A) | **SHUT** |
| Known open defects | C1 lineage-dependent offset (reduced 42%, not resolved), C3 purity inversion, C4 OV n=10, OOD rule did not flag CNS |

Read [the pipeline schematic](docs/23_PIPELINE_SCHEMATIC.md) for where each stage stands.

## Project map

| Area | Contents |
|---|---|
| docs/ | [Executive critique](docs/00_EXECUTIVE_PLAN.md), [literature](docs/01_LITERATURE_EVIDENCE.md), [new-paper reconciliation](docs/15_NEW_LITERATURE_RECONCILIATION.md), [method matrix](docs/02_METHOD_DECISION.md), [ground truth/CNV](docs/10_GROUND_TRUTH_AND_CNV.md) |
| docs/ | [Pre-event checklist](docs/04_PRE_HACKATHON_CHECKLIST.md), [runbook](docs/05_THREE_DAY_HACKATHON_RUNBOOK.md), [feature bridge](docs/17_450K_EPIC_FEATURE_BRIDGE.md), [PBTP readiness](docs/18_PBTP_READINESS.md), [validation](docs/06_VALIDATION_STRATEGY.md), [risk register](docs/11_RISK_REGISTER.md) |
| docs/ | [Synthetic data](docs/07_SYNTHETIC_DATA.md), [n-of-1](docs/08_N_OF_1_ROADMAP.md), [post-event](docs/09_POST_HACKATHON_PLAN.md), [environment](docs/12_ENVIRONMENT.md) |
| docs/ | [Blocker resolution plan](docs/21_BLOCKER_RESOLUTION_PLAN.md), [confounding and pan-cancer validity](docs/22_CONFOUNDING_AND_PANCANCER_PLAN.md) |
| docs/ | [Run 01 results](docs/24_RESULTS_LOCO_RUN01.md), [C1/C2/C3 investigation](docs/25_C1_C2_C3_INVESTIGATION.md), [C1 workaround](docs/26_C1_N_OF_1_WORKAROUND.md), [V3-abs pre-registration](docs/27_V3_PREREGISTRATION.md), [candidate selection and CNS protocol](docs/28_CANDIDATE_SELECTION_AND_CNS_PROTOCOL.md), [pre-unlock audit record](docs/29_PRE_UNLOCK_RECORD.md) |
| config/ | Public manifests, curated evidence records, analysis protocol |
| scripts/, R/ | Acquisition/preparation/model/engineering utilities |
| data/raw, data/interim, data/processed | Ignored input, working and prepared data |
| results/ | Ignored tables/models/figures/run provenance |
| app/, tests/ | Approved-output demo and smoke tests |
| project-management/ | Original team/checklist files retained |

Prior root plans and papers remain unchanged outside this clone. The original template README is preserved in `docs/ORIGINAL_TEMPLATE_README.md`; existing contributor metadata/license/team organization remain in place. Work is now tracked in Git commits in this clone; data, models and everything under `results/` remain ignored except the small committed summary tables and figures, and no protected-data publication has occurred.

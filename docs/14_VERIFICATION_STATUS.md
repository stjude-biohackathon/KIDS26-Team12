# Verification and handoff status - updated 2026-09-22 (last experiment 2026-09-18)

## 2026-09-13 readiness update

The team reports that data acquisition is complete. A recursive audit of the shared KIDS26 workspace found the verified label/development/annotation resources described below, but found no file larger than 1 GB and only the 280,257-byte historical-matrix header fixture. Before model-matrix extraction, record the actual 41,541,692,788-byte payload path and verify its checksum. This is a location/provenance gate, not a request to download a second copy.

## Executed successfully

- Python 3.12.14: all Python source files parsed; 17 unit tests passed. Tests cover specimen identity, primary-only duplicate handling, metadata fallback/conflicts, component order/sum, missing/excluded QC, balanced per-project acquisition selection, download corruption/path rejection, CN interval overlap, length-weighted features, centromere boundaries and total-CN non-identifiability of copy-neutral LOH.
- Four PanCanAtlas/PanImmune label/subtype/QC/purity files downloaded: 10,249,654 bytes. All sizes/MD5 verified; SHA256 provenance sidecars retained.
- Six current-GDC open 450K beta files downloaded: 78,530,106 bytes; checksums verified. They are six BRCA patients, not a representative cohort. Matching produced six eligible unique specimens and a 1,000-probe engineering matrix. Current beta format was observed to be headerless two-column data.
- Complete HRD table audited: 10,647 unique sample IDs; zero LOH+LST+TAI sum discrepancies.
- Actual interrupted download resumed from a 1,024-byte partial file and verified. A second run reported verified-existing. The downloader accommodates GDC's observed Content-Range header without the usual bytes prefix while checking the range bounds.
- Historical 450K matrix: bounded 1 MiB HTTP 206 read verified its real header, containing 9,664 sample columns. Only the header was retained as an ignored fixture. Indexing exposed and corrected mixed patient/sample subtype identifiers.
- Historical metadata resolution: all 8,868 unique matrix participants were batch-resolved through the GDC cases API to project IDs. Re-audit retains 7,707/9,664 columns across 32 cancer types, including 642 GBM/LGG locked-CNS candidates. Exclusions in precedence order: 1,363 non-primary, 319 missing/ambiguous HRD, 235 published quality exclusions, and 40 ambiguous multiple-primary candidates. Missing cancer type is now zero; authoritative GDC project metadata fills subtype-table gaps while true conflicts remain excluded. The prior 4,397/22-cancer count is invalid: it was driven by incomplete subtype coverage and a duplicate rule that allowed normal/recurrent files to disqualify a unique primary tumor. The new count is still metadata-level, before full beta extraction and array QC. See `docs/16_COHORT_QC.md`.
- Simulator ran: 17 truth segments and 2,160 noisy probe observations. Its simplified truth is an engineering fixture, not canonical scarHRD validation.
- Literature artifacts contain 70 entries, with 60 mapped local PDFs: 35 previously reviewed, 24 newly available entries from the original list, and one additional chordoma paper absent from the 69-row source. The 25-paper differential review found no exact duplicate or supplementary-only PDF and records explicit absent methods. See `docs/15_NEW_LITERATURE_RECONCILIATION.md`.
- HM450 and EPIC-v1 archived 2022-09 annotation manifests downloaded and SHA-256 verified (77,899,903 bytes combined). Deterministic intersection retains 384,640 shared autosomal `MASK_general=FALSE` CpGs with identical hg19 coordinates; allowlist SHA-256 is `f359e43bc171c544e55e60000905bb0bba1c44ee0b37f3144c8ec733e052e5fc`. PBTP EPIC generation remains to be confirmed.
- Public tier manifests generated without downloading the large tiers: historical beta 41,541,692,788 bytes; two CNV files 419,299,555 bytes combined. Current masked-intensity discovery identified 1,586 eligible files; a three-patient/six-file manifest totals 48,571,534 bytes. Intensities were not downloaded. A zero-result raw-intensity query is not proof of universal unavailability.

## Implemented but not executed here

*(Section heading retained for continuity. As of 2026-09-18 most of what follows
HAS executed on the cluster; see the 2026-09-18 execution update below.)*

R elastic-net nested LOCO fitting, R smoke tests, frozen inference and Shiny are implemented. Static inspection confirms outer-cancer isolation and training-only variance filtering, imputation, scaling, feature selection, alpha/lambda tuning, null mean and calibration reservation. The training script now exports final coefficients.

**UPDATED 2026-09-17 — the R runtime caveat below is resolved.** The previous
text read: *"R is not installed on Windows or the available WSL image ... so no R
runtime fit, package installation, lockfile, or rendered Shiny validation is
claimed."* R 4.5.0 is available on the cluster, `tests/smoke_model.R` passes
there (jobs 323078428, 323078739), and a full 30-fold nested LOCO fit has now
executed end-to-end (array `323078995`, 30/30 DONE, merged 2026-09-17 06:20).
Leakage assertions pass executably, not just by inspection.

What is now claimed on runtime evidence: nested LOCO fitting, training-only
preprocessing, lambda-path derivation (0/30 folds at a boundary), conformal
interval construction, and the inference provenance gate
(`tests/test_provenance_gate.R`). Results: `docs/24_RESULTS_LOCO_RUN01.md`.

What is still **not** claimed: rendered Shiny validation, a lockfile, any CNS or
pediatric result. The locked CNS partition remains unopened.

**SUPERSEDED 2026-09-18 — the CNS clause above is no longer true.** The locked
CNS partition was opened exactly once on 2026-09-18 under LSF `323264626`. The
pediatric/PBTP and Shiny/lockfile clauses still stand. See the 2026-09-18
execution update below.

The full historical matrix has not been located or streamed end-to-end by this repository. Probe-mask construction is complete for HM450 to EPIC v1, but exact PBTP platform confirmation and empirical preprocessing harmonization remain. Raw-IDAT QC/conumee2 execution, pinned canonical scarHRD comparison, broader platform/allele simulation and PBTP integration remain data-dependent.

**CORRECTED 2026-09-18.** The sentence that previously ended this paragraph —
*"No biological model has been trained; no transfer accuracy or calibrated
HRD-high probability exists."* — was already contradicted by the 2026-09-17
update above it, and is now wrong in two of its three clauses. The accurate
statement is: **a biological model has been trained, frozen and externally
tested.** No calibrated HRD-high probability exists, and none has been fitted or
validated; the continuous prediction is not a probability. Transfer accuracy now
exists as a measured quantity, and it is weak — see §2026-09-18 below.

## 2026-09-18 execution update — what has now actually run

Every item below executed on the cluster under LSF and left artifacts on disk.
Each is claimed on runtime evidence, not static inspection.

| Job | What it was | Outcome |
|---|---|---|
| `323220006` | V3-abs **sentinel**, 4 pre-registered tissues (BRCA, KICH, THCA, UCEC) | **7 of 7 gates passed**; gates frozen in `docs/27` §5 at 10:27 while `results/v3_sentinel/folds/` was verified empty. Scorecard `results/v3_sentinel/gate_scorecard.tsv` |
| `323242800` | V3-abs **full 30-fold LOCO array** | 30/30 folds complete; scored **7 of 7** against selection criteria S1–S7, frozen in `docs/28` §4 at 11:58 before the first fold returned. Outputs `results/v3_loco_full/` |
| `323241956` | **V4 lineage-penalized sentinel** | **FAILED G3, the mechanism gate**: tissue R² of predictions 0.5219 → **0.5361** against a ceiling of 0.4919, while G1, G2, G4, G5, G6, G7 all passed. Branch terminated; recorded not deleted. Scorecard `results/v4_sentinel/gate_scorecard.tsv` |
| `323234897` | **Freeze, Candidate A** (run 01, `feature_rank=pooled`) | `results/frozen_2026-09-18/frozen_nonCNS.rds`, sha256 `0e975cd6bc7caededffecf34332a3afd250e01bfea54fccab82266021f325763`, frozen 17:06:38 UTC from git commit `28bb913`. alpha 0.1, lambda 1.4201, 774 non-zero coefficients, conformal q = 21.73 |
| `323242801` | **Freeze, Candidate B** (V3-abs, `feature_rank=within_tissue`) | `results/frozen_v3_2026-09-18/frozen_nonCNS.rds`, sha256 `df9f7e82b0b371b81ecca6a1d99a1a50128a95e78daa63307d5f3b0633a4aa72`, frozen 17:53:28 UTC from git commit `28bb913`. alpha 0.1, lambda 1.3598, 905 non-zero coefficients, 5,664 training + 1,401 calibration samples, conformal q = 21.10. **This is the shipped model** |
| `323264626` | **Locked CNS evaluation**, executed exactly once | 642 samples (135 GBM, 507 LGG) scored with both frozen artifacts. Two append-only rows in `results/LOCKED_EVALUATION_LEDGER.tsv` (2026-09-18T20:05:08Z primary, 20:06:47Z sensitivity) |

Both freezes preceded the primary-candidate designation, and both preceded any
CNS label being read; the auditable pre-unlock record is
`docs/29_PRE_UNLOCK_RECORD.md`, committed while
`results/LOCKED_EVALUATION_LEDGER.tsv` did not yet exist.

### What the CNS evaluation established

Primary, V3-abs, zero-shot: GBM (n=135) Pearson 0.320 [0.173, 0.461], Spearman
0.279, MAE 8.86, bias +7.99, calibration slope 0.358; LGG (n=507) Pearson 0.329
[0.234, 0.421], Spearman 0.297, MAE 5.13, bias +1.62, slope 0.578; pooled
Pearson 0.261, MAE 5.92. Pre-declared sensitivity, Candidate A: GBM r 0.293,
LGG r 0.371, MAE 9.79 / 5.95, bias +9.11 / +3.69 — CIs overlap the primary
throughout and neither candidate is clearly better. Both are reported, as
`docs/28` §5 clause 3 requires.

Against the **oracle** within-lineage tissue-mean null, skill is **−1.009** in
GBM, **−0.048** in LGG, **−0.233** pooled: the model does not beat knowing the
lineage's own mean. Placed in the 30-tissue source LOCO distribution, GBM
Pearson sits at the 13.3rd percentile (28th of 32) and LGG at the 16.7th (26th
of 32), and GBM's |bias| ranks 29th of 32. CNS is a poor relation, not a typical
held-out tissue. Pre-declared outcome A of `docs/28` §7: **rank transfers
weakly, absolute cross-lineage calibration does not; the deployment gate stays
shut.** Full record: `results/cns_eval_2026-09-18/analysis_v3/cns_analysis_summary.md`
and `.../analysis_A/cns_analysis_summary.md`.

Note that the 42.0% reduction in excess lineage structure V3-abs achieved on the
source tissues (`ΔR²_excess = 0.0913`, CI [0.0824, 0.0996],
`results/tables/tableA2_A_vs_V3_full30.md`) **did not translate into better CNS
calibration.**

### Defects this evaluation exposed

- **The OOD reportability rule did not fire.** It flagged only 8 of 642 CNS
  samples (1.2%) and gave no indication that an entire unseen lineage was out of
  distribution. The frozen `ood_score`/`reportable` mechanism does not currently
  detect lineage shift. Recorded as a defect.
- **The conformal interval is not informative at single-sample resolution.**
  ±21.10 HRD units at 95% coverage for V3-abs (±21.73 for Candidate A) is wider
  than the label's own interquartile range of roughly 4 to 28.
- **HRDsum ≥ 42 is vacuous in CNS.** GBM has 0 samples above it and LGG exactly
  1, so no binary claim is possible in this lineage. Any AUC quoted at that
  cutoff in CNS rests on a single positive.

### Still not claimed, as of 2026-09-22

Rendered Shiny validation; a renv lockfile; any pediatric or PBTP result; any
calibrated HRD-high probability; any absolute zero-shot HRDsum claim in an
unseen lineage; any clinical or therapy-selection claim. The locked CNS
partition is now **consumed**: it has been scored and can never again serve as
an untouched held-out set. There is no second locked adult cohort.

## First reproducibility commands

Run inside KIDS26-Team12 with Python >=3.10 and R installed as needed:

```sh
python -m unittest discover -s tests -p "test_*.py"
python scripts/acquire_tcga.py published --tier beta --download
python scripts/resolve_tcga_cases.py
python scripts/index_publication.py --case-project-map config/published_case_projects.tsv
python scripts/build_master.py --metadata config/published_beta_metadata.tsv --out data/processed/historical
Rscript scripts/setup.R
Rscript tests/smoke_model.R
```

Inspect the historical matching audit and finalize the technical shared-probe allowlist before extracting model inputs. See the acquisition guide for the complete commands. Keep current-GDC and historical outputs separate.

## Scope and repository hygiene

Original parent plans/papers and team contributor/license files were retained. The template README is preserved verbatim as a historical artifact; its original relative links may no longer resolve from its archive location. The original docs/ai-guidance.md already refers to a nonexistent clone-level AGENTS.md; supplied project instructions reside in the parent directory. Newly authored handoff links were checked. Data, extracted text, models and large results remain ignored; `/results/` is gitignored apart from the small summary tables, figures, CNS analysis summaries and the locked-evaluation ledger that were added deliberately. Work is now committed to Git in this clone. No protected-data publication occurred.

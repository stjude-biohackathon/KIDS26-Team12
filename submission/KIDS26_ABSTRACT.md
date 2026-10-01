# KIDS26 Poster Abstract Submission — Team 12

## Title

**KIDS26-Team12: Can DNA Methylation Predict HRD Genomic-Scar Burden in Unseen Cancer Lineages? A Leakage-Safe Pan-Cancer Model with a Locked CNS Test**

## Authors

KIDS26 Biohackathon Team 12 — Evan Savage, Susanna Downing, Kayode Raheem, and team (see `project-management/team.md`; finalize order/affiliations before submission).

## Abstract (≈2,250 characters incl. spaces; limit 2,500)

**Background:** Homologous recombination deficiency (HRD) leaves genomic scars (HRDsum = HRD-LOH + LST + TAI) that are measured from allele-aware genomic data. DNA methylation arrays are already routine for pediatric CNS tumor classification, but whether they carry transferable HRD information is unknown.

**Objective:** To test whether a frozen, leakage-safe methylation model can predict independently measured HRDsum in a cancer lineage it has never seen, using adult glioma (GBM/LGG) as a proxy for pediatric high-grade glioma.

**Methods:** We assembled 7,707 TCGA primary tumors with 450K methylation and PanCanAtlas HRD labels, restricted to a 384,640-probe HM450/EPIC-v1 bridge. Elastic-net regression on 5,000 CpGs selected inside each training fold was evaluated by nested leave-one-cancer-out cross-validation across 30 non-CNS cancer types (n=7,065) against tissue-mean nulls and a tissue-identity permutation control. A pre-registered filter ranking CpGs by within-tissue variance (V3-abs) targeted lineage confounding. 642 GBM/LGG tumors were locked, then scored once with the frozen model under an append-only ledger. An R Shiny app enables interactive exploration of HRD scores.

**Results:** Within-tissue prediction was positive in 29/30 held-out tissues (Pearson 0.61, permutation p=0.001). Pooled AUC for HRDsum≥42 was 0.87, but tissue identity alone reached 0.72; the defensible within-tissue AUC was 0.78. V3-abs removed 42% (95% CI 38–47%) of excess lineage structure without measurable ranking loss. In the locked CNS cohort, ranking transferred weakly (GBM r=0.32 [0.17–0.46]; LGG r=0.33 [0.23–0.42]) but calibration failed (GBM bias +8.0, slope 0.36), and a tissue-mean oracle outperformed the model. Few-shot calibration with 10 labeled cases recovered 58% of the oracle gain.

**Conclusions:** Methylation carries HRD-associated signal that survives a lineage shift in rank but not in absolute scale, and reducing lineage imprinting in-distribution did not improve out-of-distribution calibration. Pediatric application will require a small labeled calibration set and a working out-of-distribution detector; this is a research prototype, not a clinical assay.

## Key numbers and sources

| Claim | Value | Source |
|---|---|---|
| Cohort | 7,707 (7,065 dev / 642 locked CNS) | `share/README.md`, `docs/30` §3 |
| Within-tissue Pearson (run 01) | 0.612, p=0.001, 29/30 | `docs/30` §6 |
| Pooled / control / within-tissue AUC | 0.866 / 0.723 / 0.777 | `docs/30` §9 |
| Lineage reduction | 42.0% (37.9–46.5%) | `results/tables/tableA2_A_vs_V3_full30.md` |
| CNS GBM / LGG Pearson | 0.320 / 0.329 | `results/tables/tableC_cns_external.md` |
| GBM bias / slope | +7.99 / 0.358 | `docs/30` §13 |
| Few-shot k=10 | 58% of oracle gain | `docs/30` §8 |

See `KIDS26_ABSTRACT_SUBMISSION_FULL.md` for the full submission document and figure.

# KIDS26 Poster Submission Package — Team 12

## 1. Title
**KIDS26-Team12: Can DNA Methylation Predict HRD Genomic-Scar Burden in Unseen Cancer Lineages? A Leakage-Safe Pan-Cancer Model with a Locked CNS Test**

Alternative (shorter): *KIDS26-Team12: Methyl-HRD — Cross-Lineage Prediction of Genomic-Scar Burden from DNA Methylation*

## 2. Abstract
See `KIDS26_ABSTRACT.md` (≈2,250 characters; structured Background / Objective / Methods / Results / Conclusions).

## 3. Figure 1 (single composite, placed under the abstract)

File: `submission/figures/KIDS26_Team12_Figure1.png` (built by `submission/build_figure1.R`). Layout is 3 × 2, landscape, 300 dpi.

| Panel | Content | Source file |
|---|---|---|
| A | Workflow: TCGA 450K + HRD labels → probe bridge → nested LOCO elastic net → frozen model → one-shot locked CNS test | `results/figures/fig01_workflow.png` |
| B | Within-tissue correlation across 30 held-out cancers (29/30 positive) | `results/figures/fig03_per_tissue_correlation.png` |
| C | Lineage imprinting: tissue R² of predictions vs truth, run 01 vs V3-abs (42% reduction) | `results/figures/fig05_lineage_imprinting.png` |
| D | ROC both ways: pooled vs tissue-identity control vs within-tissue | `results/figures/fig10_roc_both_ways_v3.png` |
| E | Locked CNS: predicted vs observed HRDsum, GBM and LGG | `results/figures/fig07_cns_pred_vs_obs.png` |
| F | Shiny app screenshot: HRD Scores page (scatter + component tabs) | `submission/figures/shiny_hrd_scores.png` (capture — see below) |

**Suggested legend:** *Figure 1. KIDS26-Team12 Methyl-HRD. (A) Leakage-safe pipeline with one-shot locked CNS evaluation. (B) Within-tissue Pearson r per held-out cancer type. (C) Within-tissue variance filtering removed 42% of excess lineage structure. (D) Pooled AUC is largely lineage recognition; within-tissue AUC = 0.78. (E) Zero-shot CNS predictions: rank transfers weakly (r≈0.32), calibration does not. (F) Interactive R Shiny explorer for HRD scores and components.*

### Capturing the Shiny screenshot (panel F)
```sh
Rscript -e "shiny::runApp('app', port = 8080)"
# open http://localhost:8080 → HRD Scores tab → Scatter Plot; screenshot at ≥1600 px wide
# save as submission/figures/shiny_hrd_scores.png
# or headless: Rscript -e "webshot2::appshot('app', 'submission/figures/shiny_hrd_scores.png', vwidth=1600, vheight=1000)"
```

## 4. Build the Word documents and figure
```sh
Rscript submission/build_figure1.R   # needs magick
Rscript submission/build_docx.R      # needs officer; writes .docx files
# or: pandoc submission/KIDS26_ABSTRACT.md -o submission/KIDS26_ABSTRACT.docx
```

## 5. Pre-submission checklist
- [ ] Confirm author list, order, affiliations
- [ ] Verify character count in the portal (spaces/headings may count)
- [ ] Capture Shiny screenshot and rebuild Figure 1
- [ ] Keep framing: research prototype, not a clinical HRD assay; adult CNS ≠ pediatric HGG

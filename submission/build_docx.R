# Build Word versions of the abstract and the full submission document.
# Usage: Rscript submission/build_docx.R
suppressPackageStartupMessages(library(officer))

title <- "KIDS26-Team12: Can DNA Methylation Predict HRD Genomic-Scar Burden in Unseen Cancer Lineages? A Leakage-Safe Pan-Cancer Model with a Locked CNS Test"
sections <- list(
  Background = "Homologous recombination deficiency (HRD) leaves genomic scars (HRDsum = HRD-LOH + LST + TAI) that are measured from allele-aware genomic data. DNA methylation arrays are already routine for pediatric CNS tumor classification, but whether they carry transferable HRD information is unknown.",
  Objective = "To test whether a frozen, leakage-safe methylation model can predict independently measured HRDsum in a cancer lineage it has never seen, using adult glioma (GBM/LGG) as a proxy for pediatric high-grade glioma.",
  Methods = "We assembled 7,707 TCGA primary tumors with 450K methylation and PanCanAtlas HRD labels, restricted to a 384,640-probe HM450/EPIC-v1 bridge. Elastic-net regression on 5,000 CpGs selected inside each training fold was evaluated by nested leave-one-cancer-out cross-validation across 30 non-CNS cancer types (n=7,065) against tissue-mean nulls and a tissue-identity permutation control. A pre-registered filter ranking CpGs by within-tissue variance (V3-abs) targeted lineage confounding. 642 GBM/LGG tumors were locked, then scored once with the frozen model under an append-only ledger. An R Shiny app enables interactive exploration of HRD scores.",
  Results = "Within-tissue prediction was positive in 29/30 held-out tissues (Pearson 0.61, permutation p=0.001). Pooled AUC for HRDsum>=42 was 0.87, but tissue identity alone reached 0.72; the defensible within-tissue AUC was 0.78. V3-abs removed 42% (95% CI 38-47%) of excess lineage structure without measurable ranking loss. In the locked CNS cohort, ranking transferred weakly (GBM r=0.32 [0.17-0.46]; LGG r=0.33 [0.23-0.42]) but calibration failed (GBM bias +8.0, slope 0.36), and a tissue-mean oracle outperformed the model. Few-shot calibration with 10 labeled cases recovered 58% of the oracle gain.",
  Conclusions = "Methylation carries HRD-associated signal that survives a lineage shift in rank but not in absolute scale, and reducing lineage imprinting in-distribution did not improve out-of-distribution calibration. Pediatric application will require a small labeled calibration set and a working out-of-distribution detector; this is a research prototype, not a clinical assay."
)
legend <- "Figure 1. KIDS26-Team12 Methyl-HRD. (A) Leakage-safe pipeline with one-shot locked CNS evaluation. (B) Within-tissue Pearson r per held-out cancer type. (C) Within-tissue variance filtering removed 42% of excess lineage structure. (D) Pooled AUC is largely lineage recognition; within-tissue AUC = 0.78. (E) Zero-shot CNS predictions: rank transfers weakly, calibration does not. (F) Interactive R Shiny explorer for HRD scores and components."

build <- function(path, with_figure) {
  doc <- read_docx()
  doc <- body_add_par(doc, title, style = "heading 1")
  doc <- body_add_par(doc, "KIDS26 Biohackathon Team 12", style = "Normal")
  for (nm in names(sections)) {
    doc <- body_add_fpar(doc, fpar(ftext(paste0(nm, ": "), fp_text(bold = TRUE)), ftext(sections[[nm]])))
  }
  chars <- nchar(paste(paste0(names(sections), ": ", unlist(sections)), collapse = " "))
  doc <- body_add_par(doc, sprintf("Character count (sections incl. headers): %d / 2500", chars), style = "Normal")
  fig <- "submission/figures/KIDS26_Team12_Figure1.png"
  if (with_figure && file.exists(fig)) {
    doc <- body_add_img(doc, fig, width = 6.5, height = 3.25)
    doc <- body_add_par(doc, legend, style = "Normal")
  }
  print(doc, target = path)
  message("Wrote ", path)
}
build("submission/KIDS26_ABSTRACT.docx", FALSE)
build("submission/KIDS26_ABSTRACT_SUBMISSION_FULL.docx", TRUE)

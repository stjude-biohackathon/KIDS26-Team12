# Build composite Figure 1 for the KIDS26-Team12 abstract submission.
# Usage: Rscript submission/build_figure1.R
suppressPackageStartupMessages(library(magick))

panels <- c(
  A = "results/figures/fig01_workflow.png",
  B = "results/figures/fig03_per_tissue_correlation.png",
  C = "results/figures/fig05_lineage_imprinting.png",
  D = "results/figures/fig10_roc_both_ways_v3.png",
  E = "results/figures/fig07_cns_pred_vs_obs.png",
  F = "submission/figures/shiny_hrd_scores.png"
)
cell_w <- 1200; cell_h <- 900

load_panel <- function(path, label) {
  img <- if (file.exists(path)) image_read(path) else
    image_annotate(image_blank(cell_w, cell_h, "grey95"),
                   paste("Missing:", path), size = 36, gravity = "center")
  img <- image_resize(img, sprintf("%dx%d", cell_w - 40, cell_h - 80))
  img <- image_extent(img, sprintf("%dx%d", cell_w, cell_h), color = "white", gravity = "center")
  image_annotate(img, label, size = 72, weight = 700, gravity = "northwest", location = "+15+5")
}

imgs <- mapply(load_panel, panels, names(panels), SIMPLIFY = FALSE)
row1 <- image_append(c(imgs$A, imgs$B, imgs$C))
row2 <- image_append(c(imgs$D, imgs$E, imgs$F))
fig  <- image_append(c(row1, row2), stack = TRUE)
fig  <- image_border(fig, "white", "20x20")

dir.create("submission/figures", showWarnings = FALSE, recursive = TRUE)
image_write(fig, "submission/figures/KIDS26_Team12_Figure1.png", density = 300)
message("Wrote submission/figures/KIDS26_Team12_Figure1.png")

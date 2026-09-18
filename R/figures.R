# =============================================================================
# R/figures.R -- shared presentation theme and helpers for KIDS26-Team12 figures
# -----------------------------------------------------------------------------
# PURPOSE
#   A single source of truth for figure look-and-feel, safe I/O, and the adopted
#   analysis conventions, so that every figure/table script renders consistently
#   and applies the same guard rails.
#
# INPUTS   : none (library only -- sourced by other scripts)
# OUTPUTS  : none (defines objects in the calling environment)
#
# PROVIDES
#   FIG_SEED                  integer seed used project-wide for figures/tables
#   PAL                       named palette (accent / warn / neutral / ...)
#   theme_kids26()            ggplot2 theme, base_size 13, clean, no chartjunk
#   scale_fill_kids26()       discrete fill scale over PAL
#   scale_colour_kids26()     discrete colour scale over PAL
#   read_required()           fread with a loud failure on a missing path
#   apply_c2_clip()           adopted C2 rule: pred_clip = pmax(pred, 0)
#   assert_development()      fails loudly if any GBM/LGG row is present
#   save_fig()                deterministic 300-dpi PNG writer
#   fmt()                     fixed-decimal number formatter for annotations
#
# PROJECT CONVENTIONS ENCODED HERE
#   * C2 (adopted): predictions are clipped at zero before ANY computation.
#   * CNS LOCK: the GBM/LGG cohort is sealed. Nothing in the figure/table layer
#     may read, score or plot it. assert_development() enforces this.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# --- seed --------------------------------------------------------------------
# Randomly drawn once from 1:10000 and then frozen for reproducibility.
FIG_SEED <- 4827L

# --- cohort lock -------------------------------------------------------------
CNS_LOCKED_TYPES <- c("GBM", "LGG")

# --- palette -----------------------------------------------------------------
# Deliberately small. Colour is used only where it carries meaning
# (a category, a sign, a highlighted subgroup) -- never as decoration.
PAL <- c(
  accent   = "#1F4E79",  # primary / model
  accent2  = "#3C8DAD",  # secondary series
  warn     = "#B02A37",  # failure mode, negative sign, flagged tissue
  good     = "#0F5132",  # positive sign, success
  amber    = "#D39E00",  # caution
  neutral  = "#6C757D",  # reference lines, de-emphasised text
  light    = "#CED4DA"
)

#' Discrete colour/fill scales drawn from PAL
scale_colour_kids26 <- function(...) {
  ggplot2::scale_colour_manual(values = unname(PAL[c("accent", "warn", "accent2",
                                                     "good", "amber", "neutral")]), ...)
}
scale_fill_kids26 <- function(...) {
  ggplot2::scale_fill_manual(values = unname(PAL[c("accent", "warn", "accent2",
                                                   "good", "amber", "neutral")]), ...)
}

# --- theme -------------------------------------------------------------------
#' Presentation theme: readable at projector distance, no chartjunk.
theme_kids26 <- function(base_size = 13, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = base_size * 1.25,
                                               margin = ggplot2::margin(b = 4)),
      plot.subtitle    = ggplot2::element_text(size = base_size * 0.92,
                                               colour = PAL[["neutral"]],
                                               margin = ggplot2::margin(b = 10)),
      plot.caption     = ggplot2::element_text(size = base_size * 0.72,
                                               colour = PAL[["neutral"]],
                                               hjust = 0,
                                               margin = ggplot2::margin(t = 10)),
      plot.title.position   = "plot",
      plot.caption.position = "plot",
      axis.title       = ggplot2::element_text(size = base_size * 0.95),
      axis.text        = ggplot2::element_text(size = base_size * 0.85, colour = "grey20"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
      panel.border     = ggplot2::element_blank(),
      axis.line        = ggplot2::element_line(colour = "grey30", linewidth = 0.4),
      legend.position  = "top",
      legend.justification = "left",
      legend.title     = ggplot2::element_text(size = base_size * 0.85),
      legend.text      = ggplot2::element_text(size = base_size * 0.85),
      legend.key.height = grid::unit(0.8, "lines"),
      strip.text       = ggplot2::element_text(face = "bold", size = base_size * 0.9),
      plot.margin      = ggplot2::margin(12, 16, 10, 12)
    )
}

# --- I/O ---------------------------------------------------------------------
#' fread() that fails loudly with an actionable message when the file is absent.
read_required <- function(path, ...) {
  if (!file.exists(path)) {
    stop("REQUIRED INPUT MISSING: '", path,
         "'\n  (cwd = ", normalizePath(".", mustWork = FALSE), ")", call. = FALSE)
  }
  dt <- data.table::fread(path, ...)
  if (nrow(dt) == 0L) {
    stop("REQUIRED INPUT EMPTY: '", path, "' has zero rows.", call. = FALSE)
  }
  message(sprintf("  read %-62s  %5d x %2d", path, nrow(dt), ncol(dt)))
  dt[]
}

#' Deterministic 300-dpi PNG writer. Returns the path invisibly.
save_fig <- function(plot, path, width, height, dpi = 300) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = path, plot = plot, width = width, height = height,
                  dpi = dpi, units = "in", bg = "white", device = "png")
  if (!file.exists(path)) stop("FAILED TO WRITE: ", path, call. = FALSE)
  message(sprintf("  wrote %-58s  %8.1f KB", path, file.size(path) / 1024))
  invisible(path)
}

# --- analysis conventions ----------------------------------------------------
#' Apply the adopted C2 rule. Adds `pred_clip` = pmax(prediction, 0).
apply_c2_clip <- function(dt, pred_col = "predicted_reference_HRDsum") {
  if (!pred_col %in% names(dt)) {
    stop("apply_c2_clip(): column '", pred_col, "' not found. Have: ",
         paste(names(dt), collapse = ", "), call. = FALSE)
  }
  dt <- data.table::copy(dt)
  dt[, pred_clip := pmax(get(pred_col), 0)]
  dt[]
}

#' Enforce the CNS lock: abort if GBM or LGG appear in a cancer-type column.
assert_development <- function(dt, type_col = "cancer_type",
                               context = deparse(substitute(dt))) {
  if (!type_col %in% names(dt)) {
    stop("assert_development(): no column '", type_col, "' in ", context, call. = FALSE)
  }
  hit <- intersect(CNS_LOCKED_TYPES, unique(as.character(dt[[type_col]])))
  if (length(hit) > 0L) {
    stop("CNS LOCK VIOLATION in ", context, ": locked cancer type(s) present -- ",
         paste(hit, collapse = ", "),
         ". The CNS cohort is sealed; figures must use development cancers only.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Filter explicitly to development cancers, then assert the lock held.
keep_development <- function(dt, type_col = "cancer_type",
                             context = deparse(substitute(dt))) {
  dt <- data.table::copy(dt)
  n0 <- nrow(dt)
  dt <- dt[!as.character(get(type_col)) %chin% CNS_LOCKED_TYPES]
  if (nrow(dt) != n0) {
    message(sprintf("  NOTE: dropped %d locked-CNS row(s) from %s", n0 - nrow(dt), context))
  }
  assert_development(dt, type_col, context)
  dt[]
}

# --- formatting --------------------------------------------------------------
fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)

# --- tissues flagged in the write-up as near-zero within-tissue signal --------
NEAR_ZERO_TISSUES <- c("THCA", "PCPG", "CHOL", "TGCT")

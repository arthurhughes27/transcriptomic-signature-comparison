# =============================================================================
# Presentation figure — aggregate-level specification robustness heatmap
# =============================================================================
# A slide-legible version of the main-text aggregate-level heatmap
# (analysis/specification_analysis/04_specification_heatmap_aggregate.R,
# R/robustness_heatmaps.R::plot_robustness_heatmap_aggregate()), built by
# layering extra theme() calls onto the existing plot rather than
# duplicating/modifying that function - the main-text figure is untouched.
#
# Kept at the same aspect ratio in spirit (width fixed at the main
# figure's 35cm) but taller (height increased - this is a standalone slide
# image, not constrained by a document page the way the main figure is),
# so the bigger axis titles / "Day X" strip text / y-axis (aggregate name)
# text have room without needing to shrink anything else.
#
# The x-axis (vaccine name) labels are the trickiest piece: bigger text at
# the main figure's existing 45-degree angle would need MORE width (the
# horizontal footprint of angled text grows with both font size and
# string length), but width is fixed here. Two changes address this
# without touching width:
#   - the x-axis labels are rotated to 90 degrees (vertical) - a vertical
#     label's horizontal footprint is just the font's character width,
#     independent of the (fairly long) vaccine name strings or how big the
#     font is, unlike angled text;
#   - the legend title is wrapped onto 3 lines (via labs(fill = ...), not
#     by changing the underlying scale) - a long one-line title at a
#     bigger font would force a wider legend column, eating into the
#     panel/x-axis-label width that's meant to stay fixed.
# XAXIS_TEXT_ANGLE below is kept as an easy-to-tweak constant if a
# different angle (e.g. 60-70 degrees) reads better once you see the PDF.
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────

library(fs)
library(tidyverse)

source(fs::path("R", "load_all.R"))

# ── Paths ─────────────────────────────────────────────────────────────────────

p_data_btm      <- fs::path("data", "BTM_processed.rds")
p_metrics       <- fs::path("output", "results", "specification_analysis", "robustness_metrics.rds")
out_dir         <- fs::path("output", "figures", "presentation")
p_fig_aggregate <- fs::path(out_dir, "robustness_heatmap_aggregate_presentation.pdf")

if (!fs::file_exists(p_metrics)) {
  stop(
    "Robustness metrics not found at ", p_metrics, ".\n",
    "Run 01_build_specification_grid.R, 02_run_raw_specifications.R and ",
    "03_apply_posthoc_and_robustness.R first."
  )
}

fs::dir_create(out_dir)

# ── Config ────────────────────────────────────────────────────────────────────

FIG_WIDTH_CM  <- 50  # same as the main-text figure - fixed per the request
FIG_HEIGHT_CM <- 32  # vs. the main-text figure's 20cm - see module comment

XAXIS_TEXT_ANGLE <- 45  # try e.g. 60-70 for a part-angled alternative

# ── Load data ─────────────────────────────────────────────────────────────────

BTM                  <- readRDS(p_data_btm)
robustness_metrics   <- readRDS(p_metrics)
robustness_annotated <- join_geneset_aggregates(robustness_metrics, BTM)

# =============================================================================
# PRESENTATION FIGURE — aggregate-level heatmap, slide-sized text
# =============================================================================

p_aggregate_presentation <- plot_robustness_heatmap_aggregate(robustness_annotated) +
  labs(fill = "Mean\nsignal-\nrobustness") +
  theme(
    plot.title    = element_text(size = 40, face = "bold"),
    axis.title    = element_text(size = 40),
    axis.text.y   = element_text(size = 20),
    axis.text.x   = element_text(
      size = 20, angle = XAXIS_TEXT_ANGLE,
      hjust = 1,
      vjust = if (XAXIS_TEXT_ANGLE == 90) 0.5 else 1
    ),
    strip.text    = element_text(size = 30, face = "bold"),
    legend.title  = element_text(size = 30, hjust = 0.5),
    legend.text   = element_text(size = 20),
    legend.key.height = unit(1.3, "cm")
  )

print(p_aggregate_presentation)

ggsave(
  filename = p_fig_aggregate, plot = p_aggregate_presentation,
  width = FIG_WIDTH_CM, height = FIG_HEIGHT_CM, units = "cm"
)
message("Saved presentation aggregate-level robustness heatmap to: ", p_fig_aggregate)

# rm(list = ls())

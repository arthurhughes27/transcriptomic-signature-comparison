# =============================================================================
# Specification analysis — 06: robustness distribution violin plot (Figure 4)
# =============================================================================
# Shows the full distribution of robustness across the 258 gene sets for
# each comparison (output/results/specification_analysis/
# robustness_metrics.rds), laid out exactly as the heatmaps in
# 04_specification_heatmaps.R: comparisons split into facets by timepoint,
# ordered within each facet by vaccine.
#
# Safe to run against partial robustness_metrics.rds results (e.g. while
# 02_run_raw_specifications.R / 03_apply_posthoc_and_robustness.R are still
# producing more of the specification grid).
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────

library(fs)
library(tidyverse)

source(fs::path("R", "load_all.R"))

# ── Paths ─────────────────────────────────────────────────────────────────────

p_metrics <- fs::path("output", "results", "specification_analysis", "robustness_metrics.rds")
out_dir   <- fs::path("output", "figures", "specification_analysis")
p_fig     <- fs::path(out_dir, "robustness_violin.pdf")

if (!fs::file_exists(p_metrics)) {
  stop(
    "Robustness metrics not found at ", p_metrics, ".\n",
    "Run 01_build_specification_grid.R, 02_run_raw_specifications.R and ",
    "03_apply_posthoc_and_robustness.R first."
  )
}

fs::dir_create(out_dir)

# ── Load data ─────────────────────────────────────────────────────────────────

robustness_metrics <- readRDS(p_metrics)

# =============================================================================
# FIGURE 4 — robustness distribution violin plot
# =============================================================================

p_violin <- plot_robustness_violin(robustness_metrics)

ggsave(filename = p_fig, plot = p_violin, width = 35, height = 15, units = "cm")
message("Saved robustness distribution violin plot to: ", p_fig)

rm(list = ls())

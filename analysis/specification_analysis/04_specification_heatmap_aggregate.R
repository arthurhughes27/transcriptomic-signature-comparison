# =============================================================================
# Specification analysis — 04: robustness heatmap (aggregate level, main figure)
# =============================================================================
# Aggregate-level visual summary of the robustness metric pi_{g,v,j}
# produced by 03_apply_posthoc_and_robustness.R (output/results/
# specification_analysis/robustness_metrics.rds) - see
# R/robustness_heatmaps.R for the plotting function itself. This is the
# main-text figure; the gene-set-level companion (one row per gene set,
# too long to fit as a main-text figure) is a supplementary figure - see
# analysis/supplementary/specification_heatmap_genesets_supplementary.R.
#
# Safe to run against partial results (e.g. while
# 02_run_raw_specifications.R / 03_apply_posthoc_and_robustness.R are still
# producing more of the specification grid) - it just plots whatever is in
# robustness_metrics.rds at the time.
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────

library(fs)
library(tidyverse)

source(fs::path("R", "load_all.R"))

# ── Paths ─────────────────────────────────────────────────────────────────────

p_data_btm      <- fs::path("data", "BTM_processed.rds")
p_metrics       <- fs::path("output", "results", "specification_analysis", "robustness_metrics.rds")
out_dir         <- fs::path("output", "figures", "specification_analysis")
p_fig_aggregate <- fs::path(out_dir, "robustness_heatmap_aggregate.pdf")

if (!fs::file_exists(p_metrics)) {
  stop(
    "Robustness metrics not found at ", p_metrics, ".\n",
    "Run 01_build_specification_grid.R, 02_run_raw_specifications.R and ",
    "03_apply_posthoc_and_robustness.R first."
  )
}

fs::dir_create(out_dir)

# ── Load data ─────────────────────────────────────────────────────────────────

BTM                  <- readRDS(p_data_btm)
robustness_metrics   <- readRDS(p_metrics)
robustness_annotated <- join_geneset_aggregates(robustness_metrics, BTM)

# =============================================================================
# MAIN FIGURE — aggregate-level heatmap
# =============================================================================

p_aggregate <- plot_robustness_heatmap_aggregate(robustness_annotated)

ggsave(
  filename = p_fig_aggregate, plot = p_aggregate,
  width = 35, height = 20, units = "cm"
)
message("Saved aggregate-level robustness heatmap to: ", p_fig_aggregate)

rm(list = ls())

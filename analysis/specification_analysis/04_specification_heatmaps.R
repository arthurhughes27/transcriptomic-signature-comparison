# =============================================================================
# Specification analysis — 04: robustness heatmaps
# =============================================================================
# Two-tier visual summary of the robustness metric pi_{g,v,j} produced by
# 03_apply_posthoc_and_robustness.R (output/results/specification_analysis/
# robustness_metrics.rds): an aggregate-level heatmap (main figure) and a
# gene-set-level heatmap (supplementary figure) - see R/robustness_heatmaps.R
# for the plotting functions themselves.
#
# With 258 gene sets, the gene-set-level heatmap is paginated (ROWS_PER_PAGE
# below) into a single multi-page PDF (save_multi_page_pdf(),
# R/plot_helpers.R) rather than one very tall page, so every page is a
# normal, consistently-sized page instead of requiring heavy zooming.
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
p_fig_genesets  <- fs::path(out_dir, "robustness_heatmap_genesets.pdf")

if (!fs::file_exists(p_metrics)) {
  stop(
    "Robustness metrics not found at ", p_metrics, ".\n",
    "Run 01_build_specification_grid.R, 02_run_raw_specifications.R and ",
    "03_apply_posthoc_and_robustness.R first."
  )
}

fs::dir_create(out_dir)

# ── Config ────────────────────────────────────────────────────────────────────

ROWS_PER_PAGE <- 45

# ── Load data ─────────────────────────────────────────────────────────────────

BTM                 <- readRDS(p_data_btm)
robustness_metrics  <- readRDS(p_metrics)
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

# =============================================================================
# SUPPLEMENTARY FIGURE — gene-set-level heatmap
# =============================================================================

p_genesets_pages <- plot_robustness_heatmap_genesets(robustness_annotated, rows_per_page = ROWS_PER_PAGE)

save_multi_page_pdf(
  p_genesets_pages, path = p_fig_genesets,
  width = 35 / 2.54, height = 20 / 2.54
)
message(sprintf(
  "Saved gene-set-level robustness heatmap (%d pages) to: %s",
  length(p_genesets_pages), p_fig_genesets
))

rm(list = ls())

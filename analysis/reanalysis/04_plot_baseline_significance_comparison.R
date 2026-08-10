# =============================================================================
# Baseline significance comparison: dearseq vs QuSAGE (Figure 5)
# =============================================================================
# Illustrates how the single largest analytical lever - DGSA method choice
# - moves results, using the baseline specification only: for each vaccine
# x timepoint comparison, the percentage of gene sets called significant
# at baseline (R/baseline_comparison_plots.R), dodged bars for dearseq vs
# QuSAGE, faceted by timepoint.
#
# Baseline results are read from the raw specification-grid outputs
# (R/baseline_results.R), so requires
# analysis/specification_analysis/01_build_specification_grid.R and
# 02_run_raw_specifications.R to have been run first.
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────

library(fs)
library(tidyverse)

source(fs::path("R", "load_all.R"))

# ── Paths ─────────────────────────────────────────────────────────────────────

p_data_btm     <- fs::path("data", "BTM_processed.rds")
spec_dir       <- fs::path("output", "results", "specification_analysis")
raw_dir        <- fs::path(spec_dir, "raw")
p_raw_grid     <- fs::path(spec_dir, "raw_specification_grid.rds")
figures_folder <- fs::path("output", "figures", "reanalysis")
p_fig          <- fs::path(figures_folder, "baseline_significance_comparison.pdf")

for (p_required in c(p_data_btm, p_raw_grid)) {
  if (!fs::file_exists(p_required)) {
    stop(
      "Required file not found: ", p_required, ".\n",
      "Run 01_build_specification_grid.R and 02_run_raw_specifications.R ",
      "(analysis/specification_analysis/) first."
    )
  }
}

fs::dir_create(figures_folder)

# ── Load baseline results ────────────────────────────────────────────────────

BTM              <- readRDS(p_data_btm)
raw_grid         <- readRDS(p_raw_grid)
baseline_results <- load_baseline_results_from_raw(
  raw_grid, raw_dir, BTM,
  conditions_order = default_conditions_order()
)

# =============================================================================
# FIGURE 5 — baseline significance, dearseq vs QuSAGE
# =============================================================================

pct_significant <- compute_pct_significant_by_comparison(baseline_results)

p_comparison <- plot_baseline_significance_comparison(pct_significant)

ggsave(filename = p_fig, plot = p_comparison, width = 35, height = 12, units = "cm")
message("Saved baseline significance comparison figure to: ", p_fig)

rm(list = ls())

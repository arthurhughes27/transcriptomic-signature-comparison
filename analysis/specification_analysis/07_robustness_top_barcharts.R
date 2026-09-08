# =============================================================================
# Specification analysis — 07: top-robustness gene set bar charts
# =============================================================================
# Visual companion to 06_robustness_tables.R's top-N table
# (build_top_robust_table()): rather than one global top-N ranking, shows
# each vaccine's own top N most robust gene sets (R/robustness_top_
# barcharts.R), separately for every timepoint. Exploratory - see the
# module-level comment in R/robustness_top_barcharts.R for the layout
# rationale.
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

# ── Config ────────────────────────────────────────────────────────────────────

TOP_N_PER_COMPARISON <- 5

# ── Paths ─────────────────────────────────────────────────────────────────────

p_data_btm  <- fs::path("data", "BTM_processed.rds")
p_metrics   <- fs::path("output", "results", "specification_analysis", "robustness_metrics.rds")
out_dir     <- fs::path("output", "figures", "specification_analysis")
p_fig       <- fs::path(out_dir, "robustness_top_barcharts.pdf")

if (!fs::file_exists(p_metrics)) {
  stop(
    "Robustness metrics not found at ", p_metrics, ".\n",
    "Run 01_build_specification_grid.R, 02_run_raw_specifications.R and ",
    "03_apply_posthoc_and_robustness.R first."
  )
}

fs::dir_create(out_dir)

# ── Load data ─────────────────────────────────────────────────────────────────

BTM                <- readRDS(p_data_btm)
robustness_metrics <- readRDS(p_metrics)

# =============================================================================
# TOP-N-PER-COMPARISON BAR CHARTS
# =============================================================================

p_barcharts <- plot_top_robust_barcharts(robustness_metrics, BTM, n = TOP_N_PER_COMPARISON)

print(p_barcharts)

# Taller than the earlier version: the y-axis headroom reserved for the
# now-vertical gene-set-name labels (plot_top_robust_barchart_day()'s
# scale_y_continuous(expand = ...)) means the 0-1 robustness bars only
# occupy a fraction of each panel's height, so the panels need to be much
# taller overall for the bars themselves to stay readable.
ggsave(
  filename = p_fig, plot = p_barcharts,
  width = 40, height = 80, units = "cm", limitsize = FALSE
)
message("Saved top-robustness gene set bar charts to: ", p_fig)

# rm(list = ls())

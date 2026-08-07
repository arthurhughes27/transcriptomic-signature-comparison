# =============================================================================
# Specification analysis — 05: robustness vs. baseline significance (Figure 3)
# =============================================================================
# Relates each gene set x vaccine x timepoint comparison's robustness metric
# pi_{g,v,j} (output/results/specification_analysis/robustness_metrics.rds)
# to that method's baseline-specification result, to check whether
# comparisons with strong baseline evidence also tend to be robust.
# Produces both figures described in R/robustness_baseline_plots.R:
# "binary" (baseline significant/not) and "continuous" (-log10 baseline
# p-value), faceted by method.
#
# Baseline results are read directly from the raw specification-grid
# outputs (R/baseline_results.R) - the same runs the robustness metric
# itself is accumulated from - rather than from separately-run
# analysis/reanalysis/dearseq_dgsa.R / qusage_dgsa.R results, which are not
# guaranteed to reproduce identically.
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

p_data_btm    <- fs::path("data", "BTM_processed.rds")
spec_dir      <- fs::path("output", "results", "specification_analysis")
raw_dir       <- fs::path(spec_dir, "raw")
p_raw_grid    <- fs::path(spec_dir, "raw_specification_grid.rds")
p_robustness  <- fs::path(spec_dir, "robustness_metrics.rds")
out_dir           <- fs::path("output", "figures", "specification_analysis")
p_fig_binary      <- fs::path(out_dir, "robustness_vs_baseline_binary.png")
p_fig_continuous  <- fs::path(out_dir, "robustness_vs_baseline_continuous.png")

for (p_required in c(p_data_btm, p_raw_grid, p_robustness)) {
  if (!fs::file_exists(p_required)) {
    stop(
      "Required file not found: ", p_required, ".\n",
      "Run the specification_analysis/ driver scripts (01-03) first."
    )
  }
}

fs::dir_create(out_dir)

# ── Load + join ───────────────────────────────────────────────────────────────

BTM                 <- readRDS(p_data_btm)
raw_grid            <- readRDS(p_raw_grid)
robustness_metrics  <- readRDS(p_robustness)
baseline_results    <- load_baseline_results_from_raw(raw_grid, raw_dir, BTM)

joined <- join_robustness_baseline(robustness_metrics, baseline_results)

# =============================================================================
# FIGURE 3 — binary and continuous modes
# =============================================================================

p_binary <- plot_robustness_vs_baseline(joined, mode = "binary")
ggsave(filename = p_fig_binary, plot = p_binary, width = 20, height = 12, units = "cm")
message("Saved binary-mode figure to: ", p_fig_binary)

p_continuous <- plot_robustness_vs_baseline(joined, mode = "continuous")
ggsave(filename = p_fig_continuous, plot = p_continuous, width = 20, height = 12, units = "cm")
message("Saved continuous-mode figure to: ", p_fig_continuous)

rm(list = ls())

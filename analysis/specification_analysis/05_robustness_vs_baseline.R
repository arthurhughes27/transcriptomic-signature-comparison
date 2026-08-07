# =============================================================================
# Specification analysis — 05: robustness vs. baseline significance (Figure 3)
# =============================================================================
# Relates each gene set x vaccine x timepoint comparison's robustness metric
# pi_{g,v,j} (output/results/specification_analysis/robustness_metrics.rds)
# to that method's baseline-specification result (Chapter 2, Section 2.1.3
# Stage 1's dearseq/QuSAGE baseline runs), to check whether comparisons with
# strong baseline evidence also tend to be robust. Produces both figures
# described in R/robustness_baseline_plots.R: "binary" (baseline
# significant/not) and "continuous" (-log10 baseline p-value), faceted by
# method.
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

p_robustness      <- fs::path("output", "results", "specification_analysis", "robustness_metrics.rds")
p_results_dearseq <- fs::path("output", "results", "reanalysis", "dearseq_dgsa_results_processed.rds")
p_results_qusage  <- fs::path("output", "results", "reanalysis", "qusage_dgsa_results_processed.rds")
out_dir           <- fs::path("output", "figures", "specification_analysis")
p_fig_binary      <- fs::path(out_dir, "robustness_vs_baseline_binary.png")
p_fig_continuous  <- fs::path(out_dir, "robustness_vs_baseline_continuous.png")

for (p_required in c(p_robustness, p_results_dearseq, p_results_qusage)) {
  if (!fs::file_exists(p_required)) {
    stop(
      "Required file not found: ", p_required, ".\n",
      "Run the specification_analysis/ and reanalysis/ driver scripts first."
    )
  }
}

fs::dir_create(out_dir)

# ── Load + join ───────────────────────────────────────────────────────────────

robustness_metrics <- readRDS(p_robustness)
baseline_results   <- dplyr::bind_rows(
  readRDS(p_results_dearseq),
  readRDS(p_results_qusage)
)

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

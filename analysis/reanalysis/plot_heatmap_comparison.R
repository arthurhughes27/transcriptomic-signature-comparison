# =============================================================================
# DGSA heatmap comparison: dearseq vs QuSAGE
# =============================================================================
# ComplexHeatmap comparison of the two methods' baseline results, analogous
# to the circos plots (plot_circos.R): mean fold-change per gene set x
# vaccine, split into column panels by timepoint, asterisk-annotated where
# significant, row strip coloured by gene-set aggregate.
#
# Cleaned up and modularised (R/heatmap_comparison.R) from a script carried
# over from an earlier repository. Two issues fixed during cleanup:
#   - the two plot_heatmap() calls feeding each "comparison" PDF used
#     IDENTICAL parameters (score_threshold = 7 for both methods, in both
#     the first and second parameter block), despite a comment claiming the
#     second was "a more relaxed sharing score threshold on QuSAGE" - and
#     both blocks wrote to the same filename, so the first PDF write was
#     always immediately overwritten by the second, functionally identical
#     one (i.e. it was dead computation, not a second real output). The
#     redundant first write is removed here; each heatmap is now built
#     once. If you did intend a relaxed QuSAGE threshold, change
#     SCORE_THRESHOLD_QUSAGE below.
#   - the single-panel "heatmap_qusage.pdf" block's overlaid title read
#     "dearseq" (copy-pasted from the block above it, drawing ht2/QuSAGE
#     under a "dearseq" label); fixed to "QuSAGE".
#
# Baseline results are read from the raw specification-grid outputs
# (R/baseline_results.R) - the same reproducible source used by every
# other baseline-dependent script (dgsa_comparison_example.R,
# plot_circos.R, plot_baseline_significance_comparison.R,
# 05_robustness_vs_baseline.R, 07_robustness_tables.R) - rather than the
# original script's output/results/{dearseq,qusage}/*_processed.rds paths,
# which don't exist in this repository's layout and are the separately-run
# baseline scripts already found not to reproduce identically.
#
# Requires analysis/specification_analysis/01_build_specification_grid.R
# and 02_run_raw_specifications.R to have been run first.
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(fs)
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

source(fs::path("R", "load_all.R"))

# ── Paths ─────────────────────────────────────────────────────────────────────

p_data_btm     <- fs::path("data", "BTM_processed.rds")
spec_dir       <- fs::path("output", "results", "specification_analysis")
raw_dir        <- fs::path(spec_dir, "raw")
p_raw_grid     <- fs::path(spec_dir, "raw_specification_grid.rds")
figures_folder <- fs::path("output", "figures", "reanalysis")

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

BTM      <- readRDS(p_data_btm)
raw_grid <- readRDS(p_raw_grid)

results_df <- load_baseline_results_from_raw(
  raw_grid, raw_dir, BTM,
  conditions_order = default_conditions_order()
) %>%
  assign_dgsa_colours(
    condition_colors = default_condition_colors(),
    aggregate_colors = default_aggregate_colors()
  ) %>%
  mutate(gs.name.description = abbreviate_geneset_label(gs.name.description))

# ── Config ────────────────────────────────────────────────────────────────────

ALL_CONDITIONS <- levels(results_df$condition)
ALL_TIMES      <- c(1, 3, 7, 10, 14, 21)
ALL_AGGREGATES <- setdiff(levels(results_df$gs.aggregate), "NA")

SCORE_THRESHOLD_DEARSEQ <- 7
SCORE_THRESHOLD_QUSAGE  <- 7  # set differently to relax the QuSAGE sharing-score filter (see header note)

common_heatmap_args <- list(
  conditions             = ALL_CONDITIONS,
  times                  = ALL_TIMES,
  aggregates             = ALL_AGGREGATES,
  fixed_row_names_width  = unit(150, "mm"),
  p_correction           = "BH",
  p_approach             = "global",
  p_threshold            = 0.05,
  scores                 = "fc.score",
  filter_mode            = "none",
  user_threshold         = 0.5,
  quantile_threshold     = 0.5,
  y_order                = "cluster",
  x_order                = "cluster",
  filter_commonDE        = "score",
  common_proportion      = 0,
  quantile_scoreclip     = 1,
  legend_max             = 1.5
)

# =============================================================================
# BUILD HEATMAPS (once each - see header note on the removed duplicate build)
# =============================================================================

ht_dearseq <- do.call(plot_dgsa_heatmap, c(
  list(results_df = results_df, method_name = "dearseq", score_threshold = SCORE_THRESHOLD_DEARSEQ),
  common_heatmap_args
))
ht_qusage <- do.call(plot_dgsa_heatmap, c(
  list(results_df = results_df, method_name = "qusage", score_threshold = SCORE_THRESHOLD_QUSAGE),
  common_heatmap_args
))

# =============================================================================
# COMBINED COMPARISON FIGURE
# =============================================================================

p_fig_comparison <- fs::path(figures_folder, "heatmap_comparison.pdf")
save_stacked_heatmaps_pdf(
  path     = p_fig_comparison,
  heatmaps = list(ht_dearseq, ht_qusage),
  titles   = c("dearseq", "QuSAGE"),
  heights  = list(unit(0.6, "npc"), unit(0.3, "npc")),
  width = 20, height = 25
)
message("Saved heatmap comparison figure to: ", p_fig_comparison)

# =============================================================================
# INDIVIDUAL METHOD FIGURES
# =============================================================================

p_fig_dearseq <- fs::path(figures_folder, "heatmap_dearseq.pdf")
save_stacked_heatmaps_pdf(
  path     = p_fig_dearseq,
  heatmaps = list(ht_dearseq),
  titles   = "dearseq",
  heights  = list(unit(1, "npc")),
  width = 28, height = 15
)
message("Saved dearseq heatmap figure to: ", p_fig_dearseq)

p_fig_qusage <- fs::path(figures_folder, "heatmap_qusage.pdf")
save_stacked_heatmaps_pdf(
  path     = p_fig_qusage,
  heatmaps = list(ht_qusage),
  titles   = "QuSAGE",
  heights  = list(unit(1, "npc")),
  width = 22.5, height = 10
)
message("Saved QuSAGE heatmap figure to: ", p_fig_qusage)

rm(list = ls())

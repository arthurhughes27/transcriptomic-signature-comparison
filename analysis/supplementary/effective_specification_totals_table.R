# =============================================================================
# Effective specification totals table (supplementary)
# =============================================================================
# LaTeX-ready table of the deduplicated ("effective") specification count
# per vaccine x timepoint comparison (R/specification_redundancy.R's
# specification_totals() - see the module-level comment there for why a
# comparison's raw specification count can be below the nominal 34: some
# comparisons have no variation in one or more covariates within their
# paired sample set - e.g. a single contributing study, or an all-"Unknown"
# gender/race - so dearseq's covariate design matrix collapses and what
# build_raw_specification_grid() counts as separate raw specifications
# produce byte-identical results for that comparison).
#
# This is a report of the SAME per-comparison totals
# 03_apply_posthoc_and_robustness.R already writes to
# output/results/specification_analysis/effective_specification_totals.csv
# - this script just re-derives them (cheap - no DGSA re-run, see
# R/specification_redundancy.R) and renders them as a LaTeX table, ordered
# by timepoint then vaccine rather than the CSV's (arbitrary) grouping
# order.
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────

library(fs)
library(tidyverse)
library(knitr)

source(fs::path("R", "load_all.R"))

# ── Paths ─────────────────────────────────────────────────────────────────────

p_data_expr <- fs::path("data", "hipc_merged_young_noNorm.rds")
spec_dir    <- fs::path("output", "results", "specification_analysis")
tables_dir  <- fs::path("output", "tables", "supplementary")

p_raw_grid     <- fs::path(spec_dir, "raw_specification_grid.rds")
p_posthoc_grid <- fs::path(spec_dir, "posthoc_specification_grid.rds")

for (p_required in c(p_data_expr, p_raw_grid, p_posthoc_grid)) {
  if (!fs::file_exists(p_required)) {
    stop(
      "Required file not found: ", p_required, ".\n",
      "Run 01_build_specification_grid.R (analysis/specification_analysis/) first."
    )
  }
}

fs::dir_create(tables_dir)

# Must match 02_run_raw_specifications.R/03_apply_posthoc_and_robustness.R's
# DAYS_TO_ANALYSE - the comparisons considered here must be exactly the
# ones actually run.
DAYS_TO_ANALYSE <- c(1, 3, 7)

# ── Load data + specification grids ─────────────────────────────────────────

raw_grid     <- readRDS(p_raw_grid)
posthoc_grid <- readRDS(p_posthoc_grid)

# Same vaccine_code/study_accession preprocessing 02_run_raw_specifications.R
# applies before calling run_dearseq_comparison() - build_covariate_matrix()
# needs study_accession as a factor to reproduce the real design matrices
# (see R/specification_redundancy.R).
hipc <- readRDS(p_data_expr) |>
  dplyr::mutate(
    vaccine_code    = as.factor(if_else(time_post_last_vax > 0, 2, 1)),
    study_accession = as.factor(study_accession)
  )

# =============================================================================
# BUILD AND SAVE THE TABLE
# =============================================================================

spec_totals_table <- specification_totals(hipc, raw_grid, posthoc_grid, days = DAYS_TO_ANALYSE) |>
  dplyr::mutate(vaccine_name = factor(vaccine_name, levels = default_conditions_order())) |>
  dplyr::arrange(day, vaccine_name) |>
  dplyr::transmute(
    `Comparison`                    = sprintf("%s (Day %s)", vaccine_name, day),
    `Effective total specifications` = n_effective_total_specs
  )

save_latex_table(
  spec_totals_table,
  path_tex = fs::path(tables_dir, "effective_specification_totals.tex"),
  path_csv = fs::path(tables_dir, "effective_specification_totals.csv"),
  caption  = paste(
    "Deduplicated (\"effective\") number of specifications considered for",
    "each vaccine x timepoint comparison. Some comparisons have no",
    "variation in one or more covariates within their paired sample set",
    "(e.g. a single contributing study), so a subset of the nominal 34 raw",
    "specifications produce identical results for that comparison and are",
    "counted once rather than separately; the total specification count is",
    "the resulting effective raw-specification count multiplied by the 81",
    "post-hoc specifications."
  ),
  label = "tab:effective-specification-totals"
)

message(sprintf(
  "Effective specification totals table (%d rows) saved to: %s",
  nrow(spec_totals_table), tables_dir
))

rm(list = ls())

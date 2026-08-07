# =============================================================================
# Specification analysis — 07: robustness summary tables
# =============================================================================
# Two LaTeX-ready tables summarising the robustness metric pi_{g,v,j}
# (output/results/specification_analysis/robustness_metrics.rds):
#   1. the TOP_N most robust gene set x comparison results.
#   2. gene set x comparison x method results that are significant at the
#      baseline specification but fall below ROBUSTNESS_THRESHOLD -
#      conclusions most fragile to analytical choices.
# Both use the full geneset.names.descriptions label (R/robustness_tables.R).
# Output written to output/tables/specification_analysis/ as .tex (ready to
# \input{} into the thesis) and .csv (for reference), matching the
# convention in analysis/descriptive/is2_appendix_descriptives.R.
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────

library(fs)
library(tidyverse)
library(knitr)

source(fs::path("R", "load_all.R"))

# ── Config ────────────────────────────────────────────────────────────────────

TOP_N                <- 20
ROBUSTNESS_THRESHOLD <- 0.5

# ── Paths ─────────────────────────────────────────────────────────────────────

p_robustness      <- fs::path("output", "results", "specification_analysis", "robustness_metrics.rds")
p_results_dearseq <- fs::path("output", "results", "reanalysis", "dearseq_dgsa_results_processed.rds")
p_results_qusage  <- fs::path("output", "results", "reanalysis", "qusage_dgsa_results_processed.rds")
p_data_btm        <- fs::path("data", "BTM_processed.rds")
tables_dir        <- fs::path("output", "tables", "specification_analysis")

for (p_required in c(p_robustness, p_results_dearseq, p_results_qusage, p_data_btm)) {
  if (!fs::file_exists(p_required)) {
    stop(
      "Required file not found: ", p_required, ".\n",
      "Run the specification_analysis/ and reanalysis/ driver scripts first."
    )
  }
}

fs::dir_create(tables_dir)

# ── Load data ─────────────────────────────────────────────────────────────────

BTM                 <- readRDS(p_data_btm)
robustness_metrics  <- readRDS(p_robustness)
baseline_results    <- dplyr::bind_rows(
  readRDS(p_results_dearseq),
  readRDS(p_results_qusage)
)

# =============================================================================
# TABLE 1 — top N most robust results
# =============================================================================

top_robust <- build_top_robust_table(robustness_metrics, BTM, n = TOP_N)

save_latex_table(
  top_robust,
  path_tex = fs::path(tables_dir, "top_robust_results.tex"),
  path_csv = fs::path(tables_dir, "top_robust_results.csv"),
  caption  = sprintf(
    "The %d most robust gene set x vaccine x timepoint results: the proportion of post-hoc specifications under which each result is called significant (robustness, pi).",
    TOP_N
  ),
  label = "tab:top-robust-results"
)
message(sprintf("Saved top-%d robust results table to: %s", TOP_N, tables_dir))

# =============================================================================
# TABLE 2 — significant at baseline but not robust
# =============================================================================

significant_nonrobust <- build_significant_nonrobust_table(
  robustness_metrics, baseline_results, BTM,
  robustness_threshold = ROBUSTNESS_THRESHOLD,
  n = 20
)

save_latex_table(
  significant_nonrobust,
  path_tex = fs::path(tables_dir, "significant_nonrobust_results.tex"),
  path_csv = fs::path(tables_dir, "significant_nonrobust_results.csv"),
  caption  = sprintf(
    "Gene set x vaccine x timepoint results significant at the baseline specification (global BH-adjusted p <= 0.05) but with robustness below %.2f: results whose significance is most sensitive to analytical choices.",
    ROBUSTNESS_THRESHOLD
  ),
  label = "tab:significant-nonrobust-results"
)
message(sprintf(
  "Saved significant-but-non-robust results table (%d rows) to: %s",
  nrow(significant_nonrobust), tables_dir
))

rm(list = ls())

# =============================================================================
# Stage 1: dearseq vs QuSAGE baseline comparison
# =============================================================================
# Compares the two DGSA methods' baseline-specification results (Chapter 2,
# Section 2.1.3, Stage 1: "reanalysing the IS2 data using an alternative
# DGSA method and comparing the results to the initial study"), using the
# concordance metrics in R/comparison_metrics.R: rank correlation of raw
# p-values and fold-change scores, and agreement on BH-adjusted
# significance calls (within-timepoint scope, alpha = 0.05), both overall
# and per vaccine x timepoint comparison.
#
# Baseline results are read from the raw specification-grid outputs
# (R/baseline_results.R), not from separately-run analysis/reanalysis/
# dearseq_dgsa.R / qusage_dgsa.R results: the latter are a distinct run of
# the same specification and were found not to reproduce identically
# (gene sets baseline-significant here but with zero robustness in the
# specification analysis, which is impossible if this were the same run).
# Requires 01_build_specification_grid.R and 02_run_raw_specifications.R
# (analysis/specification_analysis/) to have been run first.
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
p_summary_out <- fs::path("output", "results", "reanalysis", "dgsa_comparison_summary.rds")

for (p_required in c(p_data_btm, p_raw_grid)) {
  if (!fs::file_exists(p_required)) {
    stop(
      "Required file not found: ", p_required, ".\n",
      "Run 01_build_specification_grid.R and 02_run_raw_specifications.R ",
      "(analysis/specification_analysis/) first."
    )
  }
}

# ── Load baseline results ────────────────────────────────────────────────────

BTM       <- readRDS(p_data_btm)
raw_grid  <- readRDS(p_raw_grid)

results_df <- load_baseline_results_from_raw(raw_grid, raw_dir, BTM)

# =============================================================================
# CONCORDANCE METRICS
# =============================================================================

concordance <- compute_concordance_metrics(
  results_df,
  methods       = c("dearseq", "qusage"),
  adj_pval_col  = "withinTime.adjPval_BH",
  alpha         = 0.05,
  score_col     = "fc.score"
)

message(sprintf(
  "\nOverall (%d gene set x comparison pairs evaluable by both methods):\n",
  concordance$overall$n_gene_sets
))
print(concordance$overall)

message("\nBy comparison (sorted by lowest significance agreement first):\n")
print(
  concordance$by_comparison %>% arrange(pct_agree),
  n = Inf
)

# ── Save ───────────────────────────────────────────────────────────────────────

saveRDS(concordance, file = p_summary_out)
message("\nComparison summary saved to: ", p_summary_out)

rm(list = ls())

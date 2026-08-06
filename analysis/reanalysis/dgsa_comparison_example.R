# =============================================================================
# Stage 1: dearseq vs QuSAGE baseline comparison
# =============================================================================
# Compares the two DGSA methods' baseline-specification results (Chapter 2,
# Section 2.1.3, Stage 1: "reanalysing the IS2 data using an alternative
# DGSA method and comparing the results to the initial study"), using the
# concordance metrics in R/comparison_metrics.R: rank correlation of raw
# p-values and fold-change scores, and agreement on BH-adjusted
# significance calls (global scope, alpha = 0.05), both overall and per
# vaccine x timepoint comparison.
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────

library(fs)
library(tidyverse)

source(fs::path("R", "load_all.R"))

# ── Paths ─────────────────────────────────────────────────────────────────────

p_results_dearseq <- fs::path("output", "results", "reanalysis", "dearseq_dgsa_results_processed.rds")
p_results_qusage  <- fs::path("output", "results", "reanalysis", "qusage_dgsa_results_processed.rds")
p_summary_out     <- fs::path("output", "results", "reanalysis", "dgsa_comparison_summary.rds")

# ── Load and combine results ─────────────────────────────────────────────────

results_df <- bind_rows(
  readRDS(p_results_dearseq),
  readRDS(p_results_qusage)
)

# =============================================================================
# CONCORDANCE METRICS
# =============================================================================

concordance <- compute_concordance_metrics(
  results_df,
  methods       = c("dearseq", "qusage"),
  adj_pval_col  = "global.adjPval_BH",
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

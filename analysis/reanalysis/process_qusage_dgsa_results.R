# =============================================================================
# Preprocessing of QuSAGE DGSA Results
# =============================================================================
# Thin driver: tidies the raw QuSAGE results list into a long-format
# dataframe and applies p-value adjustment (R/postprocessing.R), then adds
# condition/gene-set-aggregate colours for plotting (R/plot_helpers.R). All
# the actual logic lives in R/.
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────

library(fs)
library(tidyverse)
library(RColorBrewer)

source(fs::path("R", "load_all.R"))

# ── Paths ─────────────────────────────────────────────────────────────────────

p_results_list <- fs::path("output", "results", "reanalysis", "qusage_dgsa_results_list.rds")
p_results_df   <- fs::path("output", "results", "reanalysis", "qusage_dgsa_results_processed.rds")
p_data_btm     <- fs::path("data", "BTM_processed.rds")

# ── Load data ─────────────────────────────────────────────────────────────────

qusage_dgsa_list <- readRDS(p_results_list)
BTM              <- readRDS(p_data_btm)

# =============================================================================
# TIDY, ADJUST, AND COLOUR
# =============================================================================

results_df <- build_tidy_dgsa_results(
  results_list     = qusage_dgsa_list,
  genesets         = BTM,
  conditions_order = default_conditions_order(),
  method           = "qusage"
) %>%
  assign_dgsa_colours(
    condition_colors = default_condition_colors(),
    aggregate_colors = default_aggregate_colors()
  ) %>%
  select(
    comparison, condition, time, condition.colour,
    gs.name, gs.description, gs.name.description, gs.aggregate, gs.colour,
    activation.score, fc.score, mean.corr, corr.mean,
    rawPval, method, everything()
  )

# ── Save ───────────────────────────────────────────────────────────────────────

saveRDS(results_df, file = p_results_df)
message("Processed results saved to: ", p_results_df)

rm(list = ls())

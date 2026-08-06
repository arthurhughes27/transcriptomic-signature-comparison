# =============================================================================
# Differential Gene-Set Analysis with QuSAGE — baseline specification
# =============================================================================
# Thin driver: runs the baseline QuSAGE specification (equal_variance =
# FALSE; sample_scope defaults to "study", reproducing Hagan et al.'s
# sample selection - see R/dgsa_qusage.R for why) across every vaccine x
# timepoint comparison, comparing post-vaccination samples against
# self-controls (pre-vaccination baseline). All the actual logic lives in
# R/dgsa_common.R and R/dgsa_qusage.R.
# Results are saved incrementally to guard against mid-run crashes.
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────

library(fs)
library(tidyverse)
library(qusage)

source(fs::path("R", "load_all.R"))

# ── Paths ─────────────────────────────────────────────────────────────────────

p_data_expr  <- fs::path("data", "hipc_merged_young_noNorm.rds")
p_data_btm   <- fs::path("data", "BTM_processed.rds")
p_results    <- fs::path("output", "results", "reanalysis", "qusage_dgsa_results_list.rds")

# ── Load data ─────────────────────────────────────────────────────────────────

hipc <- readRDS(p_data_expr)
BTM  <- readRDS(p_data_btm)

gene_names <- hipc %>% select(a1cf:zzz3) %>% colnames()

# =============================================================================
# MAIN LOOP
# =============================================================================

# Load any previously saved results so a restart continues rather than reruns
results_list <- if (file_exists(p_results)) readRDS(p_results) else list()

comparisons <- list_valid_comparisons(hipc)

for (i in seq_len(nrow(comparisons))) {

  vax <- comparisons$vaccine_name[i]
  day <- comparisons$day[i]

  comparison_name <- sprintf("%s vs Control - Day %s", vax, day)

  # Skip comparisons already completed in a previous run
  if (comparison_name %in% names(results_list)) {
    message("Skipping (already complete): ", comparison_name)
    next
  }

  message("Running: ", comparison_name)

  result <- tryCatch(
    run_qusage_comparison(vax, day, hipc, BTM, gene_names),
    error = function(e) {
      message("  ERROR: ", conditionMessage(e))
      NULL
    }
  )

  results_list[[comparison_name]] <- result

  # Save after every comparison so a crash loses at most one result
  saveRDS(object = results_list,
          file = p_results)

  # Free memory before the next iteration
  gc()
}

message("Analysis complete. Results saved to: ", p_results)

rm(list = ls())

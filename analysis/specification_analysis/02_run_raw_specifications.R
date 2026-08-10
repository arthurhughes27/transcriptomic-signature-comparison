# =============================================================================
# Specification analysis — 02: run raw specifications
# =============================================================================
# Runs each of the 66 raw specifications (output of 01_build_specification_
# grid.R) across every vaccine x DAYS_TO_ANALYSE timepoint comparison (days
# 1, 3, and 7 by default - edit DAYS_TO_ANALYSE below, or set it to NULL for
# every timepoint present in the data), using run_dearseq_comparison() /
# run_qusage_comparison(). Each specification's results are checkpointed to
# their own file (output/results/specification_analysis/raw/{spec_label}.rds),
# and comparisons already present in that file are skipped - so an
# interrupted run resumes rather than restarts, both across specifications
# and within one.
#
# THIS IS THE EXPENSIVE STEP: 66 specifications x every valid comparison,
# many involving dearseq's permutation test (1000 permutations for small
# samples). Expect this to take a long time on a laptop.
#
# Set SMOKE_TEST to TRUE below to restrict to a handful of comparisons
# first, to confirm the whole pipeline (this script -> 03 -> robustness
# metrics) runs end-to-end before committing to the full run. Re-run with
# SMOKE_TEST <- FALSE afterwards for the real analysis - this picks up
# wherever the smoke-test run left off (only the comparisons run under
# SMOKE_TEST will already be marked complete; delete
# output/results/specification_analysis/raw/ first if you'd rather start
# clean).
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────

library(fs)
library(tidyverse)
library(dearseq)
library(qusage)

source(fs::path("R", "load_all.R"))

# ── Timepoints to analyse ────────────────────────────────────────────────────
# Restrict to a subset of post-vaccination days for faster runs. Set to NULL
# to include every timepoint present in the data.

DAYS_TO_ANALYSE <- c(1, 3, 7)

# ── Smoke-test switch ────────────────────────────────────────────────────────

SMOKE_TEST <- FALSE
SMOKE_TEST_N_COMPARISONS <- 3

# ── Paths ─────────────────────────────────────────────────────────────────────

p_data_expr <- fs::path("data", "hipc_merged_young_noNorm.rds")
p_data_btm  <- fs::path("data", "BTM_processed.rds")
out_dir     <- fs::path("output", "results", "specification_analysis")
raw_dir     <- fs::path(out_dir, "raw")

fs::dir_create(raw_dir)

# ── Load data + specification grid ──────────────────────────────────────────

hipc <- readRDS(p_data_expr)
BTM  <- readRDS(p_data_btm)

gene_names <- hipc %>% select(a1cf:zzz3) %>% colnames()

# dearseq needs vaccine_code (pre/post indicator) and study_accession as a
# factor; harmless for qusage, which uses neither.
hipc <- hipc %>%
  mutate(
    vaccine_code    = as.factor(if_else(time_post_last_vax > 0, 2, 1)),
    study_accession = as.factor(study_accession)
  )

raw_grid <- readRDS(fs::path(out_dir, "raw_specification_grid.rds"))

comparisons <- list_valid_comparisons(hipc, days = DAYS_TO_ANALYSE)
if (SMOKE_TEST) {
  comparisons <- comparisons %>% slice_head(n = SMOKE_TEST_N_COMPARISONS)
  message("SMOKE_TEST = TRUE: restricting to ", nrow(comparisons), " comparison(s).")
}

comparison_names <- sprintf("%s vs Control - Day %s", comparisons$vaccine_name, comparisons$day)

# =============================================================================
# HELPER
# =============================================================================

# Dispatches to the right modular DGSA function for one raw specification's
# hyperparameters.
run_one_comparison <- function(spec, vax, day) {
  if (spec$method == "dearseq") {
    run_dearseq_comparison(
      vax, day, hipc, BTM, gene_names,
      covariates         = spec$covariates[[1]],
      which_weights       = spec$which_weights,
      gene_based_weights   = spec$gene_based_weights
    )
  } else {
    run_qusage_comparison(
      vax, day, hipc, BTM, gene_names,
      equal_variance = spec$equal_variance
    )
  }
}

# =============================================================================
# MAIN LOOP: EVERY RAW SPECIFICATION x EVERY COMPARISON
# =============================================================================

for (i in seq_len(nrow(raw_grid))) {

  spec        <- raw_grid[i, ]
  p_spec_file <- fs::path(raw_dir, paste0(spec$spec_label, ".rds"))

  results_list <- if (file_exists(p_spec_file)) readRDS(p_spec_file) else list()

  if (all(comparison_names %in% names(results_list))) {
    message(sprintf("[%d/%d] Skipping (already complete): %s", i, nrow(raw_grid), spec$spec_label))
    next
  }

  message(sprintf("[%d/%d] Running specification: %s", i, nrow(raw_grid), spec$spec_label))

  for (j in seq_len(nrow(comparisons))) {

    comparison_name <- comparison_names[j]

    if (comparison_name %in% names(results_list)) next

    message(sprintf("  [%d/%d] %s", j, nrow(comparisons), comparison_name))

    result <- tryCatch(
      run_one_comparison(spec, comparisons$vaccine_name[j], comparisons$day[j]),
      error = function(e) {
        message("    ERROR: ", conditionMessage(e))
        NULL
      }
    )

    results_list[[comparison_name]] <- result

    # Save after every comparison so a crash loses at most one result
    saveRDS(results_list, file = p_spec_file)
    gc()
  }
}

message("All raw specifications complete. Results saved to: ", raw_dir)

rm(list = ls())

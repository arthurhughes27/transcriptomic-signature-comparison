# =============================================================================
# Specification analysis — 03: apply post-hoc grid + accumulate robustness
# =============================================================================
# Loads each raw specification's results (output/results/specification_
# analysis/raw/, produced by 02_run_raw_specifications.R), tidies and
# p-value-adjusts them (R/postprocessing.R), and folds significance counts
# across the full post-hoc grid into the robustness metric pi_{g,v,j}
# (R/robustness_metrics.R) - without ever materialising the full
# 35,640-specification table (see the module-level comment at the top of
# R/robustness_metrics.R for why this is possible).
#
# Checkpointed at the level of "which raw specifications have been folded
# into the accumulator so far" (output/results/specification_analysis/
# robustness_accumulator_state.rds), so an interrupted run resumes rather
# than restarts. Safe to re-run at any point, including while
# 02_run_raw_specifications.R is still producing more raw specification
# results - specifications without a results file yet are skipped with a
# message and picked up on the next run of this script.
#
# Only the specification's own adjustment_method values (posthoc_methods,
# below) are folded into the robustness accumulator - build_tidy_dgsa_
# results() always computes adjusted p-values for all 6 of p.adjust()'s
# methods (for other consumers that want the full menu), not just the 3
# actually enumerated in the post-hoc grid, so counting every {scope}.
# adjPval_{method} column present would silently include 3 unintended
# extra methods (see the module-level comment in R/robustness_metrics.R).
#
# IMPORTANT: if you have an existing robustness_accumulator_state.rds from
# before this fix, its already-accumulated counts were computed with the
# unfiltered (18-column) method set and are NOT corrected by re-running
# this script as-is - resuming only processes NEW raw specifications, it
# doesn't re-derive ones already folded in. Delete
# output/results/specification_analysis/robustness_accumulator_state.rds
# (and robustness_metrics.rds) before re-running to get a fully correct
# result.
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────

library(fs)
library(tidyverse)

source(fs::path("R", "load_all.R"))

# ── Paths ─────────────────────────────────────────────────────────────────────

p_data_btm <- fs::path("data", "BTM_processed.rds")
out_dir    <- fs::path("output", "results", "specification_analysis")
raw_dir    <- fs::path(out_dir, "raw")
p_state    <- fs::path(out_dir, "robustness_accumulator_state.rds")
p_metrics  <- fs::path(out_dir, "robustness_metrics.rds")

# ── Load data + specification grids ─────────────────────────────────────────

BTM          <- readRDS(p_data_btm)
raw_grid     <- readRDS(fs::path(out_dir, "raw_specification_grid.rds"))
posthoc_grid <- readRDS(fs::path(out_dir, "posthoc_specification_grid.rds"))

alphas          <- sort(unique(posthoc_grid$alpha))
fc_thresholds   <- sort(unique(posthoc_grid$fc_threshold))
posthoc_methods <- unique(posthoc_grid$adjustment_method)

# ── Load or initialise checkpoint state ─────────────────────────────────────

state <- if (file_exists(p_state)) {
  readRDS(p_state)
} else {
  list(accumulator = NULL, processed_spec_ids = integer(0))
}

# =============================================================================
# FOLD EACH RAW SPECIFICATION'S RESULTS INTO THE ACCUMULATOR
# =============================================================================

for (i in seq_len(nrow(raw_grid))) {

  spec <- raw_grid[i, ]

  if (spec$raw_spec_id %in% state$processed_spec_ids) {
    message(sprintf("[%d/%d] Skipping (already accumulated): %s", i, nrow(raw_grid), spec$spec_label))
    next
  }

  p_spec_file <- fs::path(raw_dir, paste0(spec$spec_label, ".rds"))

  if (!file_exists(p_spec_file)) {
    message(sprintf(
      "[%d/%d] No results file yet for: %s - skipping (run 02_run_raw_specifications.R first)",
      i, nrow(raw_grid), spec$spec_label
    ))
    next
  }

  message(sprintf("[%d/%d] Accumulating: %s", i, nrow(raw_grid), spec$spec_label))

  raw_results_list <- readRDS(p_spec_file)

  tidy_df <- tryCatch(
    build_tidy_dgsa_results(raw_results_list, BTM, method = spec$method),
    error = function(e) {
      message("    ERROR tidying results (skipping this specification): ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(tidy_df)) next

  state$accumulator <- accumulate_robustness_counts(
    state$accumulator, tidy_df,
    alphas = alphas, fc_thresholds = fc_thresholds, methods = posthoc_methods
  )
  state$processed_spec_ids <- c(state$processed_spec_ids, spec$raw_spec_id)

  saveRDS(state, file = p_state)
}

# ── Compute and save (possibly partial) robustness metrics ─────────────────

n_processed <- length(state$processed_spec_ids)

if (n_processed == 0) {
  stop("No raw specifications have been accumulated yet - run 02_run_raw_specifications.R first.")
}

if (n_processed < nrow(raw_grid)) {
  warning(sprintf(
    "Only %d/%d raw specifications have been accumulated so far - robustness metrics below are PARTIAL. Re-run this script after 02_run_raw_specifications.R completes for the full result.",
    n_processed, nrow(raw_grid)
  ))
}

robustness_metrics <- compute_robustness_metric(state$accumulator)

saveRDS(robustness_metrics, file = p_metrics)
message(sprintf(
  "Robustness metrics (%d/%d specifications) saved to: %s",
  n_processed, nrow(raw_grid), p_metrics
))

rm(list = ls())

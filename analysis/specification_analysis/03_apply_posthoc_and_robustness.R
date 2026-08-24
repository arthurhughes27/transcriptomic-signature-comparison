# =============================================================================
# Specification analysis — 03: apply post-hoc grid + accumulate robustness
# =============================================================================
# Loads each raw specification's results (output/results/specification_
# analysis/raw/, produced by 02_run_raw_specifications.R), tidies and
# p-value-adjusts them (R/postprocessing.R), and folds significance counts
# across the full post-hoc grid into the robustness metric pi_{g,v,j}
# (R/robustness_metrics.R) - without ever materialising the full
# 2,754-specification table (see the module-level comment at the top of
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
# IMPORTANT: whenever raw_specification_grid.rds changes - e.g. after the
# fix above, or after removing gene_based_weights as an investigated axis
# in R/specifications.R's build_raw_specification_grid() (dropping the raw
# grid from 66 to 34 rows) - any existing
# robustness_accumulator_state.rds/robustness_metrics.rds was accumulated
# against the OLD raw grid and will NOT be corrected by simply re-running
# this script: resuming only processes raw specifications not already
# folded in, it never re-derives ones that are, and it also never removes
# a since-dropped specification's contribution. Delete
# output/results/specification_analysis/robustness_accumulator_state.rds
# (and robustness_metrics.rds) before re-running whenever the raw grid
# itself has changed, to get a fully correct result. Re-running
# 02_run_raw_specifications.R is NOT required in this case - it reuses
# whichever already-computed output/results/specification_analysis/raw/
# {spec_label}.rds files match the (possibly smaller) new raw grid.
#
# COVARIATE-REDUNDANCY DEDUPLICATION (R/specification_redundancy.R): some
# comparisons (vaccine x timepoint) have no variation in one or more
# covariates within their paired sample set - e.g. a single contributing
# study, or an all-"Unknown" gender/race - so dearseq's covariate design
# matrix collapses and two of the 16 covariate subsets produce
# byte-identical results for that comparison, even though
# build_raw_specification_grid() counts them as separate raw
# specifications. This is detected purely from the already-loaded sample
# data (build_covariate_matrix() calls only - no re-run of dearseq/qusage
# model fitting) and used below to drop the redundant duplicate's
# contribution per comparison, so it isn't double-counted in n_evaluated/
# n_significant. Like a raw-grid change, this is a change to what counts
# as "already accumulated": an existing robustness_accumulator_state.rds/
# robustness_metrics.rds predates this fix and must be deleted before
# re-running this script to get corrected results (re-running
# 02_run_raw_specifications.R is NOT required - only already-cached raw
# results are reprocessed).
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────

library(fs)
library(tidyverse)

source(fs::path("R", "load_all.R"))

# ── Paths ─────────────────────────────────────────────────────────────────────

p_data_btm  <- fs::path("data", "BTM_processed.rds")
p_data_expr <- fs::path("data", "hipc_merged_young_noNorm.rds")
out_dir     <- fs::path("output", "results", "specification_analysis")
raw_dir     <- fs::path(out_dir, "raw")
p_state     <- fs::path(out_dir, "robustness_accumulator_state.rds")
p_metrics   <- fs::path(out_dir, "robustness_metrics.rds")
p_spec_totals <- fs::path(out_dir, "effective_specification_totals.csv")

# Must match 02_run_raw_specifications.R's DAYS_TO_ANALYSE - the redundancy
# detection below (R/specification_redundancy.R) needs the exact same set
# of comparisons that were actually run.
DAYS_TO_ANALYSE <- c(1, 3, 7)

# ── Load data + specification grids ─────────────────────────────────────────

BTM          <- readRDS(p_data_btm)
raw_grid     <- readRDS(fs::path(out_dir, "raw_specification_grid.rds"))
posthoc_grid <- readRDS(fs::path(out_dir, "posthoc_specification_grid.rds"))

alphas          <- sort(unique(posthoc_grid$alpha))
fc_thresholds   <- sort(unique(posthoc_grid$fc_threshold))
posthoc_methods <- unique(posthoc_grid$adjustment_method)

# ── Covariate-redundancy detection (no DGSA re-run - see the module-level
# comment above and R/specification_redundancy.R) ───────────────────────────

# Same vaccine_code/study_accession preprocessing 02_run_raw_specifications.R
# applies before calling run_dearseq_comparison() - build_covariate_matrix()
# needs study_accession as a factor to reproduce the real design matrices.
hipc <- readRDS(p_data_expr) |>
  dplyr::mutate(
    vaccine_code    = as.factor(if_else(time_post_last_vax > 0, 2, 1)),
    study_accession = as.factor(study_accession)
  )

effective_specs <- effective_raw_specifications(hipc, raw_grid, days = DAYS_TO_ANALYSE)

# One row per (raw_spec_id, comparison) pair whose contribution should be
# dropped for that comparison specifically - the canonical member of its
# covariate-equivalence class is kept, so the pair isn't double-counted.
redundant_lookup <- effective_specs |>
  dplyr::filter(!is_canonical) |>
  dplyr::select(raw_spec_id, vaccine_name, day)

message(sprintf(
  "Covariate-redundancy detection: %d (raw specification x comparison) pair(s) flagged as redundant duplicates.",
  nrow(redundant_lookup)
))

spec_totals <- specification_totals(hipc, raw_grid, posthoc_grid, days = DAYS_TO_ANALYSE)
write.csv(spec_totals, file = p_spec_totals, row.names = FALSE)
message("Deduplicated per-comparison specification totals saved to: ", p_spec_totals)

rm(hipc)

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

  # Coerce unconditionally (not just when redundant_here has rows below) so
  # `condition`/`time`'s types are consistent across every contribution
  # folded into the accumulator, regardless of which raw specifications
  # happen to have a redundant duplicate. build_tidy_dgsa_results()
  # (R/postprocessing.R) returns `condition` and `time` as factors;
  # order_robustness_comparisons() (R/robustness_heatmaps.R) re-derives
  # both from character/numeric downstream anyway, so this is safe.
  tidy_df <- dplyr::mutate(
    tidy_df,
    condition = as.character(condition),
    time      = as.numeric(as.character(time))
  )

  # Drop this raw specification's rows for any comparison where it's a
  # covariate-redundant duplicate (see the module-level comment above) -
  # the canonical member of its equivalence class carries that
  # comparison's contribution instead, so it's counted exactly once.
  redundant_here <- dplyr::filter(redundant_lookup, raw_spec_id == spec$raw_spec_id)
  if (nrow(redundant_here) > 0) {
    n_before <- nrow(tidy_df)
    tidy_df <- dplyr::anti_join(
      tidy_df, redundant_here,
      by = c("condition" = "vaccine_name", "time" = "day")
    )
    message(sprintf(
      "    Dropped %d row(s) for %d covariate-redundant comparison(s).",
      n_before - nrow(tidy_df), nrow(redundant_here)
    ))
  }

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

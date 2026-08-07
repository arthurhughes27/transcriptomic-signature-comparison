# =============================================================================
# Baseline-specification results (single source of truth)
# =============================================================================
# The "baseline" DGSA result for each method (the bolded options in Table
# 2.1) is one row of the raw specification grid (R/specifications.R,
# is_baseline == TRUE), already run and saved by
# analysis/specification_analysis/02_run_raw_specifications.R as part of
# the specification-analysis loop. Everything in the repository that needs
# "the baseline result" (the Stage 1 dearseq-vs-QuSAGE comparison, Figure 3,
# the robustness tables) should read it from there, rather than from a
# separately-run baseline-only script (analysis/reanalysis/dearseq_dgsa.R /
# qusage_dgsa.R) - the latter is a distinct run of the same specification
# and is not guaranteed to reproduce identically (observed in practice:
# gene sets with baseline-significant p-values but zero robustness, which
# is impossible if the baseline run were among the specifications folded
# into the robustness accumulator).
# =============================================================================

#' Load and tidy the baseline-specification results from the raw
#' specification-grid outputs
#'
#' @param raw_grid Raw specification grid (output of
#'   [build_raw_specification_grid()], or the saved
#'   `raw_specification_grid.rds`), with `is_baseline`, `spec_label`, and
#'   `method` columns.
#' @param raw_dir Directory containing each raw specification's results
#'   file (`output/results/specification_analysis/raw/`, as written by
#'   `02_run_raw_specifications.R`).
#' @param genesets Gene set list (as in `data/BTM_processed.rds` /
#'   `data/BG3M_processed.rds`), passed to [build_tidy_dgsa_results()].
#' @param conditions_order Optional vaccine-condition factor ordering,
#'   passed through to [build_tidy_dgsa_results()] for every baseline
#'   specification loaded, so `condition` has the same factor levels
#'   (and order) across methods when the results are row-bound - matters
#'   for anything that relies on `levels(result$condition)` downstream
#'   (e.g. `analysis/reanalysis/plot_circos.R`). NULL (the default) uses
#'   each specification's own conditions, sorted alphabetically.
#'
#' @return A tidy, p-value-adjusted results tibble (`dplyr::bind_rows()`
#'   across every `is_baseline == TRUE` raw specification - typically one
#'   per method), in the same shape as [build_tidy_dgsa_results()]'s output
#'   for a single specification.
load_baseline_results_from_raw <- function(raw_grid, raw_dir, genesets, conditions_order = NULL) {
  baseline_specs <- dplyr::filter(raw_grid, is_baseline)

  if (nrow(baseline_specs) == 0) {
    stop("No is_baseline == TRUE rows in raw_grid.")
  }

  purrr::map_dfr(seq_len(nrow(baseline_specs)), function(i) {
    spec        <- baseline_specs[i, ]
    p_spec_file <- fs::path(raw_dir, paste0(spec$spec_label, ".rds"))

    if (!fs::file_exists(p_spec_file)) {
      stop(
        "Baseline raw specification results not found at ", p_spec_file, ".\n",
        "Run 02_run_raw_specifications.R first (it must include the baseline specifications)."
      )
    }

    raw_results_list <- readRDS(p_spec_file)
    build_tidy_dgsa_results(
      raw_results_list, genesets,
      conditions_order = conditions_order, method = spec$method
    )
  })
}

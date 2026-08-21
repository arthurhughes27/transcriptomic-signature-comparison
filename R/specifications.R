# =============================================================================
# Specification grid generation (Chapter 2, Table 2.1)
# =============================================================================
# The 2,754 specifications considered factor into two independent tiers:
#
#   - "raw" specifications (34 total: 32 dearseq + 2 QuSAGE) require actually
#     re-running a DGSA method (method; covariate set; weight method for
#     dearseq; equal-variance assumption for QuSAGE).
#   - "post-hoc" specifications (81) are cheap, vectorised operations
#     applied to an already-computed raw run's p-values and fold-changes
#     (p-value adjustment scope; adjustment method; significance level;
#     fold-change threshold).
#
# 34 raw runs x 81 post-hoc combinations = 2,754. Running DGSA 34 times (not
# 2,754 times) and applying the post-hoc grid afterwards
# (R/postprocessing.R, R/robustness_metrics.R) is what makes the full
# specification analysis computationally tractable.
#
# The weight-estimation LEVEL hyperparameter (gene-level vs
# observation-level, dearseq only) was dropped from the raw grid entirely -
# only the observation-level baseline is run - rather than investigated as
# a specification axis; `gene_based_weights` is kept as an always-FALSE
# column (not removed) so spec_label naming and downstream schema stay
# stable.
# =============================================================================

#' Build the raw specification grid (34 rows)
#'
#' One row per specification that requires actually running a DGSA method:
#' 32 dearseq specifications (16 covariate subsets x 2 weight methods,
#' observation-level weight estimation only) and 2 QuSAGE specifications
#' (equal-variance assumption).
#'
#' @return A tibble with columns `raw_spec_id`, `method` ("dearseq" or
#'   "qusage"), `covariates` (list-column of character vectors; NULL for
#'   QuSAGE), `which_weights`, `gene_based_weights` (dearseq only, always
#'   FALSE - observation-level weight estimation is no longer investigated
#'   as a specification axis; NA for QuSAGE), `equal_variance` (QuSAGE
#'   only; NA for dearseq), `is_baseline` (TRUE for the bolded options in
#'   Table 2.1), and `spec_label` (a human-readable identifier suitable for
#'   file naming).
build_raw_specification_grid <- function() {
  covariate_pool <- c("age_imputed", "gender", "study_accession", "race")
  covariate_sets <- unlist(
    lapply(0:length(covariate_pool), function(k) utils::combn(covariate_pool, k, simplify = FALSE)),
    recursive = FALSE
  )

  dearseq_grid <- tidyr::expand_grid(
    covariates          = covariate_sets,
    which_weights        = c("loclin", "voom"),
    gene_based_weights    = FALSE
  ) |>
    dplyr::mutate(
      method          = "dearseq",
      equal_variance   = NA,
      is_baseline       = purrr::map_lgl(covariates, setequal, covariate_pool) &
                          which_weights == "loclin",
      spec_label         = purrr::pmap_chr(
        list(covariates, which_weights, gene_based_weights),
        function(cov, wt, lvl) {
          cov_label <- if (length(cov) == 0) "none" else paste(cov, collapse = "-")
          sprintf("dearseq_cov-%s_wt-%s_lvl-%s", cov_label, wt, if (lvl) "gene" else "obs")
        }
      )
    )

  qusage_grid <- tibble::tibble(equal_variance = c(FALSE, TRUE)) |>
    dplyr::mutate(
      method               = "qusage",
      covariates            = list(NULL),
      which_weights          = NA_character_,
      gene_based_weights      = NA,
      is_baseline              = !equal_variance,
      spec_label                = sprintf("qusage_eqvar-%s", equal_variance)
    )

  dplyr::bind_rows(dearseq_grid, qusage_grid) |>
    dplyr::mutate(raw_spec_id = dplyr::row_number()) |>
    dplyr::select(raw_spec_id, method, covariates, which_weights, gene_based_weights,
                  equal_variance, is_baseline, spec_label)
}

#' Build the post-hoc specification grid (81 rows)
#'
#' One row per specification applied after a raw DGSA run: p-value
#' adjustment scope, adjustment method, significance level, and mean
#' absolute log2-fold-change threshold.
#'
#' @return A tibble with columns `posthoc_spec_id`, `adjustment_scope`
#'   ("global", "withinTime", or "withinComparison" - matching the column
#'   naming produced by [apply_pvalue_adjustments()]), `adjustment_method`,
#'   `alpha`, `fc_threshold`, and `is_baseline` (TRUE for the bolded options
#'   in Table 2.1: within timepoint, BH, 0.05, 0.0).
build_posthoc_specification_grid <- function() {
  tidyr::expand_grid(
    adjustment_scope   = c("global", "withinTime", "withinComparison"),
    # adjustment_method    = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY"),
    adjustment_method    = c("holm","BH", "BY"),
    # alpha                  = c(0.0001, 0.001, 0.01, 0.05, 0.1),
    alpha                    =  c(0.01, 0.05, 0.1),
    # fc_threshold             = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
    fc_threshold             = c(0, 0.5, 1)
  ) |>
    dplyr::mutate(
      is_baseline = adjustment_scope == "withinTime" &
                    adjustment_method == "BH" &
                    alpha == 0.05 &
                    fc_threshold == 0,
      posthoc_spec_id = dplyr::row_number()
    ) |>
    dplyr::select(posthoc_spec_id, adjustment_scope, adjustment_method, alpha, fc_threshold, is_baseline)
}

#' Build the full specification grid (raw x post-hoc cross join)
#'
#' Materialises every one of the 2,754 specifications considered as a
#' single row. Intended for validation/testing and small-scale inspection
#' only - the main specification-analysis pipeline (R/robustness_metrics.R)
#' deliberately avoids materialising this at full scale (2,754
#' specifications x every gene set x comparison would be very large), instead
#' running each raw specification once and accumulating post-hoc results.
#'
#' @param raw_grid Output of [build_raw_specification_grid()].
#' @param posthoc_grid Output of [build_posthoc_specification_grid()].
#'
#' @return A tibble with `nrow(raw_grid) * nrow(posthoc_grid)` rows, `spec_id`,
#'   every column from `raw_grid` and `posthoc_grid`, and `is_baseline`
#'   (TRUE only for the single specification matching both grids' baselines).
build_full_specification_grid <- function(raw_grid = build_raw_specification_grid(),
                                          posthoc_grid = build_posthoc_specification_grid()) {
  tidyr::expand_grid(
    raw_spec_id     = raw_grid$raw_spec_id,
    posthoc_spec_id = posthoc_grid$posthoc_spec_id
  ) |>
    dplyr::left_join(raw_grid, by = "raw_spec_id") |>
    dplyr::left_join(posthoc_grid, by = "posthoc_spec_id", suffix = c("_raw", "_posthoc")) |>
    dplyr::mutate(
      is_baseline = is_baseline_raw & is_baseline_posthoc,
      spec_id      = dplyr::row_number()
    ) |>
    dplyr::select(-is_baseline_raw, -is_baseline_posthoc) |>
    dplyr::relocate(spec_id)
}

#' Get the baseline raw specifications
#'
#' Convenience accessor for the two specifications (one dearseq, one
#' QuSAGE) matching the bolded options in Table 2.1 - used by the Stage 1
#' (baseline method comparison) driver.
#'
#' @param raw_grid Output of [build_raw_specification_grid()].
#'
#' @return A tibble with one row per method, filtered to `is_baseline`.
baseline_raw_specifications <- function(raw_grid = build_raw_specification_grid()) {
  dplyr::filter(raw_grid, is_baseline)
}

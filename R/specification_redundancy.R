# =============================================================================
# Raw-specification redundancy detection (per comparison)
# =============================================================================
# For a given vaccine x timepoint comparison, dearseq's covariate design
# matrix (build_covariate_matrix(), R/dgsa_common.R) collapses whenever a
# covariate doesn't vary within that comparison's paired sample set - e.g. a
# single contributing study, an all-"Unknown" gender/race, or (rarer) a
# constant age_imputed. R's model.matrix() drops a single-level factor
# entirely (0 columns), and build_covariate_matrix()'s QR-based
# dependent-column removal drops any resulting rank-deficient column too
# (e.g. a constant numeric covariate, or any subtler collinearity between
# covariates) - so two of the 16 dearseq covariate subsets that reduce to
# the same design-matrix columns for a given comparison produce
# byte-identical DGSA results for that comparison, even though
# build_raw_specification_grid() counts them as separate raw
# specifications. This degeneracy is comparison-specific (it depends on
# which participants/studies happen to fall into that particular vaccine x
# timepoint's paired sample set), not global.
#
# Detected here purely from already-loaded sample data - calling
# build_covariate_matrix() directly (cheap linear algebra) rather than
# dearseq::dgsa_seq()/qusage model fitting - so it never requires
# re-running DGSA. Used by 03_apply_posthoc_and_robustness.R to avoid
# double-counting these redundant pairs when accumulating the robustness
# metric, and to report the true (deduplicated) number of specifications
# considered per comparison.
# =============================================================================

#' Group a comparison's dearseq covariate subsets by design-matrix
#' equivalence
#'
#' Builds the dearseq covariate design matrix (build_covariate_matrix()) for
#' every element of `covariate_sets`, against one comparison's paired sample
#' set, and groups subsets that produce the same resulting columns (implying
#' identical DGSA results for THIS comparison specifically).
#'
#' @param df One comparison's paired sample set (output of
#'   [filter_paired_samples()], with the same `vaccine_code`/`study_accession`
#'   preprocessing 02_run_raw_specifications.R applies before calling
#'   [run_dearseq_comparison()]).
#' @param covariate_sets List of character vectors, one per covariate
#'   subset to check (as in [build_raw_specification_grid()]'s `covariates`
#'   column) - may contain duplicates (e.g. once per `which_weights`); each
#'   element is checked independently.
#'
#' @return An integer vector aligned to `covariate_sets`: the equivalence
#'   class id of each subset (subsets sharing an id produce identical
#'   design-matrix columns, and therefore identical DGSA results, for this
#'   comparison).
#' @keywords internal
covariate_subset_equivalence <- function(df, covariate_sets) {
  signatures <- vapply(covariate_sets, function(covs) {
    x <- suppressMessages(build_covariate_matrix(df, covariates = covs))
    paste(sort(colnames(x)), collapse = "|")
  }, character(1))

  match(signatures, unique(signatures))
}

#' Flag, per dearseq raw specification x comparison, whether it's the
#' canonical (counted) member of its comparison-specific equivalence class
#'
#' For each comparison, groups the 16 dearseq covariate subsets (per
#' `which_weights`, which doesn't affect the covariate matrix and so is
#' handled independently) into equivalence classes via
#' [covariate_subset_equivalence()], then picks one canonical member per
#' class - preferring the nominal baseline covariate set
#' (`raw_grid$is_baseline`) when it belongs to the class, else the member
#' with the most covariates, else alphabetically by `spec_label` - purely
#' for interpretability, since every member of a class is guaranteed to
#' produce identical results for that comparison.
#'
#' @param hipc Merged clinical/expression tibble, with the same
#'   `vaccine_code`/`study_accession` preprocessing
#'   02_run_raw_specifications.R applies (as.factor(study_accession), and a
#'   `vaccine_code` pre/post indicator).
#' @param raw_grid Output of [build_raw_specification_grid()].
#' @param days Optional vector restricting to specific timepoints (passed
#'   to [list_valid_comparisons()]) - pass the same `DAYS_TO_ANALYSE` used
#'   to actually run the raw specifications.
#'
#' @return A tibble with one row per (dearseq raw_spec_id x comparison)
#'   pair with data (comparisons with zero paired samples are skipped, same
#'   as [run_dearseq_comparison()] would skip them), columns `raw_spec_id`,
#'   `vaccine_name`, `day`, `equivalence_class`, `is_canonical`. QuSAGE raw
#'   specs are not included (no covariate axis, so no redundancy applies to
#'   them) - treat every QuSAGE raw_spec_id as canonical for every
#'   comparison when using this table downstream.
effective_raw_specifications <- function(hipc, raw_grid, days = NULL) {
  dearseq_grid <- dplyr::filter(raw_grid, method == "dearseq")
  comparisons  <- list_valid_comparisons(hipc, days = days)

  purrr::pmap_dfr(comparisons, function(vaccine_name, day) {
    df <- filter_paired_samples(hipc, vaccine_name, day)
    if (nrow(df) == 0) return(NULL)

    dearseq_grid |>
      dplyr::mutate(
        vaccine_name       = vaccine_name,
        day                = day,
        equivalence_class  = covariate_subset_equivalence(df, covariates),
        n_covariates       = lengths(covariates)
      ) |>
      dplyr::group_by(which_weights, equivalence_class) |>
      dplyr::arrange(
        dplyr::desc(is_baseline), dplyr::desc(n_covariates), spec_label,
        .by_group = TRUE
      ) |>
      dplyr::mutate(is_canonical = dplyr::row_number() == 1) |>
      dplyr::ungroup() |>
      dplyr::select(raw_spec_id, vaccine_name, day, equivalence_class, is_canonical)
  })
}

#' Deduplicated total-specification counts per comparison
#'
#' @param hipc Merged clinical/expression tibble (see
#'   [effective_raw_specifications()] for required preprocessing).
#' @param raw_grid Output of [build_raw_specification_grid()].
#' @param posthoc_grid Output of [build_posthoc_specification_grid()].
#' @param days Optional vector restricting to specific timepoints.
#'
#' @return A tibble with one row per comparison: `vaccine_name`, `day`,
#'   `n_effective_raw_specs` (deduplicated raw-specification count -
#'   canonical dearseq members only, plus every QuSAGE spec, since QuSAGE
#'   has no covariate-redundancy axis) and `n_effective_total_specs`
#'   (`n_effective_raw_specs * nrow(posthoc_grid)`) - the number of
#'   DISTINCT specifications actually considered for that comparison. A raw
#'   specification producing no result at all for a comparison (e.g. a
#'   method-specific sample-size filter) can still further reduce a given
#'   gene set x comparison's `n_evaluated` below this figure - see the
#'   module-level comment in R/robustness_metrics.R.
specification_totals <- function(hipc, raw_grid, posthoc_grid, days = NULL) {
  effective <- effective_raw_specifications(hipc, raw_grid, days = days)
  n_qusage  <- sum(raw_grid$method == "qusage")

  effective |>
    dplyr::group_by(vaccine_name, day) |>
    dplyr::summarise(n_effective_raw_specs = sum(is_canonical) + n_qusage, .groups = "drop") |>
    dplyr::mutate(n_effective_total_specs = n_effective_raw_specs * nrow(posthoc_grid))
}

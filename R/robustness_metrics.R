# =============================================================================
# Robustness metric computation
# =============================================================================
# Computes the robustness metric pi_{g,v,j} from Chapter 2, Section 2.2.4:
# the proportion of specifications under which a given gene set x comparison
# is called significant.
#
# For a fixed adjusted-p-value column, whether a row passes significance
# level alpha and fold-change threshold tau factors as two independent
# conditions, so the count of (alpha, tau) pairs passed by a row is simply
# (number of alphas passed) x (number of tau passed) - this lets us count
# how many of the 540 post-hoc specifications (18 scope x method columns x
# 5 alphas x 6 thresholds) call each row significant without ever
# constructing the full row x specification matrix. Accumulating these
# counts across the 66 raw runs (R/specifications.R,
# R/dgsa_dearseq.R/dgsa_qusage.R) gives pi_{g,v,j} without materialising the
# full 35,640-specification table.
# =============================================================================

#' Count, for each row, how many post-hoc specifications call it significant
#'
#' For one raw DGSA run's already-adjusted results (see
#' [build_tidy_dgsa_results()]), counts, for every row (gene set x
#' comparison), how many post-hoc specifications - every combination of the
#' adjusted p-value columns present in `tidy_df`, `alphas`, and
#' `fc_thresholds` - call that row significant. A specification calls a row
#' significant when its adjusted p-value is <= alpha AND its absolute
#' fold-change score is >= the fold-change threshold.
#'
#' @param tidy_df Output of [build_tidy_dgsa_results()] (or any tibble with
#'   one or more `{scope}.adjPval_{method}` columns and an `fc.score`
#'   column).
#' @param alphas Numeric vector of significance levels.
#' @param fc_thresholds Numeric vector of mean absolute log2-fold-change
#'   thresholds.
#'
#' @return An integer vector, one value per row of `tidy_df`, in
#'   `[0, n_adjpval_cols * length(alphas) * length(fc_thresholds)]`.
count_significant_specifications <- function(tidy_df,
                                             alphas        = c(0.0001, 0.001, 0.01, 0.05, 0.1),
                                             fc_thresholds  = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) {

  stopifnot("fc.score" %in% colnames(tidy_df))

  adjpval_cols <- grep("\\.adjPval_", colnames(tidy_df), value = TRUE)
  if (length(adjpval_cols) == 0) {
    stop("tidy_df has no {scope}.adjPval_{method} columns; run apply_pvalue_adjustments() first.")
  }

  passes <- function(values, thresholds, comparator) {
    rowSums(outer(values, thresholds, function(v, t) !is.na(v) & comparator(v, t)))
  }

  n_fc_pass <- passes(abs(tidy_df$fc.score), fc_thresholds, `>=`)

  counts <- rep(0L, nrow(tidy_df))
  for (col in adjpval_cols) {
    n_alpha_pass <- passes(tidy_df[[col]], alphas, `<=`)
    counts <- counts + as.integer(n_alpha_pass * n_fc_pass)
  }

  counts
}

#' Accumulate per-specification significance counts across raw DGSA runs
#'
#' Maintains a running total, keyed by (gene set, vaccine condition,
#' timepoint), of how many post-hoc specifications call each gene set x
#' comparison significant, and how many post-hoc specifications were
#' actually evaluated for it. The latter may differ across gene set x
#' comparisons when a raw specification does not produce a result for some
#' comparison (e.g. a method-specific sample-size filter), so `n_evaluated`
#' is tracked explicitly rather than assumed to equal the full grid size.
#'
#' Call this once per raw specification's tidy, adjusted results (output of
#' [build_tidy_dgsa_results()]), threading the returned accumulator into the
#' next call; pass `accumulator = NULL` to start.
#'
#' @param accumulator NULL to start a new accumulator, or the tibble
#'   returned by a previous call to this function.
#' @param tidy_df Output of [build_tidy_dgsa_results()] for one raw
#'   specification.
#' @param alphas,fc_thresholds Passed to
#'   [count_significant_specifications()].
#'
#' @return A tibble with columns `gs.name`, `condition`, `time`,
#'   `n_significant` (running total of significant-specification counts),
#'   and `n_evaluated` (running total of specifications evaluated).
accumulate_robustness_counts <- function(accumulator, tidy_df,
                                         alphas        = c(0.0001, 0.001, 0.01, 0.05, 0.1),
                                         fc_thresholds  = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) {

  n_adjpval_cols  <- length(grep("\\.adjPval_", colnames(tidy_df)))
  n_posthoc_specs <- n_adjpval_cols * length(alphas) * length(fc_thresholds)

  contribution <- tibble::tibble(
    gs.name       = tidy_df$gs.name,
    condition     = tidy_df$condition,
    time          = tidy_df$time,
    n_significant = count_significant_specifications(tidy_df, alphas, fc_thresholds),
    n_evaluated   = n_posthoc_specs
  )

  if (is.null(accumulator)) return(contribution)

  dplyr::bind_rows(accumulator, contribution) |>
    dplyr::group_by(gs.name, condition, time) |>
    dplyr::summarise(
      n_significant = sum(n_significant),
      n_evaluated   = sum(n_evaluated),
      .groups = "drop"
    )
}

#' Compute the robustness metric from accumulated significance counts
#'
#' @param accumulator Output of (repeated calls to)
#'   [accumulate_robustness_counts()].
#'
#' @return `accumulator` with an added `robustness` column: the proportion
#'   of evaluated specifications under which each gene set x comparison was
#'   called significant (pi_{g,v,j} in Chapter 2, Section 2.2.4).
compute_robustness_metric <- function(accumulator) {
  dplyr::mutate(accumulator, robustness = n_significant / n_evaluated)
}

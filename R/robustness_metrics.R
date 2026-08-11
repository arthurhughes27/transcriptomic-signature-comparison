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
# how many post-hoc specifications (scope x method columns x alphas x
# thresholds) call each row significant without ever constructing the full
# row x specification matrix. Accumulating these counts across the 66 raw
# runs (R/specifications.R, R/dgsa_dearseq.R/dgsa_qusage.R) gives
# pi_{g,v,j} without materialising the full specification table.
#
# build_tidy_dgsa_results() (R/postprocessing.R) always computes adjusted
# p-values for all 6 of p.adjust()'s methods (for other consumers, e.g.
# the circos plots and heatmap comparison, that want the full menu), not
# just the 3 methods actually enumerated in the post-hoc specification
# grid (R/specifications.R's build_posthoc_specification_grid() - holm,
# BH, BY). Both functions below therefore take an explicit `methods`
# argument restricting which `{scope}.adjPval_{method}` columns count as
# "post-hoc specifications" - 03_apply_posthoc_and_robustness.R passes
# the post-hoc grid's own adjustment_method values. Without this filter,
# every one of the 3 unintended extra methods (hochberg, hommel,
# bonferroni) present in tidy_df would silently be folded into the
# robustness computation too, inflating n_evaluated (and, non-uniformly,
# n_significant).
# =============================================================================

# Restrict tidy_df's {scope}.adjPval_{method} columns to just `methods`
# (NULL = every such column present, the old, unfiltered behaviour).
select_adjpval_columns <- function(tidy_df, methods = NULL) {
  adjpval_cols <- grep("\\.adjPval_", colnames(tidy_df), value = TRUE)
  if (!is.null(methods)) {
    method_pattern <- paste0("\\.adjPval_(", paste(methods, collapse = "|"), ")$")
    adjpval_cols   <- grep(method_pattern, adjpval_cols, value = TRUE)
  }
  adjpval_cols
}

#' Count, for each row, how many post-hoc specifications call it significant
#'
#' For one raw DGSA run's already-adjusted results (see
#' [build_tidy_dgsa_results()]), counts, for every row (gene set x
#' comparison), how many post-hoc specifications - every combination of the
#' selected adjusted p-value columns, `alphas`, and `fc_thresholds` - call
#' that row significant. A specification calls a row significant when its
#' adjusted p-value is <= alpha AND its absolute fold-change score is >=
#' the fold-change threshold.
#'
#' @param tidy_df Output of [build_tidy_dgsa_results()] (or any tibble with
#'   one or more `{scope}.adjPval_{method}` columns and an `fc.score`
#'   column).
#' @param alphas Numeric vector of significance levels.
#' @param fc_thresholds Numeric vector of mean absolute log2-fold-change
#'   thresholds.
#' @param methods Optional character vector restricting which
#'   `{scope}.adjPval_{method}` columns to count (matched against the
#'   `_{method}` suffix) - pass the post-hoc specification grid's own
#'   `adjustment_method` values so only the intended post-hoc grid is
#'   counted, not every adjustment method `tidy_df` happens to carry. NULL
#'   (the default) counts every `{scope}.adjPval_{method}` column present.
#'
#' @return An integer vector, one value per row of `tidy_df`, in
#'   `[0, n_adjpval_cols * length(alphas) * length(fc_thresholds)]`.
count_significant_specifications <- function(tidy_df,
                                             alphas        = c(0.0001, 0.001, 0.01, 0.05, 0.1),
                                             fc_thresholds  = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
                                             methods         = NULL) {

  stopifnot("fc.score" %in% colnames(tidy_df))

  adjpval_cols <- select_adjpval_columns(tidy_df, methods)
  if (length(adjpval_cols) == 0) {
    stop("tidy_df has no {scope}.adjPval_{method} columns matching `methods`; run apply_pvalue_adjustments() first.")
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
#' @param alphas,fc_thresholds,methods Passed to
#'   [count_significant_specifications()] - `methods` should be the post-hoc
#'   specification grid's own `adjustment_method` values, so `n_evaluated`
#'   reflects only the intended post-hoc grid rather than every adjustment
#'   method `tidy_df` happens to carry (see the module-level comment above).
#'
#' @return A tibble with columns `gs.name`, `condition`, `time`,
#'   `n_significant` (running total of significant-specification counts),
#'   and `n_evaluated` (running total of specifications evaluated).
accumulate_robustness_counts <- function(accumulator, tidy_df,
                                         alphas        = c(0.0001, 0.001, 0.01, 0.05, 0.1),
                                         fc_thresholds  = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
                                         methods         = NULL) {

  adjpval_cols    <- select_adjpval_columns(tidy_df, methods)
  n_posthoc_specs <- length(adjpval_cols) * length(alphas) * length(fc_thresholds)

  contribution <- tibble::tibble(
    gs.name       = tidy_df$gs.name,
    condition     = tidy_df$condition,
    time          = tidy_df$time,
    n_significant = count_significant_specifications(tidy_df, alphas, fc_thresholds, methods),
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

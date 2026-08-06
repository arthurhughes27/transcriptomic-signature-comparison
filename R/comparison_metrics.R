# =============================================================================
# Stage 1: dearseq vs QuSAGE baseline comparison metrics
# =============================================================================
# Summarises agreement between two DGSA methods' results for the same
# comparisons and gene sets - the Chapter 2, Section 2.1.3 "Stage 1"
# reanalysis (compare an alternative DGSA method's results to the original
# study's). This is deliberately separate from R/robustness_metrics.R, which
# computes the specification-analysis robustness metric pi_{g,v,j} (a single
# method averaged across many specifications), not a two-method comparison.
# =============================================================================

#' Cohen's kappa for two binary vectors
#'
#' @param x,y Logical (or coercible-to-logical) vectors of equal length.
#'
#' @return A single numeric value (NA if undefined, e.g. no non-missing
#'   pairs, or a degenerate case with zero expected agreement variance).
#' @keywords internal
cohens_kappa <- function(x, y) {
  x <- as.logical(x)
  y <- as.logical(y)
  keep <- !is.na(x) & !is.na(y)
  x <- x[keep]
  y <- y[keep]
  n <- length(x)
  if (n == 0) return(NA_real_)

  po <- mean(x == y)
  px <- mean(x)
  py <- mean(y)
  pe <- px * py + (1 - px) * (1 - py)

  if (pe == 1) return(NA_real_)
  (po - pe) / (1 - pe)
}

#' Reshape a combined multi-method results table to one row per comparison x
#' gene set, with separate columns per method
#'
#' @param results_df Tidy results tibble containing multiple methods' rows
#'   (e.g. `dplyr::bind_rows()` of two [build_tidy_dgsa_results()] outputs),
#'   with a `method` column.
#' @param value_cols Character vector of column names to pivot wide by
#'   method (e.g. "rawPval", "fc.score", an adjusted p-value column).
#'
#' @return A tibble with one row per (comparison, condition, time,
#'   gs.name), and `{value_col}_{method}` columns for every combination.
pivot_methods_wide <- function(results_df, value_cols) {
  results_df |>
    dplyr::select(comparison, condition, time, gs.name, method, dplyr::all_of(value_cols)) |>
    tidyr::pivot_wider(
      names_from  = method,
      values_from = dplyr::all_of(value_cols),
      names_glue  = "{.value}_{method}"
    )
}

#' Compute dearseq vs QuSAGE concordance metrics
#'
#' Summarises agreement between two DGSA methods' results for the same
#' comparisons and gene sets: Spearman rank correlation of raw p-values and
#' activation/fold-change scores, and agreement on significance calls at a
#' given adjusted-p-value threshold (percent agreement, Cohen's kappa, and
#' counts of both/either/neither significant).
#'
#' @param results_df Tidy results tibble containing both methods' rows.
#' @param methods Character vector of exactly 2 method labels to compare
#'   (must match values in `results_df$method`).
#' @param adj_pval_col Adjusted p-value column to use for significance
#'   calls (see [apply_pvalue_adjustments()] for the naming convention).
#' @param alpha Significance threshold applied to `adj_pval_col`.
#' @param score_col Score column to correlate ("activation.score" or
#'   "fc.score").
#'
#' @return A list with `method1`/`method2` (the two method labels, defining
#'   which side `n_sig_method1_only`/`n_sig_method2_only` refer to),
#'   `by_comparison` (a tibble, one row per comparison, with `n_gene_sets`,
#'   `pval_cor`, `score_cor`, `pct_agree`, `kappa`, `n_both_sig`,
#'   `n_sig_method1_only`, `n_sig_method2_only`, `n_neither_sig`), and
#'   `overall` (the same metrics pooled across all comparisons).
compute_concordance_metrics <- function(results_df,
                                        methods      = c("dearseq", "qusage"),
                                        adj_pval_col  = "global.adjPval_BH",
                                        alpha          = 0.05,
                                        score_col       = "fc.score") {
  stopifnot(length(methods) == 2)

  m1 <- methods[1]
  m2 <- methods[2]

  wide <- pivot_methods_wide(results_df, value_cols = c("rawPval", score_col, adj_pval_col))

  pval_col1  <- paste0("rawPval_", m1)
  pval_col2  <- paste0("rawPval_", m2)
  score_col1 <- paste0(score_col, "_", m1)
  score_col2 <- paste0(score_col, "_", m2)
  adj_col1   <- paste0(adj_pval_col, "_", m1)
  adj_col2   <- paste0(adj_pval_col, "_", m2)

  wide <- wide |>
    dplyr::filter(!is.na(.data[[pval_col1]]), !is.na(.data[[pval_col2]])) |>
    dplyr::mutate(
      sig1 = .data[[adj_col1]] < alpha,
      sig2 = .data[[adj_col2]] < alpha
    )

  summarise_group <- function(df) {
    tibble::tibble(
      n_gene_sets         = nrow(df),
      pval_cor             = suppressWarnings(stats::cor(df[[pval_col1]], df[[pval_col2]],
                                                         method = "spearman", use = "pairwise.complete.obs")),
      score_cor             = suppressWarnings(stats::cor(df[[score_col1]], df[[score_col2]],
                                                          method = "spearman", use = "pairwise.complete.obs")),
      pct_agree               = mean(df$sig1 == df$sig2, na.rm = TRUE),
      kappa                     = cohens_kappa(df$sig1, df$sig2),
      n_both_sig                  = sum(df$sig1 & df$sig2, na.rm = TRUE),
      n_sig_method1_only             = sum(df$sig1 & !df$sig2, na.rm = TRUE),
      n_sig_method2_only                = sum(!df$sig1 & df$sig2, na.rm = TRUE),
      n_neither_sig                        = sum(!df$sig1 & !df$sig2, na.rm = TRUE)
    )
  }

  by_comparison <- wide |>
    dplyr::group_by(comparison, condition, time) |>
    dplyr::group_modify(~ summarise_group(.x)) |>
    dplyr::ungroup()

  overall <- summarise_group(wide)

  list(method1 = m1, method2 = m2, by_comparison = by_comparison, overall = overall)
}

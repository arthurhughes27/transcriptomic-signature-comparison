# =============================================================================
# Robustness summary tables (Figure/Table 6 suggestion)
# =============================================================================
# Two tabular summaries of the robustness metric pi_{g,v,j}, complementing
# the heatmaps/violin plot: the top-N most robust gene set x comparison
# results, and gene set x comparison results that are significant at
# baseline but NOT robust across the specification grid - the cases where
# the original study's conclusion is most fragile to analytical choices.
# Both use the full geneset.names.descriptions label rather than the short
# gs.name code, and are written out as LaTeX tables (knitr::kable(),
# matching the convention in analysis/descriptive/is2_appendix_descriptives.R)
# ready to \input{} directly into the thesis.
# =============================================================================

#' Format a numeric vector in scientific notation with 3 significant figures
#'
#' Used for the `Robustness` column instead of simple rounding: robustness
#' values can be very small (e.g. a gene set significant under only 1 of
#' several thousand post-hoc specifications) while still corresponding to
#' a significant baseline p-value, and `round(x, 3)` collapses those to a
#' misleading `0`.
#'
#' @param x Numeric vector.
#'
#' @return Character vector, e.g. `0.0005678` -> `"5.68e-04"`.
format_scientific <- function(x) {
  formatC(x, format = "e", digits = 2)
}

#' Attach full gene-set name/description labels to a robustness-like table
#'
#' @param df Tibble with a `gs.name` column.
#' @param genesets Gene set list (as in data/BTM_processed.rds /
#'   data/BG3M_processed.rds), with `geneset.names` and
#'   `geneset.names.descriptions`.
#'
#' @return `df` with a `gs.label` column added (the full
#'   `geneset.names.descriptions` entry for that gene set).
attach_geneset_full_names <- function(df, genesets) {
  name_lookup <- tibble::tibble(
    gs.name  = genesets$geneset.names,
    gs.label = genesets$geneset.names.descriptions
  )
  dplyr::left_join(df, name_lookup, by = "gs.name")
}

#' Build the top-N most robust gene set x comparison table
#'
#' @param robustness_df Output of [compute_robustness_metric()] (`gs.name`,
#'   `condition`, `time`, `robustness`).
#' @param genesets Gene set list, see [attach_geneset_full_names()].
#' @param n Number of rows to keep (highest robustness first).
#'
#' @return A tibble with columns `Gene set`, `Vaccine`, `Timepoint`,
#'   `Robustness`, ready for [save_latex_table()].
build_top_robust_table <- function(robustness_df, genesets, n = 20) {
  robustness_df |>
    attach_geneset_full_names(genesets) |>
    dplyr::arrange(dplyr::desc(robustness), gs.label) |>
    utils::head(n) |>
    dplyr::transmute(
      `Gene set`  = gs.label,
      `Vaccine`    = as.character(condition),
      `Timepoint`  = as.character(time),
      # Raw notation (not scientific): the most-robust results are, by
      # construction, close to 1, so unlike the significant-but-non-robust
      # table there's no risk of a tiny value collapsing to a misleading
      # "0" here.
      `Robustness` = round(robustness, 3)
    )
}

#' Build the "significant at baseline but not robust" table
#'
#' Gene set x comparison x method results that were called significant at
#' the baseline specification (see [join_robustness_baseline()]) but have
#' robustness below `robustness_threshold` - the cases where the original
#' conclusion is most sensitive to analytical choices.
#'
#' @param robustness_df Output of [compute_robustness_metric()].
#' @param baseline_df One or more methods' baseline-specification results,
#'   see [join_robustness_baseline()].
#' @param genesets Gene set list, see [attach_geneset_full_names()].
#' @param pval_col,alpha Baseline significance definition (adjusted
#'   p-value column and threshold).
#' @param robustness_threshold Rows with `robustness < robustness_threshold`
#'   are included (ordered least robust first).
#' @param n Optional cap on the number of rows kept (NULL keeps all
#'   matching rows).
#'
#' @return A tibble with columns `Gene set`, `Vaccine`, `Timepoint`,
#'   `Method`, `Baseline p-value`, `Robustness`, ready for
#'   [save_latex_table()].
build_significant_nonrobust_table <- function(robustness_df, baseline_df, genesets,
                                              pval_col            = "withinTime.adjPval_BH",
                                              alpha                = 0.05,
                                              robustness_threshold = 0.5,
                                              n                     = NULL) {
  joined <- join_robustness_baseline(robustness_df, baseline_df)
  stopifnot(pval_col %in% colnames(joined))

  result <- joined |>
    dplyr::filter(!is.na(.data[[pval_col]]), !is.na(robustness)) |>
    dplyr::filter(.data[[pval_col]] <= alpha, robustness < robustness_threshold) |>
    attach_geneset_full_names(genesets) |>
    dplyr::arrange(robustness, .data[[pval_col]])

  if (!is.null(n)) result <- utils::head(result, n)

  result |>
    dplyr::transmute(
      `Gene set`         = gs.label,
      `Vaccine`           = as.character(condition),
      `Timepoint`         = as.character(time),
      `Method`            = method,
      `Baseline p-value` = formatC(.data[[pval_col]], format = "e", digits = 2),
      `Robustness`        = format_scientific(robustness)
    )
}

#' Save a table as a LaTeX-ready .tex file (and a .csv for reference)
#'
#' Mirrors the convention in
#' analysis/descriptive/is2_appendix_descriptives.R: `knitr::kable()` with
#' `format = "latex"`, `booktabs = TRUE`, `longtable = TRUE`.
#'
#' @param df Tibble to render (typically the output of
#'   [build_top_robust_table()] or [build_significant_nonrobust_table()]).
#' @param path_tex Path to write the `.tex` file to.
#' @param caption,label LaTeX table caption and `\label{}`.
#' @param path_csv Optional path to also write a plain `.csv` copy of `df`
#'   (NULL to skip).
#'
#' @return `path_tex`, invisibly.
save_latex_table <- function(df, path_tex, caption, label, path_csv = NULL) {
  if (!is.null(path_csv)) {
    utils::write.csv(df, file = path_csv, row.names = FALSE)
  }

  tex <- knitr::kable(
    df,
    format    = "latex",
    booktabs  = TRUE,
    longtable = TRUE,
    caption   = caption,
    label     = label
  )

  writeLines(tex, con = path_tex)
  invisible(path_tex)
}

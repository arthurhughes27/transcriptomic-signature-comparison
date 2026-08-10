# =============================================================================
# Baseline significance comparison: dearseq vs QuSAGE (Figure 5)
# =============================================================================
# Illustrates how the single largest analytical lever - DGSA method choice
# - moves results, using the baseline specification only (not the full
# robustness grid): for each vaccine x timepoint comparison and method,
# the percentage of gene sets called significant at baseline (withinTime
# BH-adjusted p <= 0.05 by default). Two dodged bars (dearseq/qusage) per
# vaccine, faceted by timepoint - same x-axis layout convention as
# R/robustness_heatmaps.R (order_robustness_comparisons()).
# =============================================================================

#' Compute the percentage of gene sets significant at baseline, per
#' comparison x method
#'
#' Percentages are of gene sets actually evaluable for that vaccine x
#' timepoint x method (`n_evaluated`), not a fixed gene-set-universe size:
#' a handful of gene set x comparison combinations fail for one method
#' (e.g. dearseq's occasional negative variance-weight estimates), so the
#' evaluable count can differ slightly between methods/comparisons.
#'
#' @param baseline_df Baseline-specification results (see
#'   [load_baseline_results_from_raw()]), with `condition`, `time`,
#'   `method`, and `pval_col`.
#' @param pval_col Adjusted p-value column defining significance (see
#'   [apply_pvalue_adjustments()] for naming). Defaults to the same
#'   `withinTime.adjPval_BH` convention used by the robustness tables.
#' @param alpha Significance threshold applied to `pval_col`.
#'
#' @return A tibble with one row per comparison x method: `condition`,
#'   `time`, `method`, `n_significant`, `n_evaluated`, `pct_significant`
#'   (0-100).
compute_pct_significant_by_comparison <- function(baseline_df,
                                                   pval_col = "withinTime.adjPval_BH",
                                                   alpha     = 0.05) {
  stopifnot(pval_col %in% colnames(baseline_df))

  baseline_df |>
    dplyr::filter(!is.na(.data[[pval_col]])) |>
    dplyr::group_by(condition, time, method) |>
    dplyr::summarise(
      n_significant = sum(.data[[pval_col]] <= alpha),
      n_evaluated   = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(pct_significant = 100 * n_significant / n_evaluated)
}

#' Plot the percentage of gene sets significant at baseline, dearseq vs
#' QuSAGE (Figure 5)
#'
#' @param pct_df Output of [compute_pct_significant_by_comparison()].
#' @param conditions_order Vaccine ordering (see [default_conditions_order()]).
#' @param method_colors Optional colours for the methods present in
#'   `pct_df$method` (see [assign_colours()]).
#' @param value_labels If TRUE (default), annotate each bar with its exact
#'   percentage and evaluated gene-set count.
#'
#' @return A ggplot object, faceted by timepoint.
plot_baseline_significance_comparison <- function(pct_df,
                                                   conditions_order = default_conditions_order(),
                                                   method_colors     = NULL,
                                                   value_labels       = TRUE) {

  plot_data <- order_robustness_comparisons(pct_df, conditions_order)

  methods            <- sort(unique(as.character(plot_data$method)))
  method_colour_map  <- assign_colours(methods, method_colors, palette = "Set1")

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = condition, y = pct_significant, fill = method)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.6) +
    ggplot2::scale_fill_manual(name = "Method", values = method_colour_map) +
    ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggh4x::facet_grid2(
      cols = ggplot2::vars(time), scales = "free_x", space = "free_x",
      labeller = ggplot2::labeller(time = function(x) paste0("Day ", x)),
      strip = day_facet_strip(plot_data$time)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor  = ggplot2::element_blank(),
      strip.text         = ggplot2::element_text(face = "bold"),
      panel.spacing.x     = grid::unit(14, "pt"),
      axis.text.x         = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::labs(
      x     = "Vaccine",
      y     = "Significant gene sets (%)",
      title = "Baseline significance by DGSA method"
    )

  if (value_labels) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.0f%%\n(n=%d)", pct_significant, n_evaluated)),
      position = ggplot2::position_dodge(width = 0.7),
      vjust = -0.2, size = 2.2
    )
  }

  p
}

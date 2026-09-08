# =============================================================================
# Top-robustness gene set bar charts (visual companion to the top-N table)
# =============================================================================
# A visual complement to build_top_robust_table() (R/robustness_tables.R,
# Table 1 in 06_robustness_tables.R's output) - rather than one global
# top-N ranking across every comparison, this shows the top N most robust
# gene sets WITHIN each vaccine x timepoint comparison specifically, as a
# grouped bar chart: vaccines on the x-axis, one small cluster of bars per
# vaccine (its top N gene sets by robustness), gene set names labelled
# above their bar at an angle.
#
# One figure per timepoint (stacked vertically via patchwork rather than
# facetted side by side, unlike plot_robustness_heatmap_aggregate()/
# plot_robustness_violin() - deliberately mirrors the circos plots'
# top-to-bottom day layout, analysis/reanalysis/02_plot_circos.R), each
# with a single day_facet_strip() (R/plot_helpers.R) spanning the whole
# plot width as a shaded, day-coloured title bar - the ggplot-native
# equivalent of the circos plots' plot_row_annotation() coloured label box.
# =============================================================================

#' Top-N most robust gene sets within each vaccine x timepoint comparison
#'
#' Unlike [build_top_robust_table()] (one global top-N ranking across every
#' comparison), this ranks WITHIN each (condition, time) group separately,
#' so every comparison contributes its own top N regardless of how its
#' robustness values compare to other comparisons'.
#'
#' @param robustness_df Output of [compute_robustness_metric()] (`gs.name`,
#'   `condition`, `time`, `robustness`).
#' @param genesets Gene set list, see [attach_geneset_full_names()].
#' @param conditions_order Vaccine ordering (see [default_conditions_order()]).
#' @param n Number of gene sets to keep per comparison.
#'
#' @return `robustness_df` restricted to the top `n` rows per (condition,
#'   time) (ties broken arbitrarily), with `gs.label` (full name, see
#'   [attach_geneset_full_names()]) and `rank` (1 = most robust) columns
#'   added, and `condition`/`time` turned into factors (see
#'   [order_robustness_comparisons()]).
build_top_robust_by_comparison <- function(robustness_df, genesets,
                                           conditions_order = default_conditions_order(),
                                           n = 5) {
  robustness_df |>
    attach_geneset_full_names(genesets) |>
    order_robustness_comparisons(conditions_order) |>
    dplyr::group_by(condition, time) |>
    dplyr::slice_max(robustness, n = n, with_ties = FALSE) |>
    dplyr::mutate(rank = dplyr::row_number(dplyr::desc(robustness))) |>
    dplyr::ungroup()
}

#' Build one timepoint's top-robustness bar chart
#'
#' @param day_data [build_top_robust_by_comparison()]'s output, pre-filtered
#'   to a single `time` level.
#' @param vaccine_colour_map Named vector, vaccine name -> hex colour (see
#'   [plot_top_robust_barcharts()]).
#'
#' @return A ggplot object: one dodged bar per gene set, grouped by vaccine,
#'   with a shaded day-coloured title strip ([day_facet_strip()]).
#' @keywords internal
plot_top_robust_barchart_day <- function(day_data, vaccine_colour_map) {
  dw <- 0.85

  ggplot2::ggplot(day_data, ggplot2::aes(x = condition, y = robustness)) +
    ggplot2::geom_col(
      ggplot2::aes(group = rank, fill = condition),
      position = ggplot2::position_dodge2(width = dw, padding = 0.1),
      colour = "black", linewidth = 0.3
    ) +
    ggplot2::geom_text(
      ggplot2::aes(group = rank, label = gs.label),
      position = ggplot2::position_dodge2(width = dw, padding = 0.1),
      angle = 60, hjust = 0, vjust = 0.3, size = 3
    ) +
    ggplot2::scale_fill_manual(values = vaccine_colour_map, guide = "none") +
    ggh4x::facet_grid2(
      cols = ggplot2::vars(time),
      labeller = ggplot2::labeller(time = function(x) paste0("Day ", x)),
      strip = day_facet_strip(day_data$time)
    ) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.25)) +
    ggplot2::coord_cartesian(ylim = c(0, 1), clip = "off") +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.border      = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.6),
      strip.text         = ggplot2::element_text(face = "bold", size = 16),
      axis.text.x         = ggplot2::element_text(size = 12),
      axis.title           = ggplot2::element_text(size = 14),
      plot.margin            = grid::unit(c(60, 10, 5, 5), "pt")
    ) +
    ggplot2::labs(x = NULL, y = "Signal-robustness")
}

#' Plot the top-N-per-comparison robustness bar charts, one per timepoint
#'
#' Three (or however many timepoints are present) bar charts stacked
#' vertically via [patchwork::wrap_plots()] - see the module-level comment
#' above for why this uses a vertical day-per-figure layout rather than the
#' side-by-side day facets used elsewhere
#' (R/robustness_heatmaps.R/R/robustness_violin.R).
#'
#' @param robustness_df Output of [compute_robustness_metric()].
#' @param genesets Gene set list, see [attach_geneset_full_names()].
#' @param conditions_order Vaccine ordering (see [default_conditions_order()]).
#' @param vaccine_colors Per-vaccine hex colours, aligned to
#'   `conditions_order` (default [default_condition_colors()] - the same
#'   palette used for `condition.colour` throughout the rest of the
#'   repository, e.g. the circos plots/heatmap comparison).
#' @param n Number of gene sets per comparison (default 5).
#'
#' @return A single patchwork object, one panel per timepoint (top to
#'   bottom, ascending day order).
plot_top_robust_barcharts <- function(robustness_df, genesets,
                                      conditions_order = default_conditions_order(),
                                      vaccine_colors = default_condition_colors(),
                                      n = 5) {
  top_data            <- build_top_robust_by_comparison(robustness_df, genesets, conditions_order, n = n)
  vaccine_colour_map  <- assign_colours(conditions_order, vaccine_colors)

  day_levels <- levels(top_data$time)

  day_plots <- lapply(day_levels, function(d) {
    plot_top_robust_barchart_day(dplyr::filter(top_data, time == d), vaccine_colour_map)
  })

  patchwork::wrap_plots(day_plots, ncol = 1) +
    patchwork::plot_annotation(
      title = sprintf("Top %d most robust gene sets per vaccine, by timepoint", n),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 22, face = "bold", hjust = 0.5))
    )
}

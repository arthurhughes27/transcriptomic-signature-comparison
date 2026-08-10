# =============================================================================
# Robustness distribution violin plot (Figure 4)
# =============================================================================
# Shows the full distribution of robustness across the 258 gene sets for
# each comparison, rather than collapsing it to a single summary (as the
# aggregate-level heatmap does). Same x-axis layout as
# R/robustness_heatmaps.R's heatmaps: comparisons split into facets by
# timepoint (so where timepoints change is unambiguous) and ordered within
# each facet by vaccine (default_conditions_order()) - reuses
# order_robustness_comparisons() from that file so the two figures' column
# ordering can never drift apart. Each day's panel gets a black outline
# (panel.border) to further separate the day blocks visually.
# =============================================================================

#' Plot the distribution of per-gene-set robustness for each comparison
#'
#' One violin per comparison (vaccine x timepoint), summarising the
#' distribution of robustness across all gene sets for that comparison,
#' with an inner boxplot for the median/IQR. Comparisons are laid out
#' exactly as in [plot_robustness_heatmap_aggregate()] /
#' [plot_robustness_heatmap_genesets()].
#'
#' @param robustness_df Tibble with one row per gene set x comparison
#'   (`condition`, `time`, `robustness`) - i.e. not pre-aggregated by
#'   gene-set aggregate, unlike [summarise_robustness_by_aggregate()]'s
#'   input.
#' @param conditions_order Vaccine ordering (see [default_conditions_order()]).
#' @param fill_colour Violin/boxplot accent colour.
#'
#' @return A ggplot object, faceted by timepoint.
plot_robustness_violin <- function(robustness_df,
                                   conditions_order = default_conditions_order(),
                                   fill_colour = "#238b45") {

  plot_data <- order_robustness_comparisons(robustness_df, conditions_order)

  ggplot2::ggplot(plot_data, ggplot2::aes(x = condition, y = robustness)) +
    ggplot2::geom_violin(fill = fill_colour, colour = fill_colour, alpha = 0.6, scale = "width", trim = TRUE) +
    ggplot2::geom_boxplot(width = 0.08, outlier.shape = NA, fill = "white", colour = "grey30", alpha = 0.8) +
    ggh4x::facet_grid2(
      cols = ggplot2::vars(time), scales = "free_x", space = "free_x",
      labeller = ggplot2::labeller(time = function(x) paste0("Day ", x)),
      strip = day_facet_strip(plot_data$time)
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor  = ggplot2::element_blank(),
      panel.border       = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.6),
      strip.text          = ggplot2::element_text(face = "bold", size = 12),
      panel.spacing.x      = grid::unit(14, "pt"),
      axis.text.x          = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title = element_text(size = 20),
      plot.title = element_text(size = 20)
    ) +
    ggplot2::labs(
      x     = "Vaccine",
      y     = "Signal-robustness",
      title = "Distribution of signal-robustness across gene sets, by comparison"
    )
}

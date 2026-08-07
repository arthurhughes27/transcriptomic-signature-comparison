# =============================================================================
# Robustness heatmaps
# =============================================================================
# Two-tier visual summary of the robustness metric pi_{g,v,j} (Chapter 2,
# Section 2.2.4): an aggregate-level heatmap (main figure) and a gene-set-
# level heatmap (supplementary figure), both structured the same way -
# columns are comparisons, split into facets by timepoint (so where
# timepoints change is unambiguous) and ordered within each facet by
# vaccine (the existing default_conditions_order()); colour runs white (0)
# to green (1).
# =============================================================================

#' Join gene-set aggregate metadata onto a robustness table
#'
#' @param robustness_df Tibble with a `gs.name` column, as returned by
#'   [compute_robustness_metric()].
#' @param genesets Gene set list (as in data/BTM_processed.rds /
#'   data/BG3M_processed.rds), with `geneset.names` and `geneset.aggregates`.
#'
#' @return `robustness_df` with a `gs.aggregate` column added, matched by
#'   `gs.name`, retaining `genesets$geneset.aggregates`'s factor levels/order.
join_geneset_aggregates <- function(robustness_df, genesets) {
  aggregate_lookup <- tibble::tibble(
    gs.name      = genesets$geneset.names,
    gs.aggregate = genesets$geneset.aggregates
  )
  dplyr::left_join(robustness_df, aggregate_lookup, by = "gs.name")
}

#' Order a robustness table's condition/time columns for consistent plotting
#'
#' Conditions absent from `conditions_order` are appended (with a warning)
#' rather than silently dropped, matching [build_tidy_dgsa_results()]'s
#' handling of the same situation.
#'
#' @param robustness_df Tibble with `condition` and `time` columns.
#' @param conditions_order Character vector giving the desired vaccine
#'   ordering (see [default_conditions_order()]).
#'
#' @return `robustness_df` with `condition` turned into a factor levelled by
#'   `conditions_order` (restricted to the conditions actually present) and
#'   `time` turned into a factor levelled by ascending numeric day.
order_robustness_comparisons <- function(robustness_df, conditions_order = default_conditions_order()) {
  observed_conditions <- unique(as.character(robustness_df$condition))
  unrecognised <- setdiff(observed_conditions, conditions_order)
  if (length(unrecognised) > 0) {
    warning("The following conditions are not in conditions_order and will be appended: ",
            paste(unrecognised, collapse = ", "))
  }
  all_levels     <- c(conditions_order, unrecognised)
  present_levels <- intersect(all_levels, observed_conditions)

  robustness_df |>
    dplyr::mutate(
      condition = factor(as.character(condition), levels = present_levels),
      time      = factor(time, levels = sort(unique(as.numeric(as.character(time)))))
    )
}

#' Summarise robustness at the gene-set-aggregate level
#'
#' For each gene-set aggregate x comparison, computes the unweighted mean of
#' the per-gene-set robustness scores - every gene set counts equally,
#' regardless of how many specifications happened to be evaluable for it
#' (see `n_evaluated` on the input if you want to investigate that
#' separately).
#'
#' @param robustness_df Output of [join_geneset_aggregates()] (must have
#'   `gs.aggregate`, `condition`, `time`, `robustness`).
#'
#' @return A tibble with one row per aggregate x comparison: `gs.aggregate`,
#'   `condition`, `time`, `mean_robustness`, `n_gene_sets`.
summarise_robustness_by_aggregate <- function(robustness_df) {
  robustness_df |>
    dplyr::group_by(gs.aggregate, condition, time) |>
    dplyr::summarise(
      mean_robustness = mean(robustness, na.rm = TRUE),
      n_gene_sets      = dplyr::n(),
      .groups = "drop"
    )
}

# Shared day-facet + colour-scale + theme layers for both heatmaps below.
robustness_heatmap_layers <- function(low_colour, high_colour) {
  list(
    ggplot2::facet_grid(
      . ~ time, scales = "free_x", space = "free_x",
      labeller = ggplot2::labeller(time = function(x) paste0("Day ", x))
    ),
    ggplot2::scale_fill_gradient(name = "Robustness", low = low_colour, high = high_colour, limits = c(0, 1)),
    ggplot2::theme_minimal(),
    ggplot2::theme(
      panel.grid       = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey90", colour = NA),
      strip.text        = ggplot2::element_text(face = "bold"),
      axis.text.x        = ggplot2::element_text(angle = 45, hjust = 1)
    )
  )
}

#' Plot the aggregate-level robustness heatmap (main figure)
#'
#' Rows are gene-set aggregates (existing BTM factor order, top = first
#' level); columns are comparisons, split into facets by timepoint and
#' ordered within each facet by vaccine. Cell colour is the unweighted mean
#' robustness across the aggregate's gene sets (white = 0, green = 1).
#' Comparisons that don't exist (a vaccine missing a timepoint) are simply
#' absent, not shown as blank/zero cells.
#'
#' @param robustness_df Output of [join_geneset_aggregates()].
#' @param conditions_order Vaccine ordering (see [default_conditions_order()]).
#' @param low_colour,high_colour Gradient endpoints.
#'
#' @return A ggplot object.
plot_robustness_heatmap_aggregate <- function(robustness_df,
                                              conditions_order = default_conditions_order(),
                                              low_colour = "white",
                                              high_colour = "#238b45") {

  ordered_df  <- order_robustness_comparisons(robustness_df, conditions_order)
  plot_data   <- summarise_robustness_by_aggregate(ordered_df)
  agg_levels  <- levels(ordered_df$gs.aggregate)

  ggplot2::ggplot(plot_data, ggplot2::aes(x = condition, y = gs.aggregate, fill = mean_robustness)) +
    ggplot2::geom_tile(colour = "grey85") +
    ggplot2::scale_y_discrete(limits = rev(agg_levels)) +
    robustness_heatmap_layers(low_colour, high_colour) +
    ggplot2::labs(x = "Vaccine", y = "Gene set aggregate", title = "Robustness by gene-set aggregate")
}

#' Plot the gene-set-level robustness heatmap (supplementary figure)
#'
#' Same comparison layout as [plot_robustness_heatmap_aggregate()], but with
#' one row per gene set instead of per aggregate: rows are grouped by
#' aggregate (existing BTM factor order), then sorted alphabetically by gene
#' set within an aggregate. A colour strip to the left of the heatmap shows
#' each row's aggregate, using [default_aggregate_colors()] (the same
#' palette used for the circos plots).
#'
#' @param robustness_df Output of [join_geneset_aggregates()].
#' @param conditions_order Vaccine ordering (see [default_conditions_order()]).
#' @param aggregate_colors Aggregate colour palette (see
#'   [default_aggregate_colors()]), aligned positionally to
#'   `levels(robustness_df$gs.aggregate)`.
#' @param drop_null_gene_sets If TRUE (default), gene sets with robustness 0
#'   (or NA) in every comparison are excluded before plotting, and the
#'   number dropped is reported via `message()`.
#' @param low_colour,high_colour Gradient endpoints.
#' @param strip_width Relative width of the annotation strip vs. the main
#'   heatmap panel.
#'
#' @return A patchwork object (aggregate strip + main heatmap, legends
#'   collected together).
plot_robustness_heatmap_genesets <- function(robustness_df,
                                             conditions_order = default_conditions_order(),
                                             aggregate_colors = default_aggregate_colors(),
                                             drop_null_gene_sets = TRUE,
                                             low_colour = "white",
                                             high_colour = "#238b45",
                                             strip_width = 1) {

  plot_data <- order_robustness_comparisons(robustness_df, conditions_order)

  if (drop_null_gene_sets) {
    always_null <- plot_data |>
      dplyr::group_by(gs.name) |>
      dplyr::summarise(all_null = all(is.na(robustness) | robustness == 0), .groups = "drop") |>
      dplyr::filter(all_null) |>
      dplyr::pull(gs.name)

    if (length(always_null) > 0) {
      message(sprintf(
        "Dropping %d gene set(s) with robustness = 0 in every comparison.",
        length(always_null)
      ))
      plot_data <- dplyr::filter(plot_data, !gs.name %in% always_null)
    }
  }

  gene_set_order <- plot_data |>
    dplyr::distinct(gs.name, gs.aggregate) |>
    dplyr::arrange(gs.aggregate, gs.name) |>
    dplyr::pull(gs.name) |>
    as.character()

  plot_data <- dplyr::mutate(plot_data, gs.name = factor(as.character(gs.name), levels = gene_set_order))

  aggregate_levels     <- levels(plot_data$gs.aggregate)
  aggregate_colour_map <- assign_colours(aggregate_levels, aggregate_colors, palette = "Spectral")

  annotation_data <- dplyr::distinct(plot_data, gs.name, gs.aggregate)

  p_strip <- ggplot2::ggplot(annotation_data, ggplot2::aes(x = 1, y = gs.name, fill = gs.aggregate)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(name = "Aggregate", values = aggregate_colour_map, drop = FALSE) +
    ggplot2::scale_y_discrete(limits = rev(gene_set_order)) +
    ggplot2::theme_void()

  p_main <- ggplot2::ggplot(plot_data, ggplot2::aes(x = condition, y = gs.name, fill = robustness)) +
    ggplot2::geom_tile() +
    ggplot2::scale_y_discrete(limits = rev(gene_set_order)) +
    robustness_heatmap_layers(low_colour, high_colour) +
    ggplot2::labs(x = "Vaccine", y = NULL, title = "Robustness by gene set") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6))

  patchwork::wrap_plots(p_strip, p_main, ncol = 2, widths = c(strip_width, 20), guides = "collect")
}

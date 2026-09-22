# =============================================================================
# Robustness heatmaps
# =============================================================================
# Two-tier visual summary of the robustness metric pi_{g,v,j} (Chapter 2,
# Section 2.2.4): an aggregate-level heatmap (main figure) and a gene-set-
# level heatmap (supplementary figure), both structured the same way -
# columns are comparisons, split into facets by timepoint (so where
# timepoints change is unambiguous - each "Day X" facet strip is coloured
# per day, see R/plot_helpers.R's assign_day_colours()/day_facet_strip())
# and ordered within each facet by vaccine (the existing
# default_conditions_order()); colour runs white (0) to green (1).
# =============================================================================

# Square-root transformation for the robustness colour scale.
# This expands low robustness values while preserving the 0-1 range
# and the legend labels on the natural scale.
power_trans <- function(power = 0.5) {
  scales::trans_new(
    name = paste0("power_", power),
    transform = function(x) x^power,
    inverse = function(x) x^(1 / power)
  )
}

#' Join gene-set aggregate metadata onto a robustness table
#'
#' @param robustness_df Tibble with a `gs.name` column, as returned by
#'   [compute_robustness_metric()].
#' @param genesets Gene set list (as in data/BTM_processed.rds /
#'   data/BG3M_processed.rds), with `geneset.names`, `geneset.aggregates`,
#'   and `geneset.names.descriptions`.
#'
#' @return `robustness_df` with `gs.aggregate` and `gs.label` columns
#'   added, matched by `gs.name`: `gs.aggregate` retains
#'   `genesets$geneset.aggregates`'s factor levels/order, and `gs.label`
#'   is the full `geneset.names.descriptions` entry (used for display -
#'   `gs.name` remains the short, unique join/sort key).
join_geneset_aggregates <- function(robustness_df, genesets) {
  aggregate_lookup <- tibble::tibble(
    gs.name      = genesets$geneset.names,
    gs.aggregate = genesets$geneset.aggregates,
    gs.label     = genesets$geneset.names.descriptions
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
# Facet strips are coloured per day (assign_day_colours()/day_facet_strip(),
# R/plot_helpers.R) via ggh4x::facet_grid2(), so it's easy to tell at a
# glance where one timepoint's block of columns ends and the next begins;
# panel.spacing.x is widened and each day's panel gets a black outline
# (panel.border) to the same end.
#
# @param colorbar_guide Guide for the Signal-robustness fill scale - lets
#   callers lengthen/resize the colourbar (e.g. the gene-set-level
#   heatmap's larger legend) without duplicating the whole scale
#   definition. Defaults to a plain [ggplot2::guide_colorbar()].
robustness_heatmap_layers <- function(low_colour, high_colour, times, day_colors = NULL,
                                      colorbar_guide = ggplot2::guide_colorbar()) {
  list(
    ggh4x::facet_grid2(
      cols = ggplot2::vars(time), scales = "free_x", space = "free_x",
      labeller = ggplot2::labeller(time = function(x) paste0("Day ", x)),
      strip = day_facet_strip(times, day_colors)
    ),
    ggplot2::scale_fill_gradient(
      name = "Signal-robustness",
      low = low_colour,
      high = high_colour,
      limits = c(0, 1),
      trans = power_trans(0.5),
      guide = colorbar_guide
    ),
    ggplot2::theme_minimal(),
    ggplot2::theme(
      panel.grid       = ggplot2::element_blank(),
      panel.border     = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.5),
      strip.text        = ggplot2::element_text(face = "bold", size = 12),
      panel.spacing.x    = grid::unit(14, "pt"),
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
    robustness_heatmap_layers(low_colour, high_colour, times = plot_data$time) + ggplot2::scale_fill_gradient(
      name = "Mean signal-robustness",
      low = low_colour,
      high = high_colour,
      limits = c(0, 1),
      trans = power_trans(0.5)
    ) +
    ggplot2::labs(x = "Vaccine", y = "Gene set aggregate", title = "Mean signal-robustness by gene-set aggregate") +
    theme(axis.title = element_text(size = 25),
          plot.title = element_text(size = 23),
          axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 14),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.key.height = unit(1, units = "cm"),
          legend.key.spacing = unit(0.75, 'cm'))
}

#' Build one page of the gene-set-level robustness heatmap
#'
#' @param plot_data Rows for this page only (already filtered).
#' @param gene_set_order This page's gene sets (full `gs.label` names), in
#'   display order.
#' @param aggregate_colour_map Named aggregate -> colour vector, shared
#'   across every page so a given aggregate's colour - and the full set of
#'   legend entries - are identical on every page, whether or not that
#'   aggregate has any gene sets on this particular page (`limits =
#'   names(aggregate_colour_map)` forces every page's legend to list all
#'   of them, not just the ones actually present).
#' @param low_colour,high_colour,strip_width,label_width,label_cutoff See
#'   [plot_robustness_heatmap_genesets()].
#' @param page_number,total_pages This page's 1-indexed position and the
#'   total page count, appended to the plot title as "(page X/N)" so a
#'   printed/exported page can always be placed back in the sequence.
#' @param row_height_cm Exact panel height (in cm) PER ROW, applied via
#'   [ggh4x::force_panelsizes()] to every sub-panel below - see
#'   [plot_robustness_heatmap_genesets()] for why this is forced rather
#'   than left to auto-stretch.
#'
#' @return A patchwork object (aggregate strip + row-label column + main
#'   heatmap, legends collected together).
#' @keywords internal
build_geneset_heatmap_page <- function(plot_data, gene_set_order, aggregate_colour_map,
                                       low_colour, high_colour, strip_width,
                                       label_width, label_cutoff,
                                       page_number, total_pages, row_height_cm) {
  annotation_data <- dplyr::distinct(plot_data, gs.label, gs.aggregate)

  # Shared sizing so the two legends (Aggregate, Signal-robustness) read
  # consistently once collected side by side by patchwork.
  legend_title_size <- 26
  legend_text_size  <- 20

  # Forced identically on every sub-panel below (not just the main heatmap)
  # so their rows stay vertically aligned with each other, and - critically
  # - so a page with fewer gene sets than rows_per_page gets a
  # proportionally SHORTER panel (leaving blank space below on the page)
  # rather than the existing rows being stretched to fill the same total
  # height as a full page. See the "exact dimensions of the rows... exact
  # same between pages" request in plot_robustness_heatmap_genesets().
  panel_height <- grid::unit(length(gene_set_order) * row_height_cm, "cm")

  p_strip <- ggplot2::ggplot(annotation_data, ggplot2::aes(x = 1, y = gs.label, fill = gs.aggregate)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(
      name = "Aggregate", values = aggregate_colour_map,
      limits = names(aggregate_colour_map), drop = FALSE,
      guide = ggplot2::guide_legend(
        keywidth = grid::unit(1, "cm"), keyheight = grid::unit(1, "cm")
      )
    ) +
    ggplot2::scale_y_discrete(limits = rev(gene_set_order)) +
    ggh4x::force_panelsizes(rows = panel_height) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.title = ggplot2::element_text(size = legend_title_size),
      legend.text  = ggplot2::element_text(size = legend_text_size)
    )

  # Row labels (truncated, see truncate_geneset_label() - R/plot_helpers.R)
  # drawn as their own fixed-width geom_text column instead of the main
  # panel's y-axis text: ggplot auto-sizes an axis-text column to fit
  # whatever's actually on that page, so on the old layout a page with
  # longer names left LESS room for the heatmap panel than a page with
  # shorter ones, misaligning the panel's size/position across pages (the
  # "dimensions ... not aligned" request). A dedicated column with a FIXED
  # relative width (label_width, set by patchwork below, same on every
  # page regardless of content) decouples the main panel's geometry from
  # label content entirely.
  label_pt <- 13
  label_mm <- label_pt / (72.27 / 25.4)

  p_labels <- ggplot2::ggplot(
    annotation_data,
    ggplot2::aes(x = 1, y = gs.label, label = truncate_geneset_label(as.character(gs.label), label_cutoff))
  ) +
    ggplot2::geom_text(hjust = 1, size = label_mm) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_y_discrete(limits = rev(gene_set_order)) +
    ggh4x::force_panelsizes(rows = panel_height) +
    ggplot2::theme_void()

  p_main <- ggplot2::ggplot(plot_data, ggplot2::aes(x = condition, y = gs.label, fill = robustness)) +
    ggplot2::geom_tile() +
    ggplot2::scale_y_discrete(limits = rev(gene_set_order)) +
    ggh4x::force_panelsizes(rows = panel_height) +
    robustness_heatmap_layers(
      low_colour, high_colour, times = plot_data$time,
      colorbar_guide = ggplot2::guide_colorbar(
        barheight = grid::unit(16, "cm"), barwidth = grid::unit(1.4, "cm")
      )
    ) +
    ggplot2::labs(
      x = "Vaccine", y = NULL,
      title = sprintf("Signal-robustness by gene set (page %d/%d)", page_number, total_pages)
    ) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_text(size = 32),
                   axis.title.y = ggplot2::element_text(size = 32),
                   plot.title = ggplot2::element_text(size = 34, face = "bold"),
                   axis.text.x = ggplot2::element_text(size = 15),
                   legend.title = ggplot2::element_text(size = legend_title_size),
                   legend.text = ggplot2::element_text(size = legend_text_size),
                  )

  patchwork::wrap_plots(
    p_strip, p_labels, p_main, ncol = 3,
    widths = c(strip_width, label_width, 20), guides = "collect"
  )
}

#' Plot the gene-set-level robustness heatmap (supplementary figure),
#' paginated
#'
#' Same comparison layout as [plot_robustness_heatmap_aggregate()], but with
#' one row per gene set instead of per aggregate: rows are labelled with the
#' full gene-set name/description (`gs.label`, from
#' [join_geneset_aggregates()]), grouped by aggregate (existing BTM factor
#' order), then sorted alphabetically by that full name within an
#' aggregate. A colour strip to the left of each page shows each row's
#' aggregate, using [default_aggregate_colors()] (the same palette used for
#' the circos plots) - every page's aggregate legend lists all aggregates
#' present anywhere in the data (not just this page's), so the legend is
#' identical across pages.
#'
#' With 258 gene sets, a single figure at a readable row height doesn't fit
#' on a page/screen without heavy zooming, so gene sets are split across
#' fixed-size pages (in the existing aggregate + alphabetical order - a
#' page's gene sets are usually, but not always, from the same aggregate,
#' since an aggregate can straddle a page boundary). Every page uses the
#' same row height and the same comparison columns, so pages are directly
#' comparable and consistently sized regardless of how large the aggregate
#' they happen to fall in is.
#'
#' Two things are forced identically across every page, rather than left to
#' ggplot's normal auto-layout, so pages stay visually comparable even when
#' their content differs (varying label lengths; a last page with fewer
#' than `rows_per_page` gene sets):
#'   - the main heatmap panel's row labels are drawn in their own
#'     fixed-relative-width column (`label_width`, `label_cutoff` - see
#'     [truncate_geneset_label()], R/plot_helpers.R) instead of as the
#'     panel's own y-axis text, whose auto-sized width would otherwise vary
#'     with the longest label on that particular page, shifting the main
#'     panel's size/position page to page;
#'   - every sub-panel's height is forced to `rows_per_page_on_this_page *
#'     row_height_cm` via [ggh4x::force_panelsizes()], so the row height
#'     itself is identical on every page - a page with fewer gene sets
#'     (typically the last one) gets a proportionally SHORTER panel with
#'     blank space below on the page, rather than its rows being stretched
#'     to fill the same total panel height as a full page.
#'
#' Save the result with [save_multi_page_pdf()] - the driver scripts
#' (04_specification_heatmap_aggregate.R and analysis/supplementary/
#' specification_heatmap_genesets_supplementary.R) pair `rows_per_page`
#' with `row_height_cm` and an A4-ratio page size so each page prints at a
#' readable row height.
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
#' @param strip_width Relative width of the aggregate-colour annotation
#'   strip vs. the main heatmap panel (the row-label column, `label_width`,
#'   sits between them - see [build_geneset_heatmap_page()]).
#' @param label_width Relative width of the row-label column - tune this
#'   together with `label_cutoff` (longer allowed labels need a wider
#'   column to avoid running into the aggregate strip on their left).
#' @param label_cutoff Passed to [truncate_geneset_label()] - the gene-set
#'   name truncation length. NULL keeps full, untruncated labels (not
#'   recommended here - an unbounded label length reintroduces the
#'   page-to-page column-width risk `label_width`'s fixed width is meant to
#'   avoid).
#' @param rows_per_page Maximum gene sets per page. `Inf` reproduces the
#'   previous single-figure behaviour (one very tall page).
#' @param row_height_cm Exact panel height per row (cm), forced identically
#'   on every page regardless of how many gene sets that page actually has
#'   - see the function-level comment above. Pick this to match the page
#'   size/`rows_per_page` the driver script actually uses (e.g. page height
#'   in cm / `rows_per_page`), so a full page's panel fills the intended
#'   height exactly.
#'
#' @return A list of patchwork objects, one per page.
plot_robustness_heatmap_genesets <- function(robustness_df,
                                             conditions_order = default_conditions_order(),
                                             aggregate_colors = default_aggregate_colors(),
                                             drop_null_gene_sets = TRUE,
                                             low_colour = "white",
                                             high_colour = "#238b45",
                                             strip_width = 1,
                                             label_width = 6,
                                             label_cutoff = 20,
                                             rows_per_page = 45,
                                             row_height_cm = 1.4) {

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
    dplyr::distinct(gs.label, gs.aggregate) |>
    dplyr::arrange(gs.aggregate, gs.label) |>
    dplyr::pull(gs.label) |>
    as.character()

  plot_data <- dplyr::mutate(plot_data, gs.label = factor(as.character(gs.label), levels = gene_set_order))

  aggregate_levels     <- levels(plot_data$gs.aggregate)
  aggregate_colour_map <- assign_colours(aggregate_levels, aggregate_colors, palette = "Spectral")

  page_gene_sets <- if (is.infinite(rows_per_page) || rows_per_page >= length(gene_set_order)) {
    list(gene_set_order)
  } else {
    unname(split(gene_set_order, ceiling(seq_along(gene_set_order) / rows_per_page)))
  }

  total_pages <- length(page_gene_sets)

  lapply(seq_along(page_gene_sets), function(page_number) {
    gene_sets_on_page <- page_gene_sets[[page_number]]
    build_geneset_heatmap_page(
      plot_data             = dplyr::filter(plot_data, as.character(gs.label) %in% gene_sets_on_page),
      gene_set_order         = gene_sets_on_page,
      aggregate_colour_map    = aggregate_colour_map,
      low_colour               = low_colour,
      high_colour               = high_colour,
      strip_width                = strip_width,
      label_width               = label_width,
      label_cutoff               = label_cutoff,
      page_number                 = page_number,
      total_pages                 = total_pages,
      row_height_cm                = row_height_cm
    )
  })
}

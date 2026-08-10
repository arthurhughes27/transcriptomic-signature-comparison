# =============================================================================
# Shared plotting/presentation helpers for DGSA results
# =============================================================================
# Colour assignment and default condition/aggregate configuration shared by
# the dearseq and QuSAGE result-processing driver scripts
# (analysis/reanalysis/process_*_dgsa_results.R), so the two don't drift
# apart. Deliberately separate from R/postprocessing.R, which is scoped to
# statistical postprocessing only (see the notes at the top of that file).
# =============================================================================

#' Assign a named colour vector to a set of labels
#'
#' Uses `provided_colours` if supplied (recycled by name, one per label);
#' otherwise generates colours from a named RColorBrewer palette, ramping if
#' there are more labels than colours in the base palette.
#'
#' @param labels Character vector of labels to assign colours to.
#' @param provided_colours Optional character vector of hex colours, one per
#'   label (same order as `labels`).
#' @param palette RColorBrewer palette name, used when `provided_colours` is
#'   NULL.
#'
#' @return A named character vector of hex colours, named by `labels`.
assign_colours <- function(labels, provided_colours = NULL, palette = "Paired") {
  n <- length(labels)
  if (!is.null(provided_colours)) {
    stopifnot(length(provided_colours) == n)
    return(stats::setNames(provided_colours, labels))
  }
  max_cols <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
  base     <- RColorBrewer::brewer.pal(min(n, max_cols), palette)
  colours  <- if (n > length(base)) grDevices::colorRampPalette(base)(n) else base
  stats::setNames(colours, labels)
}

#' Add condition and gene-set-aggregate colour columns to a tidy DGSA table
#'
#' @param results_df Tidy results tibble (output of
#'   [build_tidy_dgsa_results()]), with `condition` and `gs.aggregate`
#'   columns.
#' @param condition_colors Optional explicit colours for
#'   `levels(results_df$condition)` (or `sort(unique(...))` if not a
#'   factor).
#' @param aggregate_colors Optional explicit colours for the `gs.aggregate`
#'   levels.
#'
#' @return `results_df` with `condition.colour` and `gs.colour` columns
#'   added.
assign_dgsa_colours <- function(results_df, condition_colors = NULL, aggregate_colors = NULL) {
  unique_conditions <- if (is.factor(results_df$condition)) {
    levels(results_df$condition)
  } else {
    sort(unique(results_df$condition))
  }
  unique_aggregates <- if (is.factor(results_df$gs.aggregate)) {
    levels(results_df$gs.aggregate)
  } else {
    sort(unique(results_df$gs.aggregate))
  }

  condition_colour_map <- assign_colours(unique_conditions, condition_colors, palette = "Paired")
  aggregate_colour_map <- assign_colours(unique_aggregates, aggregate_colors, palette = "Spectral")

  results_df |>
    dplyr::mutate(
      condition.colour = condition_colour_map[as.character(condition)],
      gs.colour         = aggregate_colour_map[as.character(gs.aggregate)]
    )
}

#' Default vaccine condition ordering
#'
#' Roughly groups vaccines by pathogen family, for consistent figure
#' ordering across scripts. Shared between the dearseq and QuSAGE
#' result-processing scripts so their outputs use the same factor levels.
default_conditions_order <- function() {
  c(
    "Tuberculosis (RVV)",
    "Varicella Zoster (LV)",
    "Yellow Fever (LV)",
    "Ebola (RVV)",
    "Hepatitis A/B (IN/RP)",
    "HIV (RVV)",
    "Influenza (IN)",
    "Influenza (LV)",
    "Malaria (RP)",
    "Meningococcus (CJ)",
    "Meningococcus (PS)",
    "Pneumococcus (PS)",
    "Smallpox (LV)"
  )
}

#' Default per-condition hex colours, aligned to [default_conditions_order()]
default_condition_colors <- function() {
  c(
    "#b94a73", "#c6aa3c", "#6f71d9", "#64c46a", "#be62c2",
    "#7d973c", "#563382", "#4ea76e", "#bc69b0", "#33d4d1",
    "#bb4c41", "#6a87d3", "#b57736"
  )
}

#' Default per-gene-set-aggregate hex colours (BTM functional categories,
#' including "NA" for genesets with no assigned aggregate)
default_aggregate_colors <- function() {
  c(
    "#7c5fcd", "#57c39d", "#c1121f", "#55c463", "#7082ca",
    "#64a332", "#45aecf", "#df9545", "#b7b238", "#a6b36c",
    "#667328", "#662d2e", "#ff8fa3", "#c05299", "#8f2d56",
    "#adb5bd"
  )
}

#' Default per-DGSA-method hex colours (muted, low-contrast pair), aligned
#' to `c("dearseq", "qusage")`
default_method_colors <- function() {
  c(dearseq = "#7C9CBF", qusage = "#C98B7B")
}

#' Canonical day (timepoint) ordering used for consistent day-colour
#' assignment across every multi-day figure - a superset of whichever
#' days happen to be present in any one analysis/figure.
default_day_order <- function() {
  c(1, 3, 7, 10, 14, 21)
}

#' Default light, "not too loud" per-day hex colours, aligned to
#' [default_day_order()]
default_day_colors <- function() {
  c(
    "#cfe2f3", # Day 1  - light blue
    "#d9ead3", # Day 3  - light green
    "#fce5cd", # Day 7  - light peach
    "#e6d6ec", # Day 10 - light lavender
    "#d4f1f0", # Day 14 - light teal
    "#fff2cc"  # Day 21 - light yellow
  )
}

#' Assign light, consistent per-day colours for figure highlighting
#'
#' Colours are assigned from the fixed [default_day_order()] /
#' [default_day_colors()] pairing, then restricted to the days actually
#' present in `times` - so a given day (e.g. "Day 7") gets the SAME
#' colour in every figure it appears in, regardless of which other days
#' happen to appear alongside it in that particular figure. Days outside
#' `default_day_order()` are appended in ascending order with colours
#' ramped from the base palette.
#'
#' @param times Vector of day values (numeric or coercible) actually
#'   present in the figure being built.
#' @param day_colors Optional explicit colours to use in place of
#'   [default_day_colors()], aligned to `default_day_order()` (plus any
#'   extra days present in `times`).
#'
#' @return A named character vector: names are `"Day X"` labels (ascending
#'   day order, restricted to the days in `times`), values are hex
#'   colours.
assign_day_colours <- function(times, day_colors = NULL) {
  present_days <- sort(unique(as.numeric(as.character(times))))
  base_days    <- default_day_order()
  base_colors  <- if (is.null(day_colors)) default_day_colors() else day_colors

  extra_days <- setdiff(present_days, base_days)
  all_days   <- c(base_days, extra_days)
  all_colors <- if (length(extra_days) > 0) {
    c(base_colors, grDevices::colorRampPalette(base_colors)(length(extra_days)))
  } else {
    base_colors
  }

  full_map <- stats::setNames(all_colors[seq_along(all_days)], paste0("Day ", all_days))
  full_map[paste0("Day ", present_days)]
}

#' Build a ggh4x-themed x-facet strip with light, per-day background
#' colours
#'
#' For use with `ggh4x::facet_grid2()` in place of `ggplot2::facet_grid()`
#' so each "Day X" facet strip gets a distinct (but light) background
#' colour instead of one flat grey, making it easier to tell at a glance
#' where one timepoint's block of columns ends and the next begins.
#'
#' @param times Vector of day values, as they will appear in the plot's
#'   `time` facetting variable.
#' @param day_colors Optional explicit colours (see [assign_day_colours()]).
#'
#' @return A `ggh4x::strip_themed()` object.
day_facet_strip <- function(times, day_colors = NULL) {
  day_colour_map <- assign_day_colours(times, day_colors)
  ggh4x::strip_themed(
    background_x = ggh4x::elem_list_rect(fill = unname(day_colour_map))
  )
}

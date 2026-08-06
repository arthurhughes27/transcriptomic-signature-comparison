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

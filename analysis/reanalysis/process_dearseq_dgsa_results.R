# =============================================================================
# Preprocessing of dearseq DGSA Results
# =============================================================================
# Extracts results from the raw dearseq output list, assembles them into a
# tidy dataframe, assigns colours, and applies multiple p-value corrections.
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────

library(fs)
library(tidyverse)
library(RColorBrewer)

# ── Paths ─────────────────────────────────────────────────────────────────────

p_results_list <- fs::path("output", "results", "reanalysis", "dearseq_dgsa_results_list.rds")
p_results_df   <- fs::path("output", "results", "reanalysis", "dearseq_dgsa_results_processed.rds")
p_data_btm     <- fs::path("data", "BTM_processed.rds")

# ── Load data ─────────────────────────────────────────────────────────────────

dearseq_dgsa_list <- readRDS(p_results_list)
BTM               <- readRDS(p_data_btm)

# =============================================================================
# ANALYSIS CONFIGURATION
# =============================================================================
# Condition ordering, colours, and aggregate colours are defined here so they
# are easy to find and update without touching function logic.

conditions_order <- c(
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

condition_colors <- c(
  "#b94a73", "#c6aa3c", "#6f71d9", "#64c46a", "#be62c2",
  "#7d973c", "#563382", "#4ea76e", "#bc69b0", "#33d4d1",
  "#bb4c41", "#6a87d3", "#b57736"
)

aggregate_colors <- c(
  "#7c5fcd", "#57c39d", "#c1121f", "#55c463", "#7082ca",
  "#64a332", "#45aecf", "#df9545", "#b7b238", "#a6b36c",
  "#667328", "#662d2e", "#ff8fa3", "#c05299", "#8f2d56",
  "#adb5bd"
)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Extract the vaccine name and timepoint stored in each result list entry,
# avoiding fragile regex parsing of list names.
extract_results <- function(results_list) {
  list(
    raw_pvals     = lapply(results_list, function(x) x[["pvals"]][["rawPval"]]),
    scores        = lapply(results_list, function(x)
      list(activation.scores = x[["score"]][["activation.scores"]],
           fc.scores         = x[["score"]][["fc.scores"]])),
    correlations  = lapply(results_list, function(x)
      list(mean.corr = x[["cor"]][["mean.corr"]],
           corr.mean = x[["cor"]][["corr.mean"]])),
    conditions    = sub("^(.*?)\\s+vs.*", "\\1", names(results_list)),
    times         = as.numeric(sub(".*Day\\s+([0-9.]+)$", "\\1", names(results_list)))
  )
}

# Assign a named colour vector to a set of labels, using provided colours if
# supplied, otherwise generating them from a named RColorBrewer palette.
assign_colours <- function(labels, provided_colours = NULL, palette = "Paired") {
  n <- length(labels)
  if (!is.null(provided_colours)) {
    stopifnot(length(provided_colours) == n)
    return(setNames(provided_colours, labels))
  }
  max_cols <- brewer.pal.info[palette, "maxcolors"]
  base     <- brewer.pal(min(n, max_cols), palette)
  colours  <- if (n > length(base)) colorRampPalette(base)(n) else base
  setNames(colours, labels)
}

# Apply a single p-value correction method at three scopes (global,
# within-timepoint, within-comparison) and append the resulting columns.
apply_pval_correction <- function(df, correction_method) {
  df %>%
    mutate(
      !!paste0("global.adjPval_", correction_method) := p.adjust(rawPval, method = correction_method)
    ) %>%
    group_by(time) %>%
    mutate(
      !!paste0("withinTime.adjPval_", correction_method) := p.adjust(rawPval, method = correction_method)
    ) %>%
    group_by(comparison) %>%
    mutate(
      !!paste0("withinComparison.adjPval_", correction_method) := p.adjust(rawPval, method = correction_method)
    ) %>%
    ungroup()
}

# =============================================================================
# MAIN PREPROCESSING FUNCTION
# =============================================================================

#' Preprocess dearseq DGSA results into a tidy dataframe
#'
#' @param results_list   Named list of raw dearseq outputs, one entry per
#'                       vaccine x timepoint comparison.
#' @param genesets       GSA.genesets object with \code{genesets},
#'                       \code{geneset.names}, \code{geneset.descriptions},
#'                       and \code{geneset.aggregates}.
#' @param conditions_order Character vector giving the desired factor ordering
#'                       of vaccine conditions.
#' @param condition_colors Optional named or unnamed character vector of hex
#'                       colours, one per unique condition.
#' @param aggregate_colors Optional named or unnamed character vector of hex
#'                       colours, one per unique geneset aggregate.
#' @param method         String label identifying the analysis method, appended
#'                       as a column (useful when combining results across
#'                       methods downstream).
#'
#' @return A tidy dataframe with one row per geneset x comparison, with
#'         scores, correlations, and multiple p-value corrections appended.

dgsa_comparison_preprocessing <- function(results_list,
                                          genesets,
                                          conditions_order  = NULL,
                                          condition_colors  = NULL,
                                          aggregate_colors  = NULL,
                                          method            = NA_character_) {
  
  # ── Extract components from raw results list ───────────────────────────────
  
  extracted     <- extract_results(results_list)
  raw_pvals     <- extracted$raw_pvals
  scores        <- extracted$scores
  correlations  <- extracted$correlations
  conditions    <- extracted$conditions
  times         <- extracted$times
  
  n_comparisons <- length(raw_pvals)
  n_gene_sets   <- length(raw_pvals[[1]])
  
  # ── Validity checks ────────────────────────────────────────────────────────
  
  if (length(genesets[["geneset.names"]]) != n_gene_sets)
    stop("Number of genesets in BTM does not match number of p-values per comparison.")
  
  if (!all(c(length(scores), length(correlations), length(conditions), length(times)) == n_comparisons))
    stop("Number of comparisons does not match across extracted components.")
  
  # ── Factor ordering for conditions and times ───────────────────────────────
  
  if (!is.null(conditions_order)) {
    unrecognised <- setdiff(conditions, conditions_order)
    if (length(unrecognised) > 0)
      warning("The following conditions are not in conditions_order and will be appended: ",
              paste(unrecognised, collapse = ", "))
    all_levels <- c(conditions_order, unrecognised)
    conditions <- factor(conditions, levels = intersect(all_levels, conditions))
  } else {
    conditions <- factor(conditions, levels = sort(unique(conditions)))
  }
  
  times <- factor(times, levels = unique(times))
  
  # ── Colour assignment ──────────────────────────────────────────────────────
  
  unique_conditions <- levels(conditions)
  unique_aggregates <- if (is.factor(genesets$geneset.aggregates)) {
    levels(genesets$geneset.aggregates)
  } else {
    sort(unique(genesets$geneset.aggregates))
  }
  
  condition_colour_map <- assign_colours(unique_conditions, condition_colors,  palette = "Paired")
  aggregate_colour_map <- assign_colours(unique_aggregates, aggregate_colors,  palette = "Spectral")
  
  # ── Assemble tidy results dataframe ───────────────────────────────────────
  
  comparison_labels <- paste0(as.character(conditions), " - Day ", as.character(times))
  
  results_df <- purrr::map_dfr(seq_along(raw_pvals), function(idx) {
    tibble(
      comparison          = comparison_labels[idx],
      condition           = conditions[idx],
      time                = times[idx],
      condition.colour    = condition_colour_map[as.character(conditions[idx])],
      gs.name             = genesets$geneset.names,
      gs.description      = genesets$geneset.descriptions,
      gs.name.description = paste0(genesets$geneset.names, " - ", genesets$geneset.descriptions),
      gs.aggregate        = genesets$geneset.aggregates,
      gs.colour           = aggregate_colour_map[as.character(genesets$geneset.aggregates)],
      activation.score    = as.numeric(unlist(scores[[idx]][["activation.scores"]])),
      fc.score            = as.numeric(unlist(scores[[idx]][["fc.scores"]])),
      mean.corr           = correlations[[idx]][["mean.corr"]],
      corr.mean           = correlations[[idx]][["corr.mean"]],
      rawPval             = raw_pvals[[idx]],
      method              = method
    )
  })
  
  # ── P-value correction across six methods and three scopes ────────────────
  
  correction_methods <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY")
  
  results_df <- reduce(correction_methods, apply_pval_correction, .init = results_df)
  
  results_df
}

# =============================================================================
# RUN PREPROCESSING
# =============================================================================

results_df <- dgsa_comparison_preprocessing(
  results_list     = dearseq_dgsa_list,
  genesets         = BTM,
  conditions_order = conditions_order,
  condition_colors = condition_colors,
  aggregate_colors = aggregate_colors,
  method           = "dearseq"
)

# ── Save ───────────────────────────────────────────────────────────────────────

saveRDS(results_df, file = p_results_df)
message("Processed results saved to: ", p_results_df)

rm(list = ls())

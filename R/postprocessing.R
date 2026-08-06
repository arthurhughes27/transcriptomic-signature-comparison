# =============================================================================
# Post-hoc DGSA result postprocessing
# =============================================================================
# Tidies a raw DGSA results list (one entry per vaccine x timepoint
# comparison, as produced by looping run_dearseq_comparison() or the
# equivalent QuSAGE function) into a single long-format results table, one
# row per gene set x comparison, and applies p-value adjustment across every
# combination of scope and method in Table 2.1 (18 columns: 3 scopes x 6
# methods). This is deliberately scoped to statistical postprocessing only -
# presentation concerns (condition/gene-set colours, factor ordering for
# plotting) are left to the plotting code that consumes this output.
# =============================================================================

#' Extract p-values, scores, correlations, and comparison labels from a raw
#' DGSA results list
#'
#' @param results_list Named list, one entry per comparison (name format
#'   "<vaccine> vs Control - Day <day>"), each entry a list with elements
#'   `pvals$rawPval`, `score$activation.scores`, `score$fc.scores`,
#'   `cor$mean.corr`, `cor$corr.mean`.
#'
#' @return A list with `raw_pvals`, `scores`, `correlations` (all named
#'   lists, one entry per comparison), and `conditions`/`times` (character/
#'   numeric vectors parsed from `names(results_list)`).
#'
#' @keywords internal
extract_dgsa_results <- function(results_list) {
  list(
    raw_pvals    = lapply(results_list, function(x) x[["pvals"]][["rawPval"]]),
    scores       = lapply(results_list, function(x)
      list(activation.scores = x[["score"]][["activation.scores"]],
           fc.scores         = x[["score"]][["fc.scores"]])),
    correlations = lapply(results_list, function(x)
      list(mean.corr = x[["cor"]][["mean.corr"]],
           corr.mean = x[["cor"]][["corr.mean"]])),
    conditions   = sub("^(.*?)\\s+vs.*", "\\1", names(results_list)),
    times        = as.numeric(sub(".*Day\\s+([0-9.]+)$", "\\1", names(results_list)))
  )
}

#' Apply one p-value adjustment method across all three scopes
#'
#' Adds three columns: `global.adjPval_{method}` (adjusted across every row
#' of `df`), `withinTime.adjPval_{method}` (adjusted within `time`), and
#' `withinComparison.adjPval_{method}` (adjusted within `comparison`, i.e.
#' within a single vaccine x timepoint).
#'
#' @keywords internal
apply_pvalue_adjustment_method <- function(df, method) {
  df |>
    dplyr::mutate(!!paste0("global.adjPval_", method) := stats::p.adjust(rawPval, method = method)) |>
    dplyr::group_by(time) |>
    dplyr::mutate(!!paste0("withinTime.adjPval_", method) := stats::p.adjust(rawPval, method = method)) |>
    dplyr::group_by(comparison) |>
    dplyr::mutate(!!paste0("withinComparison.adjPval_", method) := stats::p.adjust(rawPval, method = method)) |>
    dplyr::ungroup()
}

#' Apply p-value adjustment across all scopes and methods
#'
#' Extends a tidy one-row-per-(gene set, comparison) results table with
#' adjusted p-value columns for every combination of the three adjustment
#' scopes (global / within timepoint / within vaccine-timepoint) and the
#' given adjustment methods, following the naming convention
#' `{scope}.adjPval_{method}`.
#'
#' @param df Tidy results tibble with (at least) `rawPval`, `time`, and
#'   `comparison` columns (one row per gene set x comparison).
#' @param methods Character vector of `p.adjust()` methods to apply.
#'
#' @return `df` with `3 * length(methods)` additional columns.
apply_pvalue_adjustments <- function(df,
                                     methods = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY")) {
  Reduce(apply_pvalue_adjustment_method, methods, .init = df)
}

#' Tidy a raw DGSA results list into a long-format results table
#'
#' Converts the named list produced by looping `run_dearseq_comparison()`
#' (or the equivalent QuSAGE function) over every vaccine x timepoint
#' comparison into a single tidy tibble, one row per gene set x comparison,
#' and applies p-value adjustment across all three scopes x six methods
#' (see [apply_pvalue_adjustments()]). Entries in `results_list` that are
#' NULL (comparisons skipped for insufficient data) are dropped.
#'
#' @param results_list Named list, one entry per comparison, as described in
#'   [extract_dgsa_results()]. May contain NULL entries.
#' @param genesets Gene set list (as in data/BTM_processed.rds /
#'   data/BG3M_processed.rds), with `geneset.names`, `geneset.descriptions`,
#'   `geneset.aggregates`.
#' @param conditions_order Optional character vector giving the desired
#'   factor ordering of vaccine conditions (`condition` column). Conditions
#'   present in the data but not listed are appended, with a warning.
#' @param method Character label identifying the DGSA method (e.g.
#'   "dearseq", "qusage"), stored in the `method` column so results from
#'   multiple methods/specifications can be row-bound downstream.
#'
#' @return A tidy tibble, one row per gene set x comparison, with columns
#'   `comparison`, `condition`, `time`, `gs.name`, `gs.description`,
#'   `gs.name.description`, `gs.aggregate`, `activation.score`, `fc.score`,
#'   `mean.corr`, `corr.mean`, `rawPval`, `method`, and one
#'   `{scope}.adjPval_{method}` column per scope x adjustment-method
#'   combination.
build_tidy_dgsa_results <- function(results_list, genesets, conditions_order = NULL, method = NA_character_) {

  results_list <- Filter(Negate(is.null), results_list)

  if (length(results_list) == 0) {
    stop("results_list contains no completed comparisons (all entries are NULL).")
  }

  extracted    <- extract_dgsa_results(results_list)
  raw_pvals    <- extracted$raw_pvals
  scores       <- extracted$scores
  correlations <- extracted$correlations
  conditions   <- extracted$conditions
  times        <- extracted$times

  n_gene_sets <- length(genesets[["geneset.names"]])
  if (!all(vapply(raw_pvals, length, integer(1)) == n_gene_sets)) {
    stop("Number of genesets in `genesets` does not match the number of p-values in some comparisons.")
  }

  if (!is.null(conditions_order)) {
    unrecognised <- setdiff(conditions, conditions_order)
    if (length(unrecognised) > 0) {
      warning("The following conditions are not in conditions_order and will be appended: ",
              paste(unrecognised, collapse = ", "))
    }
    all_levels <- c(conditions_order, unrecognised)
    conditions <- factor(conditions, levels = intersect(all_levels, conditions))
  } else {
    conditions <- factor(conditions, levels = sort(unique(conditions)))
  }
  times <- factor(times, levels = sort(unique(times)))

  comparison_labels <- paste0(as.character(conditions), " - Day ", as.character(times))

  results_df <- purrr::map_dfr(seq_along(raw_pvals), function(idx) {
    tibble::tibble(
      comparison           = comparison_labels[idx],
      condition             = conditions[idx],
      time                   = times[idx],
      gs.name                 = genesets$geneset.names,
      gs.description            = genesets$geneset.descriptions,
      gs.name.description        = paste0(genesets$geneset.names, " - ", genesets$geneset.descriptions),
      gs.aggregate                 = genesets$geneset.aggregates,
      activation.score               = as.numeric(unlist(scores[[idx]][["activation.scores"]])),
      fc.score                        = as.numeric(unlist(scores[[idx]][["fc.scores"]])),
      mean.corr                        = correlations[[idx]][["mean.corr"]],
      corr.mean                         = correlations[[idx]][["corr.mean"]],
      rawPval                            = raw_pvals[[idx]],
      method                               = method
    )
  })

  apply_pvalue_adjustments(results_df)
}

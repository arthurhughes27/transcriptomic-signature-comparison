# =============================================================================
# Robustness vs. baseline significance (Figure 3)
# =============================================================================
# Relates the robustness metric pi_{g,v,j} (R/robustness_metrics.R) to each
# method's baseline-specification result (Chapter 2, Section 2.1.3 Stage
# 1's dearseq/QuSAGE baseline runs) for the same gene set x vaccine x
# timepoint comparison, to check whether comparisons with strong baseline
# evidence also tend to be robust across the specification grid. Two
# modes: "binary" (baseline significant/not, at a chosen adjusted-p-value
# column and alpha) and "continuous" (-log10 of the baseline p-value).
# Faceted by method, since the baseline specification and its result
# differ between dearseq and QuSAGE.
# =============================================================================

#' Join robustness metrics onto baseline DGSA results
#'
#' @param robustness_df Output of [compute_robustness_metric()] (`gs.name`,
#'   `condition`, `time`, `robustness`).
#' @param baseline_df One or more methods' baseline-specification results
#'   (e.g. `dplyr::bind_rows()` of `dearseq_dgsa_results_processed.rds` and
#'   `qusage_dgsa_results_processed.rds`), with `gs.name`, `condition`,
#'   `time`, `method`, `rawPval`, and `{scope}.adjPval_{method}` columns.
#'
#' @return An inner join of the two on `gs.name`/`condition`/`time`: one row
#'   per gene set x comparison x method, with `robustness` alongside every
#'   baseline column.
join_robustness_baseline <- function(robustness_df, baseline_df) {
  dplyr::inner_join(
    baseline_df, robustness_df,
    by = c("gs.name", "condition", "time")
  )
}

#' Plot robustness against baseline significance/evidence (Figure 3)
#'
#' @param df Output of [join_robustness_baseline()].
#' @param mode "binary" (x-axis is baseline significant TRUE/FALSE, a
#'   boxplot + jittered points per level) or "continuous" (x-axis is
#'   -log10 of the baseline p-value, one point per gene set x comparison,
#'   with a loess trend line).
#' @param pval_col Which baseline p-value column to use: an adjusted
#'   p-value column (see [apply_pvalue_adjustments()] for naming) for
#'   `mode = "binary"`'s significance call, or any p-value column (raw or
#'   adjusted) for `mode = "continuous"`'s x-axis. Defaults to the same
#'   `global.adjPval_BH` convention used for the Stage 1 comparison
#'   (R/comparison_metrics.R).
#' @param alpha Significance threshold applied to `pval_col`, `mode =
#'   "binary"` only.
#'
#' @return A ggplot object, faceted by `method`.
plot_robustness_vs_baseline <- function(df,
                                        mode     = c("binary", "continuous"),
                                        pval_col = "global.adjPval_BH",
                                        alpha    = 0.05) {
  mode <- match.arg(mode)
  stopifnot(
    pval_col %in% colnames(df),
    "robustness" %in% colnames(df),
    "method" %in% colnames(df)
  )

  plot_data <- dplyr::filter(df, !is.na(.data[[pval_col]]), !is.na(robustness))

  shared_layers <- list(
    ggplot2::facet_wrap(~method),
    ggplot2::coord_cartesian(ylim = c(0, 1)),
    ggplot2::theme_minimal(),
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background  = ggplot2::element_rect(fill = "grey90", colour = NA),
      strip.text         = ggplot2::element_text(face = "bold")
    )
  )

  if (mode == "binary") {
    plot_data <- dplyr::mutate(
      plot_data,
      baseline_significant = factor(
        ifelse(.data[[pval_col]] <= alpha, "Significant", "Not significant"),
        levels = c("Not significant", "Significant")
      )
    )

    ggplot2::ggplot(plot_data, ggplot2::aes(x = baseline_significant, y = robustness)) +
      ggplot2::geom_boxplot(fill = "grey90", outlier.shape = NA, width = 0.5) +
      ggplot2::geom_jitter(width = 0.15, alpha = 0.35, size = 0.8, colour = "#238b45") +
      shared_layers +
      ggplot2::labs(
        x     = sprintf("Baseline significant (%s <= %s)", pval_col, alpha),
        y     = "Robustness",
        title = "Robustness vs. baseline significance"
      )
  } else {
    plot_data <- dplyr::mutate(plot_data, neg_log10_p = -log10(pmax(.data[[pval_col]], 1e-300)))

    ggplot2::ggplot(plot_data, ggplot2::aes(x = neg_log10_p, y = robustness)) +
      ggplot2::geom_point(alpha = 0.35, size = 0.8, colour = "grey40") +
      ggplot2::geom_smooth(method = "loess", formula = y ~ x, colour = "#238b45", fill = "#238b45", alpha = 0.15) +
      shared_layers +
      ggplot2::labs(
        x     = sprintf("-log10(baseline %s)", pval_col),
        y     = "Robustness",
        title = "Robustness vs. baseline evidence"
      )
  }
}

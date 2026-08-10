# =============================================================================
# DGSA heatmap comparison: dearseq vs QuSAGE
# =============================================================================
# ComplexHeatmap comparison of the two methods' baseline results, analogous
# to the circos plots (analysis/reanalysis/plot_circos.R): mean fold-change
# score per gene set x vaccine, split into column panels by timepoint,
# asterisk-annotated where significant, with a row strip coloured by
# gene-set aggregate.
#
# Cleaned up and modularised from a script carried over from an earlier
# repository (analysis/reanalysis/plot_heatmap_comparison.R). The original
# monolithic plot_heatmap() function is split into the smaller helpers
# below, each handling one step of the original's five-step pipeline:
#   1. determine_heatmap_significance() - flag rows significant
#   2. select_heatmap_comparisons()     - restrict to chosen conditions/
#                                          times/aggregates
#   3. filter_common_de()               - optional "commonly
#                                          differentially-expressed" gene
#                                          set filtering
#   4. clip_scores()                    - clip extreme scores for colour
#                                          scaling
#   5. build_dgsa_heatmap()             - assemble the ComplexHeatmap
# plot_dgsa_heatmap() orchestrates all five, preserving the original
# function's external contract. filter_common_de() no longer closes over a
# package-level `results_df` (the original relied on a global variable for
# its "score" mode) - it now takes the full results tibble as an explicit
# argument.
# =============================================================================

#' Abbreviate long gene-set description text for heatmap row labels
#'
#' @param x Character vector of gene-set descriptions.
#'
#' @return `x` with common long phrases replaced by shorter abbreviations.
abbreviate_geneset_label <- function(x) {
  x <- gsub("regulation", "reg.", x)
  x <- gsub("activation", "act.", x)
  x <- gsub("signaling", "sign.", x)
  x <- gsub("transcription", "transcr.", x)
  x <- gsub("differentiation", "diff.", x)
  x <- gsub("extracellular matrix", "ECM", x)
  x <- gsub("mitotic cell cycle", "mit. cell cycle", x)
  x <- gsub("cell division", "cell div.", x)
  x <- gsub("antigen presentation", "antigen pres.", x)
  x <- gsub("immune response", "imm. resp.", x)
  x <- gsub("transporters", "transp.", x)
  x <- gsub("T cells", "T-cells", x)
  x <- gsub("B cells", "B-cells", x)
  x <- gsub("plasma membrane", "PM", x)
  x <- gsub("E2F transcription factor network", "E2F TF network", x)
  x <- gsub("C-MYC transcriptional network", "C-MYC network", x)
  x <- gsub("Rho GTPase cycle", "Rho cycle", x)
  x <- gsub("lipid metabolism", "lipid metab.", x)
  x <- gsub("chemokines and receptors", "chemokines/receptors.", x)
  x <- gsub("adhesion and migration", "adhesion/migration", x)
  x <- gsub("growth factor induced", "GF induced", x)
  x <- gsub("viral sensing & immunity", "viral sens./imm.", x)
  x <- gsub("network", "netw.", x)
  x <- gsub("presentation", "pres.", x)
  x <- gsub("complement & other receptors", "complement/other receptors", x)
  x
}

#' Flag rows as significant for heatmap cell annotation
#'
#' @param df Tidy DGSA results (one row per gene set x comparison), with
#'   `condition`, `time`, `score_col`, and the `{p_approach}.adjPval_
#'   {p_correction}` column.
#' @param p_approach,p_correction,p_threshold Adjusted p-value column
#'   (see [apply_pvalue_adjustments()] for naming) and threshold.
#' @param filter_mode "none" (p-value only), "user" (p-value AND absolute
#'   score above `user_threshold`), or "data" (p-value AND absolute score
#'   above a `quantile_threshold` computed from the data - globally, within
#'   `time`, or within `condition` x `time`, depending on `p_approach`).
#' @param score_col Score column used for the effect-size condition.
#' @param user_threshold,quantile_threshold See `filter_mode`.
#'
#' @return `df` with a `significant` factor column added (levels
#'   `c(TRUE, FALSE)`).
determine_heatmap_significance <- function(df,
                                           p_approach,
                                           p_correction,
                                           p_threshold,
                                           filter_mode         = c("none", "user", "data"),
                                           score_col            = "fc.score",
                                           user_threshold        = 0.5,
                                           quantile_threshold     = 0.5) {
  filter_mode <- match.arg(filter_mode)

  adj_pval_col <- paste0(p_approach, ".adjPval_", p_correction)
  stopifnot(adj_pval_col %in% colnames(df), score_col %in% colnames(df))

  base_significant <- df[[adj_pval_col]] < p_threshold

  effect_pass <- switch(
    filter_mode,
    none = TRUE,
    user = abs(df[[score_col]]) > user_threshold,
    data = if (p_approach == "global") {
      threshold <- stats::quantile(abs(df[[score_col]]), quantile_threshold, na.rm = TRUE)
      abs(df[[score_col]]) > threshold
    } else if (p_approach == "withinTime") {
      thresholds <- df |>
        dplyr::group_by(time) |>
        dplyr::summarise(
          threshold_abs = stats::quantile(abs(.data[[score_col]]), quantile_threshold, na.rm = TRUE),
          .groups = "drop"
        )
      df |>
        dplyr::left_join(thresholds, by = "time") |>
        dplyr::mutate(pass = abs(.data[[score_col]]) > threshold_abs) |>
        dplyr::pull(pass)
    } else {
      df |>
        dplyr::group_by(condition, time) |>
        dplyr::mutate(pass = abs(.data[[score_col]]) > stats::quantile(abs(.data[[score_col]]), quantile_threshold, na.rm = TRUE)) |>
        dplyr::ungroup() |>
        dplyr::pull(pass)
    }
  )

  dplyr::mutate(df, significant = factor(base_significant & effect_pass, levels = c(TRUE, FALSE)))
}

#' Restrict to the chosen conditions/times/aggregates for one heatmap
#'
#' @param df Tidy DGSA results with `gs.aggregate`, `condition`, `time`.
#' @param conditions,times,aggregates Values to keep.
#'
#' @return The filtered subset of `df`, with a `time_label` column added.
select_heatmap_comparisons <- function(df, conditions, times, aggregates) {
  df |>
    dplyr::filter(gs.aggregate %in% aggregates, condition %in% conditions, time %in% times) |>
    dplyr::mutate(
      time_label = factor(
        paste0(condition, " - Day ", time),
        levels = unique(paste0(condition, " - Day ", time))
      )
    )
}

#' Optionally restrict to "commonly differentially-expressed" gene sets
#'
#' @param df Already-significance-flagged, comparison-selected results
#'   (output of [select_heatmap_comparisons()] after
#'   [determine_heatmap_significance()]).
#' @param full_df The FULL baseline results (both methods, every
#'   comparison) - needed by `filter_commonDE = "score"`, which computes
#'   each gene set's per-method "sharing score" (the number of vaccines it
#'   is significant, in its predominant direction, for) from the complete
#'   data rather than just the rows selected for this particular heatmap.
#' @param method_name Method label used to restrict `full_df` for the
#'   `"score"` mode.
#' @param filter_commonDE "none" (no filtering), "global" (keep gene sets
#'   significant in >= `common_proportion` of the selected rows), "within
#'   Time" (keep gene sets significant in >= `common_proportion` of rows
#'   for at least one selected time), or "score" (keep gene sets with a
#'   sharing score >= `score_threshold`).
#' @param common_proportion,score_threshold See `filter_commonDE`.
#' @param adj_pval_col,score_col,p_threshold Passed through to the sharing
#'   score computation (`filter_commonDE = "score"` only).
#'
#' @return The filtered subset of `df`.
filter_common_de <- function(df, full_df, method_name,
                             filter_commonDE    = c("none", "global", "withinTime", "score"),
                             common_proportion    = 0,
                             score_threshold        = 0,
                             adj_pval_col            = NULL,
                             score_col                = "fc.score",
                             p_threshold                = 0.05) {
  filter_commonDE <- match.arg(filter_commonDE)

  if (filter_commonDE == "none") return(df)

  if (filter_commonDE == "global") {
    keep_sets <- df |>
      dplyr::group_by(gs.name) |>
      dplyr::summarise(frac_sig = sum(as.logical(significant)) / dplyr::n(), .groups = "drop") |>
      dplyr::filter(frac_sig >= common_proportion) |>
      dplyr::pull(gs.name)
    return(dplyr::filter(df, gs.name %in% keep_sets))
  }

  if (filter_commonDE == "withinTime") {
    keep_sets <- df |>
      dplyr::group_by(gs.name, time) |>
      dplyr::summarise(frac_sig = sum(significant == "TRUE") / dplyr::n(), .groups = "drop") |>
      dplyr::filter(frac_sig >= common_proportion) |>
      dplyr::pull(gs.name) |>
      unique()
    return(dplyr::filter(df, gs.name %in% keep_sets))
  }

  # filter_commonDE == "score": sharing score - for each gene set, the
  # number of unique vaccines (among those selected for this heatmap) for
  # which it is significant in its predominant direction, computed from
  # the full per-method results rather than just the selected rows.
  stopifnot(!is.null(adj_pval_col))

  sharing_score_df <- full_df |>
    dplyr::filter(condition %in% unique(df$condition), method == method_name)

  direction_counts <- sharing_score_df |>
    dplyr::mutate(
      direction      = ifelse(.data[[score_col]] > 0, "up", "down"),
      is_significant = !is.na(.data[[adj_pval_col]]) & .data[[adj_pval_col]] < p_threshold
    ) |>
    dplyr::select(condition, gs.name, direction, is_significant) |>
    dplyr::distinct() |>
    dplyr::group_by(gs.name, direction) |>
    dplyr::summarise(n = sum(is_significant), .groups = "drop") |>
    dplyr::group_by(gs.name) |>
    dplyr::slice_max(n, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  keep_sets <- dplyr::filter(direction_counts, n >= score_threshold)$gs.name

  dplyr::filter(df, gs.name %in% keep_sets)
}

#' Clip a score column at a quantile of its absolute value
#'
#' @param df Tidy results tibble.
#' @param score_col Score column to clip.
#' @param quantile_scoreclip Quantile of `abs(score_col)` to clip at (1 =
#'   no clipping).
#'
#' @return `df` with `score_col` clipped to `[-threshold, threshold]`.
clip_scores <- function(df, score_col, quantile_scoreclip) {
  clip_threshold <- stats::quantile(abs(df[[score_col]]), quantile_scoreclip, na.rm = TRUE)
  dplyr::mutate(df, !!score_col := pmin(pmax(.data[[score_col]], -clip_threshold), clip_threshold))
}

#' Assemble the ComplexHeatmap object for one method's heatmap
#'
#' @param df Cleaned, filtered, clipped results for this heatmap (output
#'   of the earlier pipeline steps).
#' @param full_results_df The full baseline results (both methods) -
#'   needed for the gene-set-aggregate factor levels/colours, which should
#'   be consistent across every heatmap regardless of which aggregates
#'   happen to be present in `df`.
#' @param score_col Score column to display as heatmap fill.
#' @param y_order "cluster" (hierarchical clustering) or "aggregate"
#'   (grouped by gene-set aggregate, alphabetical within).
#' @param x_order "cluster" or "set" (as given).
#' @param fixed_row_names_width `grid::unit()` for the row-name column
#'   width.
#' @param legend_max Fixed colour-scale limit (NULL to use the data's own
#'   max absolute value).
#'
#' @return A `ComplexHeatmap::Heatmap` object.
build_dgsa_heatmap <- function(df, full_results_df, score_col,
                               y_order                = c("cluster", "aggregate"),
                               x_order                = c("cluster", "set"),
                               fixed_row_names_width  = grid::unit(150, "mm"),
                               legend_max              = NULL) {
  y_order <- match.arg(y_order)
  x_order <- match.arg(x_order)

  df         <- dplyr::mutate(df, time = droplevels(time))
  time_order <- levels(df$time)
  sep        <- "___"

  df_wide <- df |>
    dplyr::mutate(col_id = paste0(condition, sep, time)) |>
    dplyr::select(gs.name.description, col_id, !!rlang::sym(score_col)) |>
    tidyr::pivot_wider(names_from = col_id, values_from = !!rlang::sym(score_col))

  mat <- as.matrix(dplyr::select(df_wide, -gs.name.description))
  rownames(mat) <- df_wide$gs.name.description
  col_full      <- colnames(mat)
  col_time_raw  <- sub(".*___", "", col_full)
  col_condition <- sub("___.*", "", col_full)

  time_to_daylabel <- function(t) {
    num <- stringr::str_extract(as.character(t), "\\d+")
    if (is.na(num)) paste0("Day ", as.character(t)) else paste0("Day ", num)
  }
  time_display_levels <- vapply(time_order, time_to_daylabel, FUN.VALUE = character(1))
  col_time_display <- factor(
    vapply(col_time_raw, time_to_daylabel, FUN.VALUE = character(1)),
    levels = time_display_levels
  )

  # Light, day-consistent colours behind each column-slice title (same
  # palette as the ggplot2 figures' facet strips - see
  # R/plot_helpers.R's assign_day_colours()), so it's clear at a glance
  # where one timepoint's block of columns ends and the next begins.
  day_title_colours <- assign_day_colours(time_order)[time_display_levels]

  colnames(mat) <- col_condition

  all_aggregates <- levels(full_results_df$gs.aggregate)
  rows_meta <- df |>
    dplyr::distinct(gs.name.description, gs.aggregate, gs.colour) |>
    dplyr::arrange(match(gs.name.description, rownames(mat)))
  row_groups_vec <- factor(rows_meta$gs.aggregate, levels = all_aggregates)
  names(row_groups_vec) <- rows_meta$gs.name.description

  row_group_cols <- full_results_df |>
    dplyr::distinct(gs.aggregate, gs.colour) |>
    dplyr::filter(gs.aggregate %in% all_aggregates)
  row_group_cols <- stats::setNames(row_group_cols$gs.colour, row_group_cols$gs.aggregate)

  cond_col_map <- df |>
    dplyr::distinct(condition, condition.colour) |>
    tibble::deframe()

  sig_wide <- df |>
    dplyr::mutate(col_id = paste0(condition, sep, time)) |>
    dplyr::select(gs.name.description, col_id, significant) |>
    tidyr::pivot_wider(names_from = col_id, values_from = significant)

  sig_mat <- as.matrix(dplyr::select(sig_wide, -gs.name.description))
  rownames(sig_mat) <- sig_wide$gs.name.description
  sig_mat <- matrix(
    as.logical(sig_mat),
    nrow = nrow(sig_mat), ncol = ncol(sig_mat), dimnames = dimnames(sig_mat)
  )
  sig_mat <- sig_mat[rownames(mat), col_full, drop = FALSE]
  colnames(sig_mat) <- col_condition
  sig_mat[is.na(sig_mat)] <- FALSE

  max_abs <- if (is.null(legend_max)) max(abs(mat[!is.na(mat)]), na.rm = TRUE) else legend_max
  col_fun <- circlize::colorRamp2(c(-max_abs, 0, max_abs), c("blue", "white", "red"))

  cluster_rows <- y_order == "cluster"
  cluster_cols <- x_order == "cluster"

  if (!cluster_rows) {
    ord            <- order(match(row_groups_vec, levels(row_groups_vec)))
    mat            <- mat[ord, , drop = FALSE]
    row_groups_vec <- row_groups_vec[ord]
    sig_mat        <- sig_mat[ord, , drop = FALSE]
  }
  show_row_dend <- cluster_rows

  row_ha <- ComplexHeatmap::rowAnnotation(
    `Geneset aggregate` = row_groups_vec,
    col = list(`Geneset aggregate` = row_group_cols),
    show_annotation_name = FALSE,
    width = grid::unit(6, "mm")
  )

  top_ha <- ComplexHeatmap::HeatmapAnnotation(
    Vaccine = col_condition,
    col = list(Vaccine = cond_col_map),
    show_annotation_name = FALSE,
    annotation_height = grid::unit(2, "mm")
  )

  ComplexHeatmap::Heatmap(
    mat,
    name                  = "Mean fold-change",
    col                   = col_fun,
    column_gap            = grid::unit(8, "mm"),
    na_col                = "grey95",
    cluster_rows          = cluster_rows,
    show_row_dend         = show_row_dend,
    row_dend_side         = "left",
    cluster_columns       = cluster_cols,
    column_split          = col_time_display,
    cluster_column_slices = FALSE,
    column_title_gp       = grid::gpar(fill = unname(day_title_colours), fontface = "bold", fontsize = 12),
    show_column_dend      = TRUE,
    column_dend_side      = "top",
    column_dend_height    = grid::unit(2, "cm"),
    row_names_gp          = grid::gpar(fontsize = 10),
    column_names_side     = "bottom",
    top_annotation        = top_ha,
    left_annotation       = row_ha,
    show_column_names     = TRUE,
    show_row_names        = TRUE,
    row_dend_width        = grid::unit(4, "cm"),
    heatmap_legend_param  = list(title = "Mean fold-change"),
    row_names_max_width   = fixed_row_names_width,
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (isTRUE(sig_mat[i, j])) {
        grid::grid.text("*", x = x, y = y, gp = grid::gpar(fontsize = 10))
      }
    }
  )
}

#' Build one method's DGSA comparison heatmap
#'
#' Orchestrates significance determination, comparison/aggregate
#' selection, optional "commonly differentially-expressed" gene set
#' filtering, score clipping, and ComplexHeatmap construction (see the
#' module-level comment above for the breakdown into steps).
#'
#' @param results_df Baseline results for BOTH methods (e.g.
#'   [load_baseline_results_from_raw()]'s output, with
#'   `assign_dgsa_colours()` applied for `condition.colour`/`gs.colour`).
#'   The "score" `filter_commonDE` mode and the aggregate colour/level
#'   lookup draw on the full dataset, not just the rows selected for this
#'   particular heatmap.
#' @param method_name "dearseq" or "qusage" (NULL to not filter by
#'   method).
#' @param conditions,times,aggregates Values to include in this heatmap.
#' @param fixed_row_names_width `grid::unit()` for the row-name column
#'   width.
#' @param p_correction,p_approach,p_threshold Baseline significance
#'   definition (see [determine_heatmap_significance()]).
#' @param filter_mode,scores,user_threshold,quantile_threshold Effect-size
#'   filtering for the significance call (see
#'   [determine_heatmap_significance()]).
#' @param y_order,x_order Row/column ordering (see [build_dgsa_heatmap()]).
#' @param filter_commonDE,common_proportion,score_threshold "Commonly
#'   differentially-expressed" gene set filtering (see
#'   [filter_common_de()]).
#' @param quantile_scoreclip Score-clipping quantile (see [clip_scores()]).
#' @param legend_max Fixed colour-scale limit (NULL to use the data's own
#'   max absolute value).
#'
#' @return A `ComplexHeatmap::Heatmap` object.
plot_dgsa_heatmap <- function(results_df,
                              method_name          = NULL,
                              conditions,
                              times,
                              aggregates,
                              fixed_row_names_width = grid::unit(150, "mm"),
                              p_correction          = c("BH", "bonferroni", "holm", "hommel", "hochberg", "BY"),
                              p_approach             = c("global", "withinTime", "withinComparison"),
                              p_threshold             = 0.05,
                              filter_mode              = c("none", "user", "data"),
                              user_threshold            = 0.5,
                              quantile_threshold         = 0.5,
                              scores                      = c("fc.score", "activation.score"),
                              y_order                      = c("cluster", "aggregate"),
                              x_order                       = c("cluster", "set"),
                              filter_commonDE                = c("none", "global", "withinTime", "score"),
                              common_proportion                = 0,
                              score_threshold                    = 0,
                              quantile_scoreclip                  = 0.995,
                              legend_max                            = NULL) {
  p_correction    <- match.arg(p_correction)
  p_approach      <- match.arg(p_approach)
  filter_mode     <- match.arg(filter_mode)
  scores          <- match.arg(scores)
  y_order         <- match.arg(y_order)
  x_order         <- match.arg(x_order)
  filter_commonDE <- match.arg(filter_commonDE)

  df <- results_df
  if (!is.null(method_name)) df <- dplyr::filter(df, method == method_name)

  adj_pval_col <- paste0(p_approach, ".adjPval_", p_correction)

  df <- determine_heatmap_significance(
    df,
    p_approach = p_approach, p_correction = p_correction, p_threshold = p_threshold,
    filter_mode = filter_mode, score_col = scores,
    user_threshold = user_threshold, quantile_threshold = quantile_threshold
  )

  df <- select_heatmap_comparisons(df, conditions, times, aggregates)

  df <- filter_common_de(
    df, full_df = results_df, method_name = method_name,
    filter_commonDE = filter_commonDE, common_proportion = common_proportion,
    score_threshold = score_threshold, adj_pval_col = adj_pval_col,
    score_col = scores, p_threshold = p_threshold
  )

  if (nrow(df) == 0) stop("No gene sets to plot under current parameters!")

  df <- clip_scores(df, score_col = scores, quantile_scoreclip = quantile_scoreclip)

  build_dgsa_heatmap(
    df, full_results_df = results_df, score_col = scores,
    y_order = y_order, x_order = x_order,
    fixed_row_names_width = fixed_row_names_width, legend_max = legend_max
  )
}

#' Save one or more heatmaps stacked vertically into a single PDF
#'
#' Consolidates the repeated `pushViewport()`/`draw()`/`grid.text()` block
#' pattern needed for every multi-heatmap PDF in the driver script into a
#' single loop over N heatmaps.
#'
#' @param path Output PDF path.
#' @param heatmaps List of `ComplexHeatmap::Heatmap` objects.
#' @param titles Character vector, same length as `heatmaps` - a bold
#'   title drawn above each heatmap panel.
#' @param heights List of `grid::unit()` values, same length as
#'   `heatmaps` - the relative height of each heatmap's row.
#' @param width,height PDF device size (inches).
#' @param spacer_height `grid::unit()` for the blank row between stacked
#'   heatmaps.
#' @param title_size,title_offset Title font size and vertical offset
#'   above each heatmap panel.
#'
#' @return `path`, invisibly.
save_stacked_heatmaps_pdf <- function(path, heatmaps, titles, heights,
                                      width, height,
                                      spacer_height = grid::unit(30, "mm"),
                                      title_size     = 32,
                                      title_offset    = grid::unit(10, "mm")) {
  stopifnot(length(heatmaps) == length(titles), length(heatmaps) == length(heights))

  n_panels <- length(heatmaps)
  n_rows   <- 2 * n_panels - 1

  row_heights <- vector("list", n_rows)
  for (i in seq_len(n_panels)) {
    row_heights[[2 * i - 1]] <- heights[[i]]
    if (i < n_panels) row_heights[[2 * i]] <- spacer_height
  }
  row_heights <- do.call(grid::unit.c, row_heights)

  grDevices::pdf(path, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)

  grid::pushViewport(grid::viewport(
    layout = grid::grid.layout(nrow = n_rows, ncol = 1, heights = row_heights)
  ))

  for (i in seq_len(n_panels)) {
    panel_row <- 2 * i - 1
    grid::pushViewport(grid::viewport(layout.pos.row = panel_row, layout.pos.col = 1))
    ComplexHeatmap::draw(
      heatmaps[[i]],
      heatmap_legend_side = "right", annotation_legend_side = "right", newpage = FALSE
    )
    grid::grid.text(
      titles[i],
      x = grid::unit(0.5, "npc"), y = grid::unit(1, "npc") + title_offset,
      gp = grid::gpar(fontsize = title_size, fontface = "bold")
    )
    grid::upViewport()
  }

  invisible(path)
}

# =============================================================================
# Shared DGSA helpers
# =============================================================================
# Generic building blocks used by both the dearseq (R/dgsa_dearseq.R) and
# QuSAGE (R/dgsa_qusage.R, forthcoming) differential gene-set analysis
# pipelines: listing valid vaccine x timepoint comparisons, pairing each
# participant's pre-/post-vaccination samples, building covariate design
# matrices, and computing gene-set-level correlations with immune response.
# =============================================================================

#' List every vaccine x timepoint comparison with usable data
#'
#' A comparison is included if at least one participant receiving that
#' vaccine has both a pre-vaccination sample (time_post_last_vax <= 0) and
#' the post-vaccination sample at that timepoint.
#'
#' @param hipc Merged clinical/expression tibble (as in
#'   data/hipc_merged_*.rds), with columns `vaccine_name`, `participant_id`,
#'   and `time_post_last_vax`.
#' @param days Optional numeric vector restricting the result to specific
#'   post-vaccination timepoints (days). NULL (default) includes every
#'   post-vaccination timepoint present in `hipc`. Days with no matching
#'   data in `hipc` are silently ignored (rather than erroring), so this can
#'   be set to a fixed "days of interest" vector shared across datasets that
#'   don't all cover every day.
#'
#' @return A tibble with one row per valid comparison, columns
#'   `vaccine_name` and `day`.
list_valid_comparisons <- function(hipc, days = NULL) {
  timepoints <- hipc |>
    dplyr::filter(time_post_last_vax > 0) |>
    dplyr::pull(time_post_last_vax) |>
    unique() |>
    sort()

  if (!is.null(days)) {
    timepoints <- intersect(timepoints, days)
  }

  purrr::map_dfr(timepoints, function(day) {
    valid_vaccines <- hipc |>
      dplyr::group_by(vaccine_name, participant_id) |>
      dplyr::summarise(
        has_post = any(time_post_last_vax == day),
        has_pre  = any(time_post_last_vax <= 0),
        .groups  = "drop"
      ) |>
      dplyr::filter(has_post & has_pre) |>
      dplyr::pull(vaccine_name) |>
      unique()

    if (length(valid_vaccines) == 0) return(NULL)

    tibble::tibble(vaccine_name = valid_vaccines, day = day)
  })
}

#' Pair each participant's pre- and post-vaccination samples
#'
#' Filters `hipc` to participants of vaccine `vax` who have both a
#' pre-vaccination sample and the post-vaccination sample at `day`. Where a
#' participant has multiple pre-vaccination samples, only the most recent
#' (by `study_time_collected`) is kept.
#'
#' @param hipc Merged clinical/expression tibble.
#' @param vax Vaccine name to filter to (matches `hipc$vaccine_name`).
#' @param day Post-vaccination timepoint (matches `hipc$time_post_last_vax`).
#'
#' @return A tibble with one row per retained sample, `participant_id`
#'   recoded to a compact factor.
filter_paired_samples <- function(hipc, vax, day) {
  hipc |>
    dplyr::filter(
      vaccine_name == vax,
      time_post_last_vax == day | time_post_last_vax <= 0
    ) |>
    dplyr::group_by(participant_id) |>
    dplyr::filter(
      any(time_post_last_vax == day),
      any(time_post_last_vax <= 0)
    ) |>
    dplyr::mutate(
      max_pre_time = max(study_time_collected[time_post_last_vax <= 0], na.rm = TRUE)
    ) |>
    dplyr::filter(time_post_last_vax == day | study_time_collected == max_pre_time) |>
    dplyr::select(-max_pre_time) |>
    dplyr::ungroup() |>
    dplyr::mutate(participant_id = as.factor(as.numeric(as.factor(participant_id))))
}

#' Build a dearseq covariate design matrix
#'
#' Builds a design matrix adjusting for `vaccine_code` (the pre/post
#' indicator) plus an arbitrary subset of baseline covariates, removes any
#' linearly dependent columns via QR decomposition, and drops the
#' `vaccine_code` columns themselves, leaving only the adjustment covariates
#' (as required by `dearseq::dgsa_seq()`, which tests `vaccine_code`
#' separately via `variables2test`).
#'
#' @param df Sample-level tibble (one comparison's worth of paired samples,
#'   as returned by [filter_paired_samples()]), with a `vaccine_code` column.
#' @param covariates Character vector of covariate column names to adjust
#'   for (may be empty, for "no adjustment").
#'
#' @return A numeric matrix of adjustment covariates (intercept + any
#'   independent covariate columns).
build_covariate_matrix <- function(df, covariates = c("age_imputed", "gender", "study_accession", "race")) {
  rhs    <- paste(c("vaccine_code", covariates), collapse = " + ")
  X_full <- stats::model.matrix(stats::as.formula(paste("~", rhs)), data = df)

  qrX     <- qr(X_full)
  X_indep <- X_full[, sort(qrX$pivot[seq_len(qrX$rank)]), drop = FALSE]
  x       <- X_indep[, !grepl("^vaccine_code", colnames(X_indep)), drop = FALSE]

  message("Adjusting for covariates: ", paste(colnames(x), collapse = ", "))
  x
}

#' Gene-set-level correlation with an immune response variable
#'
#' For each gene set, computes two summaries of its correlation with
#' `response` across samples: the mean of per-gene Spearman correlations
#' (`mean.corr`), and the Spearman correlation of a standardised
#' per-sample gene-set meta-score (`corr.mean`).
#'
#' @param y Gene x sample expression matrix (genes in rownames).
#' @param response Numeric vector, one value per sample (column of `y`).
#' @param GSA Gene set list (as in data/BTM_processed.rds), with
#'   `genesets` and `geneset.names`.
#'
#' @return A list with named numeric vectors `mean.corr` and `corr.mean`,
#'   one value per gene set.
calculate_gs_correlation <- function(y, response, GSA) {
  response <- as.matrix(response)

  n <- ncol(y)
  if (length(response) != n) {
    stop("Number of response values given does not match number of samples in y")
  }

  mean.corr <- sapply(seq_along(GSA[["genesets"]]), function(i) {
    gs   <- GSA[["genesets"]][[i]]
    keep <- intersect(gs, rownames(y))
    if (length(keep) == 0) return(NA_real_)
    cors <- apply(y[keep, , drop = FALSE], 1, function(gene_expr)
      stats::cor(gene_expr, response, method = "spearman", use = "pairwise.complete.obs"))
    mean(cors, na.rm = TRUE)
  })
  names(mean.corr) <- GSA[["geneset.names"]]

  corr.mean <- sapply(seq_along(GSA[["genesets"]]), function(i) {
    gs   <- GSA[["genesets"]][[i]]
    keep <- intersect(gs, rownames(y))
    if (length(keep) == 0) return(NA_real_)

    y_sub <- y[keep, , drop = FALSE]
    zsub  <- t(scale(t(y_sub), center = TRUE, scale = TRUE))
    meta  <- colMeans(zsub, na.rm = TRUE)

    stats::cor(meta, response, method = "spearman", use = "pairwise.complete.obs")
  })
  names(corr.mean) <- GSA[["geneset.names"]]

  list(mean.corr = mean.corr, corr.mean = corr.mean)
}

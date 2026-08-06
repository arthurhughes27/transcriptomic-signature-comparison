# =============================================================================
# Modular QuSAGE differential gene-set analysis
# =============================================================================
# run_qusage_comparison() runs one vaccine x timepoint comparison with
# qusage::qusage() + qusage::combinePDFs() (Yaari et al. 2013; Meng et al.
# 2019 multi-study meta-analysis), replacing the monolithic, Hagan-et-al.
# derived analysis/reanalysis/qusage_dgsa.R with a function matching the
# same contract as run_dearseq_comparison() (R/dgsa_dearseq.R), so both
# methods can be driven by the same comparison loop and specification grid.
#
# Design notes / deliberate deviations from the original script:
#
#   - Sample pairing uses the shared filter_paired_samples() helper
#     (R/dgsa_common.R) - the same one dearseq uses - rather than the
#     original script's separate keep_most_recent_prevax_time_point()
#     applied globally across the whole dataset. This guarantees both
#     methods analyse literally the same set of samples for a given
#     comparison, as intended by Chapter 2, Section 2.2.2.
#   - Genes with any missing value are dropped per comparison (across all
#     studies contributing to that vaccine x timepoint), not once globally
#     across the entire dataset. This still guarantees every per-study
#     QSarray tests the same gene sets (required for combinePDFs()), while
#     not discarding genes for comparisons unaffected by NAs elsewhere.
#   - Correlation with immune response reuses calculate_gs_correlation()
#     (R/dgsa_common.R), the same helper dearseq uses, rather than the
#     original script's separate weighted multi-study correlation pipeline
#     (which fed the circos plots' correlation ring, not the DGSA test
#     itself). This is not part of the QuSAGE algorithm's own output, so it
#     is not covered by the reproducibility check against the old script.
# =============================================================================

#' Keep only specified symbols in a given list of gene sets
#'
#' Restricts every gene set to genes present in `symbols2keep`, then drops
#' any gene set left with fewer than `min_size` genes. (Unchanged from the
#' original Hagan-et-al.-derived implementation.)
#'
#' @param genesets.symbol Named list of character vectors (gene symbols).
#' @param symbols2keep Character vector of gene symbols to retain.
#' @param min_size Minimum number of genes required to keep a gene set.
#'
#' @return The filtered, still-named list.
#' @keywords internal
keep_symbols <- function(genesets.symbol, symbols2keep, min_size = 2) {
  stopifnot(is.character(symbols2keep))
  stopifnot(is.list(genesets.symbol))
  common_symbols <- intersect(Reduce(union, genesets.symbol), symbols2keep)
  genesets.symbol |>
    purrr::map(~ .x[.x %in% common_symbols]) |>
    purrr::discard(~ length(.x) < min_size)
}

#' Run one QuSAGE comparison (vaccine x timepoint)
#'
#' Filters to paired pre-/post-vaccination samples for `vax` at `day`
#' ([filter_paired_samples()]), restricts to genes complete across those
#' samples and gene sets retaining at least 2 such genes
#' ([keep_symbols()]), runs `qusage::qusage()` independently within each
#' contributing study, and meta-analyses the resulting per-study QSarray
#' objects across studies via `qusage::combinePDFs()` (a single study's
#' QSarray is used directly with no meta-analysis needed). Also computes
#' gene-set correlations with immune response ([calculate_gs_correlation()]).
#'
#' @param vax Vaccine name.
#' @param day Post-vaccination timepoint (days).
#' @param hipc Merged clinical/expression tibble.
#' @param BTM Gene set list (genesets + geneset.names).
#' @param gene_names Character vector of gene expression column names in
#'   `hipc`.
#' @param equal_variance QuSAGE `var.equal` argument: whether to assume
#'   equal variance between the pre- and post-vaccination groups. FALSE
#'   (default) is the baseline specification (Table 2.1).
#'
#' @return A list with elements `pvals` (rawPval, one value per BTM gene
#'   set, NA for gene sets dropped by the completeness filter), `score`
#'   (activation.scores/fc.scores - QuSAGE has no separate fold-change
#'   score, so both equal the QSarray `path.mean`), and `cor` (gene-set
#'   correlations with immune response); or NULL if there is insufficient
#'   data for this comparison.
run_qusage_comparison <- function(vax, day, hipc, BTM, gene_names, equal_variance = FALSE) {

  df <- filter_paired_samples(hipc, vax, day)

  if (nrow(df) == 0) {
    message("  No data after filtering — skipping.")
    return(NULL)
  }

  # Genes complete across every sample in this comparison, so every study's
  # QSarray tests an identical set of gene sets and combinePDFs() is valid.
  exprmat_all    <- t(as.matrix(dplyr::select(df, dplyr::all_of(gene_names))))
  complete_genes <- rownames(exprmat_all)[stats::complete.cases(exprmat_all)]

  genesets <- BTM[["genesets"]]
  names(genesets) <- BTM[["geneset.names"]]
  genesets <- keep_symbols(genesets, complete_genes, min_size = 2)

  if (length(genesets) == 0) {
    message("  No gene sets with >= 2 complete genes — skipping.")
    return(NULL)
  }

  studies <- sort(unique(as.character(df$study_accession)))

  qsarray_list <- list()
  for (study in studies) {
    df_study <- df[as.character(df$study_accession) == study, ]

    labels <- ifelse(df_study$time_post_last_vax == day, "post", "pre")
    if (length(unique(labels)) < 2 || dplyr::n_distinct(df_study$participant_id) < 2) next

    exprmat <- t(as.matrix(dplyr::select(df_study, dplyr::all_of(complete_genes))))
    colnames(exprmat) <- as.character(df_study$participant_id)

    result <- tryCatch(
      qusage::qusage(
        eset        = exprmat,
        labels      = labels,
        contrast    = "post-pre",
        geneSets    = genesets,
        pairVector  = as.character(df_study$participant_id),
        var.equal   = equal_variance
      ),
      error = function(e) {
        message("    QuSAGE failed for study ", study, ": ", conditionMessage(e))
        NULL
      }
    )

    if (!is.null(result)) {
      qsarray_list[[study]] <- result
    }
  }

  if (length(qsarray_list) == 0) {
    message("  No study produced a valid QuSAGE result — skipping.")
    return(NULL)
  }

  # Meta-analyse across studies (Meng et al. 2019); a single contributing
  # study needs no combination.
  meta_result <- if (length(qsarray_list) > 1) {
    qusage::combinePDFs(QSarrayList = qsarray_list)
  } else {
    qsarray_list[[1]]
  }

  raw_pvals         <- stats::setNames(qusage::pdf.pVal(meta_result), colnames(meta_result$path.PDF))
  activation_scores <- stats::setNames(meta_result$path.mean, colnames(meta_result$path.PDF))

  # Align to the full BTM gene set list (gene sets dropped above get NA), so
  # output is directly comparable to run_dearseq_comparison()'s
  # BTM-length p-value/score vectors.
  full_names  <- BTM[["geneset.names"]]
  pvals_full  <- stats::setNames(rep(NA_real_, length(full_names)), full_names)
  scores_full <- stats::setNames(rep(NA_real_, length(full_names)), full_names)
  pvals_full[names(raw_pvals)]   <- raw_pvals
  scores_full[names(activation_scores)] <- activation_scores

  exprmat_post <- t(na.omit(
    dplyr::filter(df, time_post_last_vax == day) |> dplyr::select(dplyr::all_of(gene_names))
  ))
  response <- df |>
    dplyr::filter(time_post_last_vax == day) |>
    dplyr::pull(immResp_MFC_anyAssay_log2_MFC)

  cor_res <- calculate_gs_correlation(y = exprmat_post, response = response, GSA = BTM)

  sig_pct <- 100 * mean(stats::p.adjust(pvals_full, method = "BH") < 0.05, na.rm = TRUE)
  message(sprintf("  %.1f%% of gene sets significant (BH-adjusted p < 0.05); %d studies combined",
                  sig_pct, length(qsarray_list)))

  list(
    pvals = list(rawPval = unname(pvals_full)),
    score = list(
      activation.scores = as.list(unname(scores_full)),
      fc.scores         = as.list(unname(scores_full))
    ),
    cor = cor_res
  )
}

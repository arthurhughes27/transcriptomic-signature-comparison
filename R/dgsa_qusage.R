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
#   - Sample pairing defaults to the shared filter_paired_samples() helper
#     (R/dgsa_common.R) - the same one dearseq uses - rather than the
#     original script's separate keep_most_recent_prevax_time_point()
#     applied globally across the whole dataset. This guarantees both
#     methods analyse literally the same set of samples for a given
#     comparison, as intended by Chapter 2, Section 2.2.2.
#
#     IMPORTANT: this is a *statistically meaningful* difference, not just a
#     bookkeeping one. The original script's per-comparison eset
#     (qusage_dgsa.R lines ~971-979) filters only by pathogen/vaccine_type/
#     study_accession and (time <= 0 | time == day) - it does not require a
#     participant to actually have a day-X sample to be included. Since
#     keep_most_recent_prevax_time_point(drop_postvax = FALSE) retains every
#     participant's baseline plus *all* their post-vaccination samples (any
#     day), a study with participants spread across several post-vaccination
#     days pools every one of those participants' baselines into *every*
#     day's comparison for that study - including participants never
#     measured at that specific day. That inflates the "pre" group with
#     unpaired, extraneous samples, which shifts its mean and inflates its
#     variance relative to a strictly-paired comparison, weakening
#     (increasing) the original script's p-values. filter_paired_samples()
#     requires both timepoints per participant, so run_qusage_comparison()
#     is systematically more powerful (smaller p-values) than the original
#     script whenever this occurs - which is common, since most studies
#     measure participants at several post-vaccination days. This is why
#     rawPval/activation.score did not match under `sample_scope = "paired"`
#     when validated against the original script's output.
#
#     Set `sample_scope = "study"` to reproduce the original script's
#     sample selection exactly (see filter_study_baseline_samples()) - this
#     is provided only to validate that the QuSAGE invocation itself
#     matches the original script once the same samples are used; the main
#     specification-analysis pipeline should always use the "paired"
#     default, which is both the methodologically correct choice and the
#     one consistent with dearseq and Chapter 2, Section 2.2.2.
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

#' Filter samples the way the original qusage_dgsa.R did
#'
#' Reproduces the sample selection of the original Hagan-et-al.-derived
#' script's per-comparison `cur_eset` (lines ~971-979): every participant's
#' most recent baseline (pre-vaccination) sample is kept, together with
#' every day = `day` post-vaccination sample, for the given vaccine -
#' regardless of whether a given participant has both. This can leave some
#' baseline samples "unpaired" (no matching day = `day` sample for that
#' participant), unlike [filter_paired_samples()] (R/dgsa_common.R), which
#' requires both and is what [run_qusage_comparison()] uses by default.
#'
#' Provided only to validate [run_qusage_comparison()] against the original
#' script's output (`sample_scope = "study"`) - see the design notes at the
#' top of this file. The main specification-analysis pipeline should always
#' use the "paired" default.
#'
#' @param hipc Merged clinical/expression tibble.
#' @param vax Vaccine name to filter to (matches `hipc$vaccine_name`).
#' @param day Post-vaccination timepoint (matches `hipc$time_post_last_vax`).
#'
#' @return A tibble with one row per retained sample, `participant_id`
#'   recoded to a compact factor.
#' @keywords internal
filter_study_baseline_samples <- function(hipc, vax, day) {
  hipc |>
    dplyr::filter(vaccine_name == vax) |>
    dplyr::group_by(participant_id) |>
    dplyr::mutate(
      max_pre_time = suppressWarnings(max(study_time_collected[time_post_last_vax <= 0], na.rm = TRUE))
    ) |>
    dplyr::filter(
      time_post_last_vax == day |
        (time_post_last_vax <= 0 & study_time_collected == max_pre_time)
    ) |>
    dplyr::select(-max_pre_time) |>
    dplyr::ungroup() |>
    dplyr::mutate(participant_id = as.factor(as.numeric(as.factor(participant_id))))
}

#' Run one QuSAGE comparison (vaccine x timepoint)
#'
#' Filters to samples for `vax` at `day` (by default, paired pre-/
#' post-vaccination samples, via [filter_paired_samples()]; see
#' `sample_scope`), restricts to genes complete across those samples and
#' gene sets retaining at least 2 such genes ([keep_symbols()]), runs
#' `qusage::qusage()` independently within each contributing study, and
#' meta-analyses the resulting per-study QSarray objects across studies via
#' `qusage::combinePDFs()` (a single study's QSarray is used directly with
#' no meta-analysis needed). Also computes gene-set correlations with
#' immune response ([calculate_gs_correlation()]).
#'
#' @param vax Vaccine name.
#' @param day Post-vaccination timepoint (days).
#' @param hipc Merged clinical/expression tibble.
#' @param BTM Gene set list (genesets + geneset.names).
#' @param gene_names Character vector of gene expression column names in
#'   `hipc`.
#' @param equal_variance QuSAGE `var.equal` argument: whether to assume
#'   equal variance between the pre- and post-vaccination groups. FALSE
#'   (default) is the baseline specification (Table 2.1). Set to NULL to
#'   omit the argument entirely and use `qusage::qusage()`'s own default,
#'   matching the original script (which never passed `var.equal`).
#' @param sample_scope "paired" (default): require each participant to have
#'   both the baseline and the day = `day` sample
#'   ([filter_paired_samples()]) - the methodologically correct choice, and
#'   the one consistent with dearseq. "study": reproduce the original
#'   script's more permissive sample selection
#'   ([filter_study_baseline_samples()]), for validating this function
#'   against the original script's output only.
#'
#' @return A list with elements `pvals` (rawPval, one value per BTM gene
#'   set, NA for gene sets dropped by the completeness filter), `score`
#'   (activation.scores/fc.scores - QuSAGE has no separate fold-change
#'   score, so both equal the QSarray `path.mean`), and `cor` (gene-set
#'   correlations with immune response); or NULL if there is insufficient
#'   data for this comparison.
run_qusage_comparison <- function(vax, day, hipc, BTM, gene_names, equal_variance = FALSE,
                                  sample_scope = c("paired", "study")) {
  sample_scope <- match.arg(sample_scope)

  df <- if (sample_scope == "paired") {
    filter_paired_samples(hipc, vax, day)
  } else {
    filter_study_baseline_samples(hipc, vax, day)
  }

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

    qusage_args <- list(
      eset       = exprmat,
      labels     = labels,
      contrast   = "post-pre",
      geneSets   = genesets,
      pairVector = as.character(df_study$participant_id)
    )
    # equal_variance = NULL omits var.equal entirely, matching the original
    # script (which never passed it, relying on qusage::qusage()'s own
    # default) - used when validating with sample_scope = "study".
    if (!is.null(equal_variance)) {
      qusage_args$var.equal <- equal_variance
    }

    result <- tryCatch(
      do.call(qusage::qusage, qusage_args),
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

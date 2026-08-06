# =============================================================================
# Modular dearseq differential gene-set analysis
# =============================================================================
# run_dearseq_comparison() runs one vaccine x timepoint comparison with
# dearseq::dgsa_seq(), parameterised by the dearseq-specific hyperparameters
# from Table 2.1 (Chapter 2): the covariate adjustment set, and the
# mean-variance weighting method/level. Depends on the helpers in
# R/dgsa_common.R (filter_paired_samples, build_covariate_matrix,
# calculate_gs_correlation).
# =============================================================================

#' Calculate gene-set activation and fold-change scores
#'
#' For each gene set, computes per-sample dearseq-weighted expression
#' scores, gene-level activation scores (difference in weighted expression
#' between the active and inactive `phi` levels) and fold-change scores
#' (paired difference when `sample_group` pairs samples one-to-one across
#' `phi` levels), and their gene-set-level averages.
#'
#' @keywords internal
calculate_scores <- function(y,
                             x = NULL,
                             phi,
                             GSA,
                             use_phi = TRUE,
                             preprocessed = TRUE,
                             gene_based = FALSE,
                             bw = c("nrd", "ucv", "SJ", "nrd0", "bcv"),
                             kernel = c(
                               "gaussian",
                               "epanechnikov",
                               "rectangular",
                               "triangular",
                               "biweight",
                               "tricube",
                               "cosine",
                               "optcosine"
                             ),
                             transform = TRUE,
                             verbose = TRUE,
                             na.rm = TRUE,
                             active_level,
                             sample_group = NULL) {
  ## dimensions & validity checks

  stopifnot(is.matrix(y))
  stopifnot(is.matrix(x))
  stopifnot(is.null(phi) | is.matrix(phi))

  g <- nrow(y)  # the number of genes measured
  n <- ncol(y)  # the number of samples measured
  qq <- ncol(x)  # the number of covariates
  stopifnot(nrow(x) == n)
  if (use_phi) {
    stopifnot(nrow(phi) == n)
  }

  # If sample_group given, check that data is exactly paired
  if (!is.null(sample_group)) {
    # Tabulate phi and sample_group to check levels for each individual
    tab <- table(sample_group, phi)
    # All elements must be 1 in this table if data is truly paired
    is_paired <- all(tab == 1)
    if (!is_paired) {
      message(
        "You have provided individuals but the data is not paired.
              Ignoring this when calculating the fold-change and proceeding in the default manner."
      )
    }
  }

  # removing genes never observed:
  observed <- which(rowSums(y, na.rm = TRUE) != 0)
  nb_g_sum0 <- length(observed) - g
  if (nb_g_sum0 > 0) {
    warning(nb_g_sum0,
            " y rows sum to 0 (i.e. are never observed)",
            "and have been removed")
  }

  if (preprocessed) {
    y_lcpm <- t(y[observed, ])
  } else {
    # transforming raw counts to log-counts per million (lcpm)
    y_lcpm <- t(apply(y[observed, ], MARGIN = 2, function(v) {
      log2((v + 0.5) / (sum(v) + 1) * 10^6)
    }))
  }
  N <- length(y_lcpm)
  p <- ncol(y_lcpm)

  # fitting OLS to the lcpm
  if (na.rm) {
    y_lcpm0 <- y_lcpm
    y_lcpm0[is.na(y_lcpm0)] <- 0
    B_ols <- solve(crossprod(x)) %*% t(x) %*% y_lcpm0
  } else {
    B_ols <- solve(crossprod(x)) %*% t(x) %*% y_lcpm
  }
  mu <- x %*% B_ols

  if (gene_based) {
    mu_avg <- colMeans(mu, na.rm = na.rm)
    mu_x <- mu_avg
  } else {
    mu_x <- mu
    mu_x[is.na(y_lcpm)] <- NA
  }

  y_T <- t(y)
  yt_mu <- y_T - mu_x

  # calculate weights
  weights <- dearseq::sp_weights(y,
                                 x,
                                 phi,
                                 use_phi,
                                 preprocessed,
                                 gene_based,
                                 bw,
                                 kernel,
                                 transform,
                                 verbose,
                                 na.rm)

  weights_all <- weights$weights

  yt_mu_standardised <- yt_mu * t(weights_all)

  expr_mat <- as.matrix(yt_mu_standardised)

  # 2) Your gene‐sets, with names
  genesets <- GSA[["genesets"]]                   # list of character vectors
  names(genesets) <- GSA[["geneset.names"]]

  # 3) Precompute, for each set, the column‐indices into expr_mat
  geneset_idx <- lapply(genesets, function(gs) {
    which(colnames(expr_mat) %in% gs)
  })

  phi_inactive <- phi != active_level
  phi_active <- phi == active_level

  gene_scores_all <-
    colMeans(expr_mat[phi_active, , drop = FALSE], na.rm = TRUE) -
    colMeans(expr_mat[phi_inactive, , drop = FALSE], na.rm = TRUE)

  # Calculate gene-level FC scores
  if (is_paired) {
    #If paired, first calculate individual-gene-level scores
    active_idx   <- which(phi == active_level)
    inactive_idx <- which(phi != active_level)

    ids_active   <- sample_group[active_idx]
    ids_inactive <- sample_group[inactive_idx]

    # 2. order each so that rows line up by individual
    ord_active   <- order(ids_active)
    ord_inactive <- order(ids_inactive)

    Y_act <- y_T[active_idx[ord_active], , drop = FALSE]
    Y_in <- y_T[inactive_idx[ord_inactive], , drop = FALSE]

    # sanity check: same individual order
    stopifnot(all(ids_active[ord_active] == ids_inactive[ord_inactive]))

    # 3. per‐individual fold‐change matrix (rows = individuals, cols = predictors)
    ind_fc_mat <- Y_act - Y_in

    # 4. predictor‐level score = average across individuals
    fc_scores_all <- colMeans(ind_fc_mat)

  } else {
    #If not paired, take difference in gene averages between binary groups
    fc_scores_all <- colMeans(y_T[phi_active, , drop = FALSE], na.rm = TRUE) -
      colMeans(y_T[phi_inactive, , drop = FALSE], na.rm = TRUE)
  }

  # === Extract per‐geneset results ===

  # 1) individual.scores: list of sub‐matrices (n_samples × |geneset|)
  individual.scores <- lapply(geneset_idx, function(idxs) {
    expr_mat[, idxs, drop = FALSE]
  })

  # 2) gene.scores: list of per‐gene score vectors
  gene.scores <- lapply(geneset_idx, function(idxs) {
    gene_scores_all[idxs]
  })

  fc_gene.scores <- lapply(geneset_idx, function(idxs) {
    fc_scores_all[idxs]
  })

  # 3) activation.scores: list of single values = mean of gene.scores
  activation.scores <- lapply(gene.scores, function(gs) {
    mean(gs, na.rm = TRUE)
  })

  fc.scores <- lapply(fc_gene.scores, function(gs) {
    mean(gs, na.rm = TRUE)
  })

  list(
    individual.scores = individual.scores,
    gene.scores       = gene.scores,
    activation.scores = activation.scores,
    fc.scores         = fc.scores
  )
}

#' Run dearseq::dgsa_seq() with shared arguments
#'
#' Switches between the permutation test (small samples) and the asymptotic
#' test (large samples), as in the original analysis. `which_weights` and
#' `gene_based_weights` are the two dearseq-specific hyperparameters varied
#' in the specification analysis (Table 2.1): the mean-variance weighting
#' method ("loclin"/"voom") and estimation level (gene- vs
#' observation-level).
#'
#' @keywords internal
run_dgsa_seq <- function(exprmat, x, phi, subject_ids, BTM,
                         which_weights = "loclin",
                         gene_based_weights = FALSE) {
  shared_args <- list(
    exprmat             = exprmat,
    covariates           = x,
    variables2test       = as.matrix(phi),
    genesets             = BTM[["genesets"]],
    sample_group         = subject_ids,
    which_weights        = which_weights,
    progressbar          = FALSE,
    parallel_comp        = TRUE,
    preprocessed         = TRUE,
    gene_based_weights   = gene_based_weights,
    transform            = TRUE,
    padjust_methods      = "BH",
    bw                   = "nrd",
    kernel               = "gaussian",
    homogen_traj         = FALSE,
    na.rm_gsaseq         = TRUE,
    verbose              = FALSE
  )

  if (ncol(exprmat) < 200) {
    shared_args <- c(shared_args, list(
      which_test              = "permutation",
      weights_var2test_condi  = FALSE,
      n_perm                  = 1000,
      adaptive                = TRUE,
      nb_cores                = 1
    ))
  } else {
    shared_args <- c(shared_args, list(
      which_test              = "asymptotic",
      weights_var2test_condi  = TRUE,
      nb_cores                = 7
    ))
  }

  do.call(dearseq::dgsa_seq, shared_args)
}

#' Run one dearseq comparison (vaccine x timepoint)
#'
#' Filters to paired pre-/post-vaccination samples for `vax` at `day`
#' ([filter_paired_samples()]), builds the covariate design matrix
#' ([build_covariate_matrix()]), runs the dearseq variance-component test
#' ([run_dgsa_seq()]), and computes gene-set activation/fold-change scores
#' and correlations with immune response.
#'
#' @param vax Vaccine name.
#' @param day Post-vaccination timepoint (days).
#' @param hipc Merged clinical/expression tibble, with `vaccine_code`
#'   (pre/post indicator) and `study_accession` (factor) already added.
#' @param BTM Gene set list (genesets + geneset.names).
#' @param gene_names Character vector of gene expression column names in
#'   `hipc`.
#' @param covariates Character vector of covariate columns to adjust for
#'   (see [build_covariate_matrix()]); default matches the baseline
#'   specification (age, sex, study, race).
#' @param which_weights dearseq mean-variance weighting method: "loclin"
#'   (default, baseline) or "voom".
#' @param gene_based_weights Weight estimation level: FALSE for
#'   observation-level (default, baseline), TRUE for gene-level.
#'
#' @return A list with elements `pvals` (dearseq p-value table), `score`
#'   (gene-set activation/fold-change scores), and `cor` (gene-set
#'   correlations with immune response); or NULL if there is insufficient
#'   data for this comparison.
run_dearseq_comparison <- function(vax, day, hipc, BTM, gene_names,
                                   covariates = c("age_imputed", "gender", "study_accession", "race"),
                                   which_weights = "loclin",
                                   gene_based_weights = FALSE) {

  df <- filter_paired_samples(hipc, vax, day)

  if (nrow(df) == 0) {
    message("  No data after filtering — skipping.")
    return(NULL)
  }

  subject_ids <- df$participant_id
  phi         <- as.numeric(df$vaccine_code)
  x           <- build_covariate_matrix(df, covariates = covariates)

  # Expression matrices: full (pre + post) for testing/scoring; post only for
  # correlation with immune response
  exprmat      <- t(na.omit(dplyr::select(df, dplyr::all_of(gene_names))))
  exprmat_post <- t(na.omit(
    dplyr::filter(df, time_post_last_vax == day) |> dplyr::select(dplyr::all_of(gene_names))
  ))
  response     <- df |>
    dplyr::filter(time_post_last_vax == day) |>
    dplyr::pull(immResp_MFC_anyAssay_log2_MFC)

  dgsa_result <- NULL
  suppressMessages(
    utils::capture.output(
      dgsa_result <- run_dgsa_seq(exprmat, x, phi, subject_ids, BTM,
                                  which_weights = which_weights,
                                  gene_based_weights = gene_based_weights),
      type = "message"
    )
  )

  score_res <- calculate_scores(
    y            = exprmat,
    x,
    GSA          = BTM,
    phi          = as.matrix(phi),
    use_phi      = TRUE,
    preprocessed = TRUE,
    gene_based   = FALSE,
    bw           = "nrd",
    active_level = 2,
    sample_group = subject_ids
  )

  cor_res <- calculate_gs_correlation(
    y        = exprmat_post,
    response = response,
    GSA      = BTM
  )

  sig_pct <- 100 * mean(dgsa_result$pvals$adjPval < 0.05)
  message(sprintf("  %.1f%% of gene sets significant (BH-adjusted p < 0.05)", sig_pct))

  list(pvals = dgsa_result$pvals, score = score_res, cor = cor_res)
}

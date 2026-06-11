# =============================================================================
# Differential Gene-Set Analysis with dearseq
# =============================================================================
# One DGSA is run per vaccine x timepoint combination, comparing post-vaccination
# samples against self-controls (pre-vaccination baseline).
# Covariates: age, sex, study (with automatic collinearity removal).
# Results are saved incrementally to guard against mid-run crashes.
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────

library(fs)
library(tidyverse)
library(dearseq)

# ── Paths ─────────────────────────────────────────────────────────────────────

p_data_expr  <- fs::path("data", "hipc_merged_young_noNorm.rds")
p_data_btm   <- fs::path("data", "BTM_processed.rds")
p_results    <- fs::path("output", "results", "reanalysis", "dearseq_dgsa_results_list.rds")

# ── Load data ─────────────────────────────────────────────────────────────────

hipc <- readRDS(p_data_expr)
BTM  <- readRDS(p_data_btm)

gene_names <- hipc %>% select(a1cf:zzz3) %>% colnames()

# Label pre- (1) and post-vaccination (2) samples; factorise study
hipc <- hipc %>%
  mutate(
    vaccine_code    = as.factor(if_else(time_post_last_vax > 0, 2, 1)),
    study_accession = as.factor(study_accession)
  )

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Calculate gene set activity scores
calculate_scores = function(y,
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
  
  # kernel <- match.arg(kernel)
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
  
  y_T = t(y)
  yt_mu = y_T - mu_x
  
  # calculate weights
  weights = sp_weights(y,
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
  
  weights_all = weights$weights
  
  yt_mu_standardised = yt_mu * t(weights_all)
  
  expr_mat <- as.matrix(yt_mu_standardised)
  
  # 2) Your gene‐sets, with names
  genesets <-
    GSA[["genesets"]]                   # list of character vectors
  names(genesets) <- GSA[["geneset.names"]]
  
  # 3) Precompute, for each set, the column‐indices into expr_mat
  geneset_idx <- lapply(genesets, function(g) {
    which(colnames(expr_mat) %in% g)
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
    fc_scores_all = colMeans(y_T[phi_active, , drop = FALSE], na.rm = TRUE) -
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
  
  result <- list(
    individual.scores = individual.scores,
    gene.scores       = gene.scores,
    activation.scores = activation.scores,
    fc.scores = fc.scores
  )
}

# Calculate gene-set correlation
calculate_gs_correlation = function(y, response, GSA) {
  ### VALIDITY CHECKS ###
  
  response = response %>% as.matrix()
  
  # Check that number of samples are the same in y and response
  n = ncol(y)
  if (length(response) != n) {
    stop("Number of response values given does not match number of samples in y")
  }
  
  mean.corr <- sapply(seq_along(GSA[["genesets"]]), function(i) {
    gs <- GSA[["genesets"]][[i]]
    # keep only genes actually in y
    keep <- intersect(gs, rownames(y))
    if (length(keep) == 0) {
      return(NA_real_)  # or 0, or whatever you’d prefer
    }
    # compute Spearman cor for each gene
    cors <- apply(y[keep, , drop = FALSE], 1, function(gene_expr)
      cor(gene_expr, response, method = "spearman", use = "pairwise.complete.obs"))
    # average them
    mean(cors, na.rm = TRUE)
  })
  
  # assign names
  names(mean.corr) <- GSA[["geneset.names"]]
  
  corr.mean <- sapply(seq_along(GSA[["genesets"]]), function(i) {
    # 1. pick your genes
    gs    <- GSA[["genesets"]][[i]]
    keep  <- intersect(gs, rownames(y))
    if (length(keep) == 0)
      return(NA_real_)
    
    # 2. subset and standardize each row (gene) across samples
    y_sub <- y[keep, , drop = FALSE]
    #    t(scale(t(y_sub))) gives genes×samples Z-scores
    zsub  <- t(scale(t(y_sub), center = TRUE, scale = TRUE))
    
    # 3. meta-gene: mean Z per sample
    meta  <- colMeans(zsub, na.rm = TRUE)
    
    # 4. Spearman correlation with response
    cor(meta, response, method = "spearman", use = "pairwise.complete.obs")
  })
  
  # assign names
  names(corr.mean) <- GSA[["geneset.names"]]
  
  return(list(mean.corr = mean.corr, corr.mean = corr.mean))
}

# Build a covariate matrix adjusting for age, gender, and study, then remove
# any linearly dependent columns and the vaccine_code columns themselves,
# leaving only the adjustment covariates for dearseq.
build_covariate_matrix <- function(df) {
  X_full  <- model.matrix(~ vaccine_code + age_imputed + gender + study_accession,
                          data = df)
  qrX     <- qr(X_full)
  X_indep <- X_full[, sort(qrX$pivot[seq_len(qrX$rank)]), drop = FALSE]
  x       <- X_indep[, !grepl("^vaccine_code", colnames(X_indep)), drop = FALSE]
  message("Adjusting for covariates: ", paste(colnames(x), collapse = ", "))
  x
}

# Run dgsa_seq with shared arguments, switching between permutation (small
# samples) and asymptotic (large samples) tests based on sample size.
run_dgsa <- function(exprmat, x, phi, subject_ids, BTM) {
  shared_args <- list(
    exprmat                = exprmat,
    covariates             = x,
    variables2test         = as.matrix(phi),
    genesets               = BTM[["genesets"]],
    sample_group           = subject_ids,
    which_weights          = "loclin",
    progressbar            = FALSE,
    parallel_comp          = TRUE,
    preprocessed           = TRUE,
    gene_based_weights     = FALSE,
    transform              = TRUE,
    padjust_methods        = "BH",
    bw                     = "nrd",
    kernel                 = "gaussian",
    homogen_traj           = FALSE,
    na.rm_gsaseq           = TRUE,
    verbose                = FALSE
  )
  
  if (ncol(exprmat) < 200) {
    shared_args <- c(shared_args, list(
      which_test             = "permutation",
      weights_var2test_condi = FALSE,
      n_perm                 = 1000,
      adaptive               = TRUE,
      nb_cores               = 1
    ))
  } else {
    shared_args <- c(shared_args, list(
      which_test             = "asymptotic",
      weights_var2test_condi = TRUE,
      nb_cores               = 7
    ))
  }
  
  do.call(dgsa_seq, shared_args)
}

# Run a full comparison for one vaccine x timepoint combination.
# Returns a named list with p-values, gene-set scores, and response correlations,
# or NULL if the participant filter leaves insufficient data.
run_comparison <- function(vax, day, hipc, BTM, gene_names) {
  
  # Keep participants with both a pre-vaccination and this post-vaccination sample;
  # for participants with multiple pre-vaccination samples, use the latest one
  df <- hipc %>%
    filter(
      vaccine_name == vax,
      time_post_last_vax == day | time_post_last_vax <= 0
    ) %>%
    group_by(participant_id) %>%
    filter(
      any(time_post_last_vax == day),
      any(time_post_last_vax <= 0)
    ) %>%
    mutate(max_pre_time = max(study_time_collected[time_post_last_vax <= 0],
                              na.rm = TRUE)) %>%
    filter(time_post_last_vax == day | study_time_collected == max_pre_time) %>%
    select(-max_pre_time) %>%
    ungroup() %>%
    mutate(participant_id = as.factor(as.numeric(as.factor(participant_id))))
  
  if (nrow(df) == 0) {
    message("  No data after filtering — skipping.")
    return(NULL)
  }
  
  subject_ids <- df$participant_id
  phi         <- as.numeric(df$vaccine_code)
  x           <- build_covariate_matrix(df)
  
  # Expression matrices: full (pre + post) for testing/scoring; post only for
  # correlation with immune response
  exprmat      <- t(na.omit(select(df, all_of(gene_names))))
  exprmat_post <- t(na.omit(filter(df, time_post_last_vax == day) %>%
                              select(all_of(gene_names))))
  response     <- df %>%
    filter(time_post_last_vax == day) %>%
    pull(immResp_MFC_anyAssay_log2_MFC)
  
  dgsa_result <- NULL
  suppressMessages(
    capture.output(
      dgsa_result <- run_dgsa(exprmat, x, phi, subject_ids, BTM),
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

# =============================================================================
# MAIN LOOP
# =============================================================================

# Load any previously saved results so a restart continues rather than reruns
results_list <- if (file_exists(p_results)) readRDS(p_results) else list()

timepoints <- hipc %>%
  filter(time_post_last_vax > 0) %>%
  pull(time_post_last_vax) %>%
  unique() %>%
  sort()

for (day in timepoints) {
  
  # Vaccines with at least one participant having both pre- and post-vax samples
  valid_vaccines <- hipc %>%
    group_by(vaccine_name, participant_id) %>%
    summarise(
      has_post = any(time_post_last_vax == day),
      has_pre  = any(time_post_last_vax <= 0),
      .groups  = "drop"
    ) %>%
    filter(has_post & has_pre) %>%
    pull(vaccine_name) %>%
    unique()
  
  if (length(valid_vaccines) == 0) next
  
  for (vax in valid_vaccines) {
    
    comparison_name <- sprintf("%s vs Control - Day %s", vax, day)
    
    # Skip comparisons already completed in a previous run
    if (comparison_name %in% names(results_list)) {
      message("Skipping (already complete): ", comparison_name)
      next
    }
    
    message("Running: ", comparison_name)
    
    result <- tryCatch(
      run_comparison(vax, day, hipc, BTM, gene_names),
      error = function(e) {
        message("  ERROR: ", conditionMessage(e))
        NULL
      }
    )
    
    results_list[[comparison_name]] <- result
    
    # Save after every comparison so a crash loses at most one result
    saveRDS(object = results_list, 
            file = p_results)
    
    # Free memory before the next iteration
    gc()
  }
}

message("Analysis complete. Results saved to: ", p_results)

rm(list = ls())
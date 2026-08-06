# These tests exercise the real dearseq::dgsa_seq() call end-to-end on the
# synthetic dataset, so they're skipped where the (Bioconductor) dearseq
# package isn't installed.

test_that("run_dearseq_comparison() runs end-to-end on synthetic data", {
  testthat::skip_if_not_installed("dearseq")

  d <- make_synthetic_is2_data(seed = 3, n_participants_per_study = 10, n_genes = 40)
  hipc <- d$hipc |>
    dplyr::mutate(
      vaccine_code    = as.factor(ifelse(time_post_last_vax > 0, 2, 1)),
      study_accession = as.factor(study_accession)
    )
  gene_names <- setdiff(colnames(hipc), c(
    "participant_id", "vaccine_name", "study_accession", "time_post_last_vax",
    "study_time_collected", "age_imputed", "gender", "race",
    "immResp_MFC_anyAssay_log2_MFC", "vaccine_code"
  ))

  result <- run_dearseq_comparison(
    vax = "Vaccine A", day = 7, hipc = hipc, BTM = d$BTM, gene_names = gene_names
  )

  expect_false(is.null(result))
  expect_named(result, c("pvals", "score", "cor"))
  expect_equal(nrow(result$pvals), length(d$BTM$geneset.names))
  expect_true(all(result$pvals$rawPval >= 0 & result$pvals$rawPval <= 1, na.rm = TRUE))
})

test_that("run_dearseq_comparison() respects the covariates argument", {
  testthat::skip_if_not_installed("dearseq")

  d <- make_synthetic_is2_data(seed = 4, n_participants_per_study = 10, n_genes = 40)
  hipc <- d$hipc |>
    dplyr::mutate(
      vaccine_code    = as.factor(ifelse(time_post_last_vax > 0, 2, 1)),
      study_accession = as.factor(study_accession)
    )
  gene_names <- setdiff(colnames(hipc), c(
    "participant_id", "vaccine_name", "study_accession", "time_post_last_vax",
    "study_time_collected", "age_imputed", "gender", "race",
    "immResp_MFC_anyAssay_log2_MFC", "vaccine_code"
  ))

  # No adjustment at all: covariate matrix should be intercept-only, and the
  # comparison should still run to completion.
  result <- run_dearseq_comparison(
    vax = "Vaccine A", day = 7, hipc = hipc, BTM = d$BTM, gene_names = gene_names,
    covariates = character(0)
  )

  expect_false(is.null(result))
  expect_equal(nrow(result$pvals), length(d$BTM$geneset.names))
})

test_that("run_dearseq_comparison() returns NULL when there is no paired data", {
  d <- make_synthetic_is2_data(seed = 6, n_participants_per_study = 5, n_genes = 40)
  hipc <- d$hipc |>
    dplyr::mutate(
      vaccine_code    = as.factor(ifelse(time_post_last_vax > 0, 2, 1)),
      study_accession = as.factor(study_accession)
    )
  gene_names <- setdiff(colnames(hipc), c(
    "participant_id", "vaccine_name", "study_accession", "time_post_last_vax",
    "study_time_collected", "age_imputed", "gender", "race",
    "immResp_MFC_anyAssay_log2_MFC", "vaccine_code"
  ))

  # day 99 does not exist in the synthetic data
  result <- run_dearseq_comparison(
    vax = "Vaccine A", day = 99, hipc = hipc, BTM = d$BTM, gene_names = gene_names
  )

  expect_null(result)
})

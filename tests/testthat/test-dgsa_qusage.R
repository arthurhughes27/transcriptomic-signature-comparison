test_that("keep_symbols() restricts gene sets and drops undersized ones", {
  genesets <- list(gs1 = c("g1", "g2", "g3"), gs2 = c("g4"), gs3 = c("g1", "g5"))
  out <- keep_symbols(genesets, symbols2keep = c("g1", "g2", "g5"), min_size = 2)

  expect_equal(names(out), c("gs1", "gs3"))
  expect_equal(out$gs1, c("g1", "g2"))
  expect_equal(out$gs3, c("g1", "g5"))
})

test_that("keep_symbols() keeps a gene set exactly at min_size", {
  genesets <- list(gs1 = c("g1", "g2"))
  out <- keep_symbols(genesets, symbols2keep = c("g1", "g2"), min_size = 2)
  expect_equal(names(out), "gs1")
})

test_that("filter_study_baseline_samples() includes unpaired baseline samples, unlike filter_paired_samples()", {
  hipc <- tibble::tibble(
    vaccine_name         = c("V1", "V1", "V1", "V1", "V1"),
    participant_id       = c("p1", "p1", "p2", "p2", "p3"),
    time_post_last_vax   = c(-1,   7,    -1,   3,    -1),
    study_time_collected = c(-1,   7,    -1,   3,    -1)
  )
  # p1: pre + day7 sample -> paired for day 7
  # p2: pre + day3 sample only (no day7) -> baseline still pooled into the
  #     day-7 comparison under "study" scope (unpaired there)
  # p3: pre sample only, no post-vax sample at all -> baseline still pooled
  #     into the day-7 comparison under "study" scope, matching the
  #     original script (which pools every participant's baseline into
  #     every day's comparison for their study, regardless of whether that
  #     participant has any post-vaccination sample)

  paired <- filter_paired_samples(hipc, "V1", 7)
  study  <- filter_study_baseline_samples(hipc, "V1", 7)

  expect_equal(nrow(paired), 2)
  expect_equal(length(unique(paired$participant_id)), 1)

  expect_equal(nrow(study), 4)
  expect_equal(length(unique(study$participant_id)), 3)
})

# The following tests exercise the real qusage::qusage() call end-to-end on
# the synthetic dataset, so they're skipped where qusage isn't installed.

qusage_gene_names <- function(hipc) {
  setdiff(colnames(hipc), c(
    "participant_id", "vaccine_name", "study_accession", "time_post_last_vax",
    "study_time_collected", "age_imputed", "gender", "race",
    "immResp_MFC_anyAssay_log2_MFC"
  ))
}

test_that("run_qusage_comparison() runs end-to-end on synthetic data", {
  testthat::skip_if_not_installed("qusage")

  d <- make_synthetic_is2_data(seed = 3, n_participants_per_study = 10, n_genes = 40)
  gene_names <- qusage_gene_names(d$hipc)

  result <- run_qusage_comparison(
    vax = "Vaccine A", day = 7, hipc = d$hipc, BTM = d$BTM, gene_names = gene_names
  )

  expect_false(is.null(result))
  expect_named(result, c("pvals", "score", "cor"))
  expect_equal(length(result$pvals$rawPval), length(d$BTM$geneset.names))
  expect_true(all(result$pvals$rawPval >= 0 & result$pvals$rawPval <= 1, na.rm = TRUE))
})

test_that("run_qusage_comparison() respects the equal_variance argument", {
  testthat::skip_if_not_installed("qusage")

  d <- make_synthetic_is2_data(seed = 4, n_participants_per_study = 10, n_genes = 40)
  gene_names <- qusage_gene_names(d$hipc)

  result <- run_qusage_comparison(
    vax = "Vaccine A", day = 7, hipc = d$hipc, BTM = d$BTM, gene_names = gene_names,
    equal_variance = TRUE
  )

  expect_false(is.null(result))
  expect_equal(length(result$pvals$rawPval), length(d$BTM$geneset.names))
})

test_that("run_qusage_comparison() runs with sample_scope = 'study' and equal_variance = NULL", {
  testthat::skip_if_not_installed("qusage")

  d <- make_synthetic_is2_data(seed = 9, n_participants_per_study = 10, n_genes = 40)
  gene_names <- qusage_gene_names(d$hipc)

  result <- run_qusage_comparison(
    vax = "Vaccine A", day = 7, hipc = d$hipc, BTM = d$BTM, gene_names = gene_names,
    sample_scope = "study", equal_variance = NULL
  )

  expect_false(is.null(result))
  expect_equal(length(result$pvals$rawPval), length(d$BTM$geneset.names))
})

test_that("run_qusage_comparison() returns NULL when there is no paired data", {
  d <- make_synthetic_is2_data(seed = 6, n_participants_per_study = 5, n_genes = 40)
  gene_names <- qusage_gene_names(d$hipc)

  # day 99 does not exist in the synthetic data
  result <- run_qusage_comparison(
    vax = "Vaccine A", day = 99, hipc = d$hipc, BTM = d$BTM, gene_names = gene_names
  )

  expect_null(result)
})

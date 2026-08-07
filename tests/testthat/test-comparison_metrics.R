test_that("cohens_kappa() returns 1 for perfect agreement and handles degenerate cases", {
  expect_equal(cohens_kappa(c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE)), 1)
  expect_true(is.na(cohens_kappa(logical(0), logical(0))))
  expect_true(is.na(cohens_kappa(c(TRUE, TRUE), c(TRUE, TRUE))))  # pe == 1 degenerate case
})

test_that("compute_concordance_metrics() computes expected agreement and correlation", {
  results_df <- tibble::tibble(
    comparison        = rep("V1 - Day 7", 4),
    condition         = "V1",
    time              = 7,
    gs.name           = c("gs1", "gs2", "gs1", "gs2"),
    method            = c("dearseq", "dearseq", "qusage", "qusage"),
    rawPval           = c(0.01, 0.5, 0.02, 0.6),
    fc.score              = c(1.0, 0.1, 0.9, 0.05),
    withinTime.adjPval_BH = c(0.01, 0.5, 0.02, 0.6)
  )

  out <- compute_concordance_metrics(results_df, alpha = 0.05)

  expect_equal(out$method1, "dearseq")
  expect_equal(out$method2, "qusage")
  expect_equal(out$overall$n_gene_sets, 2)
  expect_equal(out$overall$pval_cor, 1)
  expect_equal(out$overall$score_cor, 1)
  expect_equal(out$overall$pct_agree, 1)
  expect_equal(out$overall$kappa, 1)
  expect_equal(out$overall$n_both_sig, 1)
  expect_equal(out$overall$n_neither_sig, 1)
  expect_equal(out$overall$n_sig_method1_only, 0)
  expect_equal(out$overall$n_sig_method2_only, 0)

  expect_equal(nrow(out$by_comparison), 1)
  expect_equal(out$by_comparison$n_gene_sets, 2)
})

test_that("compute_concordance_metrics() detects method-specific significance", {
  results_df <- tibble::tibble(
    comparison        = rep("V1 - Day 7", 2),
    condition         = "V1",
    time              = 7,
    gs.name           = c("gs1", "gs1"),
    method            = c("dearseq", "qusage"),
    rawPval           = c(0.01, 0.5),
    fc.score              = c(1.0, 0.1),
    withinTime.adjPval_BH = c(0.01, 0.5)
  )

  out <- compute_concordance_metrics(results_df, alpha = 0.05)

  expect_equal(out$overall$n_sig_method1_only, 1)
  expect_equal(out$overall$n_sig_method2_only, 0)
  expect_equal(out$overall$pct_agree, 0)
})

test_that("pivot_methods_wide() produces one row per comparison x gene set", {
  results_df <- tibble::tibble(
    comparison = c("V1 - Day 7", "V1 - Day 7"),
    condition  = "V1",
    time       = 7,
    gs.name    = c("gs1", "gs1"),
    method     = c("dearseq", "qusage"),
    rawPval    = c(0.01, 0.02)
  )

  out <- pivot_methods_wide(results_df, value_cols = "rawPval")

  expect_equal(nrow(out), 1)
  expect_true(all(c("rawPval_dearseq", "rawPval_qusage") %in% colnames(out)))
  expect_equal(out$rawPval_dearseq, 0.01)
  expect_equal(out$rawPval_qusage, 0.02)
})

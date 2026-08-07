test_that("compute_pct_significant_by_comparison() computes percentages relative to n evaluated, excluding NA", {
  baseline_df <- tibble::tibble(
    condition             = c("V1", "V1", "V1", "V1", "V1"),
    time                   = c(1, 1, 1, 1, 1),
    method                 = c("dearseq", "dearseq", "dearseq", "dearseq", "qusage"),
    withinTime.adjPval_BH = c(0.01, 0.5, 0.2, NA, 0.01)
  )

  out <- compute_pct_significant_by_comparison(baseline_df)

  dearseq_row <- out[out$method == "dearseq", ]
  qusage_row  <- out[out$method == "qusage", ]

  expect_equal(dearseq_row$n_evaluated, 3)  # NA row excluded
  expect_equal(dearseq_row$n_significant, 1)
  expect_equal(dearseq_row$pct_significant, 100 / 3)

  expect_equal(qusage_row$n_evaluated, 1)
  expect_equal(qusage_row$n_significant, 1)
  expect_equal(qusage_row$pct_significant, 100)
})

test_that("compute_pct_significant_by_comparison() respects a custom pval_col and alpha", {
  baseline_df <- tibble::tibble(
    condition = "V1", time = 1, method = "dearseq",
    rawPval    = c(0.001, 0.2, 0.3)
  )

  out <- compute_pct_significant_by_comparison(baseline_df, pval_col = "rawPval", alpha = 0.01)

  expect_equal(out$n_significant, 1)
  expect_equal(out$n_evaluated, 3)
})

test_that("compute_pct_significant_by_comparison() errors when pval_col is missing", {
  baseline_df <- tibble::tibble(condition = "V1", time = 1, method = "dearseq")
  expect_error(compute_pct_significant_by_comparison(baseline_df))
})

test_that("plot_baseline_significance_comparison() returns a ggplot with bars matching pct_significant", {
  testthat::skip_if_not_installed("ggplot2")

  pct_df <- tibble::tibble(
    condition        = c("Vaccine A", "Vaccine A", "Vaccine B", "Vaccine B"),
    time              = c(1, 1, 1, 1),
    method            = c("dearseq", "qusage", "dearseq", "qusage"),
    n_significant      = c(10, 20, 5, 5),
    n_evaluated         = c(100, 100, 100, 100),
    pct_significant      = c(10, 20, 5, 5)
  )

  p <- plot_baseline_significance_comparison(
    pct_df, conditions_order = c("Vaccine A", "Vaccine B")
  )

  expect_s3_class(p, "ggplot")
  expect_equal(sort(p$data$pct_significant), sort(pct_df$pct_significant))
})

test_that("plot_baseline_significance_comparison() can suppress value labels", {
  testthat::skip_if_not_installed("ggplot2")

  pct_df <- tibble::tibble(
    condition = "Vaccine A", time = 1, method = c("dearseq", "qusage"),
    n_significant = c(1, 2), n_evaluated = c(10, 10), pct_significant = c(10, 20)
  )

  p_with_labels    <- plot_baseline_significance_comparison(pct_df, conditions_order = "Vaccine A", value_labels = TRUE)
  p_without_labels <- plot_baseline_significance_comparison(pct_df, conditions_order = "Vaccine A", value_labels = FALSE)

  expect_equal(length(p_with_labels$layers), length(p_without_labels$layers) + 1)
})

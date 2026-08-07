test_that("join_robustness_baseline() keeps one row per gene set x comparison x method", {
  robustness_df <- tibble::tibble(
    gs.name   = c("gs1", "gs1", "gs2"),
    condition = c("V1", "V1", "V1"),
    time      = c(1, 1, 1),
    robustness = c(0.2, 0.2, 0.8)
  )
  baseline_df <- tibble::tibble(
    gs.name              = c("gs1", "gs1", "gs2"),
    condition            = c("V1", "V1", "V1"),
    time                 = c(1, 1, 1),
    method               = c("dearseq", "qusage", "dearseq"),
    global.adjPval_BH    = c(0.01, 0.2, 0.9)
  )

  out <- join_robustness_baseline(robustness_df, baseline_df)

  expect_equal(nrow(out), 3)
  expect_true(all(c("robustness", "method", "global.adjPval_BH") %in% colnames(out)))
  expect_equal(out$robustness[out$gs.name == "gs1" & out$method == "dearseq"], 0.2)
  expect_equal(out$robustness[out$gs.name == "gs1" & out$method == "qusage"], 0.2)
})

test_that("join_robustness_baseline() drops rows with no robustness match", {
  robustness_df <- tibble::tibble(
    gs.name = "gs1", condition = "V1", time = 1, robustness = 0.5
  )
  baseline_df <- tibble::tibble(
    gs.name = c("gs1", "gs2"), condition = c("V1", "V1"), time = c(1, 1),
    method = c("dearseq", "dearseq"), global.adjPval_BH = c(0.01, 0.02)
  )

  out <- join_robustness_baseline(robustness_df, baseline_df)

  expect_equal(nrow(out), 1)
  expect_equal(out$gs.name, "gs1")
})

test_that("plot_robustness_vs_baseline() binary mode returns a ggplot and classifies significance correctly", {
  testthat::skip_if_not_installed("ggplot2")

  df <- tibble::tibble(
    gs.name               = c("gs1", "gs2", "gs3", "gs4"),
    method                = c("dearseq", "dearseq", "qusage", "qusage"),
    robustness             = c(0.9, 0.1, 0.8, 0.2),
    withinTime.adjPval_BH = c(0.001, 0.5, 0.02, 0.9)
  )

  p <- plot_robustness_vs_baseline(df, mode = "binary")

  expect_s3_class(p, "ggplot")
  expect_equal(
    as.character(p$data$baseline_significant),
    c("Significant", "Not significant", "Significant", "Not significant")
  )
})

test_that("plot_robustness_vs_baseline() continuous mode returns a ggplot with a -log10(p) column", {
  testthat::skip_if_not_installed("ggplot2")

  df <- tibble::tibble(
    gs.name               = c("gs1", "gs2", "gs3"),
    method                = c("dearseq", "dearseq", "qusage"),
    robustness             = c(0.9, 0.1, 0.5),
    withinTime.adjPval_BH = c(0.001, 0.5, 0)
  )

  p <- plot_robustness_vs_baseline(df, mode = "continuous")

  expect_s3_class(p, "ggplot")
  expect_true(all(is.finite(p$data$neg_log10_p)))
  expect_equal(p$data$neg_log10_p[1], -log10(0.001))
})

test_that("plot_robustness_vs_baseline() errors informatively when pval_col is missing", {
  df <- tibble::tibble(gs.name = "gs1", method = "dearseq", robustness = 0.5)
  expect_error(plot_robustness_vs_baseline(df, mode = "binary", pval_col = "not_a_column"))
})

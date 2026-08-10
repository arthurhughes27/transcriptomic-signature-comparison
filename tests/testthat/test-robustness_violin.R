test_that("plot_robustness_violin() returns a ggplot without erroring", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("ggh4x")

  set.seed(1)
  robustness_df <- tidyr::expand_grid(
    gs.name   = paste0("gs", 1:10),
    condition = c("Vaccine A", "Vaccine B"),
    time      = c(1, 7)
  ) |>
    dplyr::mutate(robustness = stats::runif(dplyr::n()))

  p <- plot_robustness_violin(robustness_df, conditions_order = c("Vaccine A", "Vaccine B"))

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$y, "Signal-robustness")
  expect_equal(p$labels$title, "Distribution of signal-robustness across gene sets, by comparison")
})

test_that("plot_robustness_violin() orders comparisons the same way as the heatmaps", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("ggh4x")

  robustness_df <- tidyr::expand_grid(
    gs.name   = paste0("gs", 1:3),
    condition = c("V2", "V1"),
    time      = c(7, 1)
  ) |>
    dplyr::mutate(robustness = 0.5)

  p <- plot_robustness_violin(robustness_df, conditions_order = c("V1", "V2"))

  expect_equal(levels(p$data$condition), c("V1", "V2"))
  expect_equal(levels(p$data$time), c("1", "7"))
})

test_that("plot_robustness_violin() appends unrecognised conditions with a warning", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("ggh4x")

  robustness_df <- tibble::tibble(
    gs.name = c("gs1", "gs1"), condition = c("V1", "V3"), time = c(1, 1),
    robustness = c(0.4, 0.6)
  )

  expect_warning(
    p <- plot_robustness_violin(robustness_df, conditions_order = c("V1")),
    "V3"
  )
  expect_equal(levels(p$data$condition), c("V1", "V3"))
})

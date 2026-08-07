test_that("abbreviate_geneset_label() applies the expected substitutions", {
  expect_equal(abbreviate_geneset_label("T cells activation network"), "T-cells act. netw.")
  expect_equal(abbreviate_geneset_label("antigen presentation"), "antigen pres.")
  expect_equal(abbreviate_geneset_label("no substitutions here"), "no substitutions here")
})

test_that("determine_heatmap_significance() filter_mode = 'none' uses the p-value only", {
  df <- tibble::tibble(
    condition          = "V1", time = 1,
    fc.score            = c(2, 0.1, -3),
    global.adjPval_BH = c(0.01, 0.01, 0.5)
  )

  out <- determine_heatmap_significance(
    df, p_approach = "global", p_correction = "BH", p_threshold = 0.05, filter_mode = "none"
  )

  expect_equal(as.logical(as.character(out$significant)), c(TRUE, TRUE, FALSE))
})

test_that("determine_heatmap_significance() filter_mode = 'user' requires p-value AND score above threshold", {
  df <- tibble::tibble(
    condition          = "V1", time = 1,
    fc.score            = c(2, 0.1, 3),
    global.adjPval_BH = c(0.01, 0.01, 0.5)
  )

  out <- determine_heatmap_significance(
    df, p_approach = "global", p_correction = "BH", p_threshold = 0.05,
    filter_mode = "user", score_col = "fc.score", user_threshold = 0.5
  )

  expect_equal(as.logical(as.character(out$significant)), c(TRUE, FALSE, FALSE))
})

test_that("determine_heatmap_significance() filter_mode = 'data' withinTime computes per-time quantile thresholds", {
  df <- tibble::tibble(
    condition          = "V1",
    time                = c(1, 1, 7, 7),
    fc.score             = c(1, 5, 1, 5),
    withinTime.adjPval_BH = 0.01
  )

  out <- determine_heatmap_significance(
    df, p_approach = "withinTime", p_correction = "BH", p_threshold = 0.05,
    filter_mode = "data", score_col = "fc.score", quantile_threshold = 0.5
  )

  # Within each time, only the higher-scoring row exceeds the median.
  expect_equal(as.logical(as.character(out$significant)), c(FALSE, TRUE, FALSE, TRUE))
})

test_that("select_heatmap_comparisons() filters and builds time_label", {
  df <- tibble::tibble(
    gs.aggregate = c("A", "A", "B"),
    condition     = c("V1", "V1", "V2"),
    time           = c(1, 7, 1)
  )

  out <- select_heatmap_comparisons(df, conditions = "V1", times = c(1, 7), aggregates = "A")

  expect_equal(nrow(out), 2)
  expect_equal(as.character(out$time_label), c("V1 - Day 1", "V1 - Day 7"))
})

test_that("filter_common_de() 'none' returns df unchanged", {
  df <- tibble::tibble(gs.name = c("gs1", "gs2"))
  expect_equal(filter_common_de(df, full_df = df, method_name = "dearseq", filter_commonDE = "none"), df)
})

test_that("filter_common_de() 'global' keeps gene sets significant in enough of the selected rows", {
  df <- tibble::tibble(
    gs.name     = c("gs1", "gs1", "gs2", "gs2"),
    significant = factor(c(TRUE, TRUE, TRUE, FALSE), levels = c(TRUE, FALSE))
  )

  out <- filter_common_de(df, full_df = df, method_name = "dearseq", filter_commonDE = "global", common_proportion = 1)

  expect_equal(unique(out$gs.name), "gs1")
})

test_that("filter_common_de() 'score' keeps gene sets with enough same-direction significant vaccines", {
  full_df <- tibble::tibble(
    gs.name           = c("gs1", "gs1", "gs1", "gs2", "gs2"),
    condition          = c("V1", "V2", "V3", "V1", "V2"),
    method              = "dearseq",
    fc.score             = c(1, 1, 1, 1, -1),
    global.adjPval_BH  = c(0.01, 0.01, 0.5, 0.01, 0.01)
  )
  df <- dplyr::filter(full_df, condition %in% c("V1", "V2", "V3"))

  out <- filter_common_de(
    df, full_df = full_df, method_name = "dearseq", filter_commonDE = "score",
    score_threshold = 2, adj_pval_col = "global.adjPval_BH", score_col = "fc.score", p_threshold = 0.05
  )

  # gs1: 2 vaccines significant "up" (V1, V2) -> sharing score 2, kept.
  # gs2: 1 vaccine "up" and 1 "down", each with sharing score 1 -> dropped.
  expect_equal(unique(out$gs.name), "gs1")
})

test_that("clip_scores() clips to the given quantile of absolute value", {
  df <- tibble::tibble(fc.score = c(-10, -1, 0, 1, 10))
  expected_threshold <- stats::quantile(abs(df$fc.score), 0.6)

  out <- clip_scores(df, score_col = "fc.score", quantile_scoreclip = 0.6)

  expect_equal(out$fc.score, pmin(pmax(df$fc.score, -expected_threshold), expected_threshold), ignore_attr = TRUE)
  expect_true(max(abs(out$fc.score)) < max(abs(df$fc.score)))
})

synthetic_heatmap_results <- function() {
  tibble::tibble(
    gs.name              = rep(c("gs1", "gs2"), each = 4),
    gs.name.description  = rep(c("Gene set one", "Gene set two"), each = 4),
    gs.aggregate          = factor(rep(c("A", "B"), each = 4), levels = c("A", "B")),
    gs.colour              = rep(c("#111111", "#222222"), each = 4),
    condition                = rep(c("V1", "V2"), times = 4),
    condition.colour          = rep(c("#aaaaaa", "#bbbbbb"), times = 4),
    time                        = rep(c(1, 1, 7, 7), times = 2),
    method                       = "dearseq",
    fc.score                      = c(1, -1, 2, -2, 0.5, -0.5, 1.5, -1.5),
    global.adjPval_BH            = 0.01
  )
}

test_that("plot_dgsa_heatmap() and save_stacked_heatmaps_pdf() build and save a real heatmap", {
  testthat::skip_if_not_installed("ComplexHeatmap")
  testthat::skip_if_not_installed("circlize")

  results_df <- synthetic_heatmap_results()

  ht <- plot_dgsa_heatmap(
    results_df,
    method_name = "dearseq",
    conditions  = c("V1", "V2"),
    times       = c(1, 7),
    aggregates  = c("A", "B"),
    p_correction = "BH", p_approach = "global", p_threshold = 0.05,
    filter_mode = "none", filter_commonDE = "none"
  )

  expect_s4_class(ht, "Heatmap")

  tmp_pdf <- tempfile(fileext = ".pdf")
  on.exit(unlink(tmp_pdf), add = TRUE)

  save_stacked_heatmaps_pdf(
    path = tmp_pdf, heatmaps = list(ht), titles = "dearseq",
    heights = list(grid::unit(1, "npc")), width = 6, height = 4
  )

  expect_true(fs::file_exists(tmp_pdf))
  expect_gt(fs::file_size(tmp_pdf), 0)
})

test_that("plot_dgsa_heatmap() errors informatively when no gene sets match", {
  testthat::skip_if_not_installed("ComplexHeatmap")

  results_df <- synthetic_heatmap_results()

  expect_error(
    plot_dgsa_heatmap(
      results_df,
      method_name = "dearseq",
      conditions  = "V1",
      times       = 999,
      aggregates  = c("A", "B")
    ),
    "No gene sets to plot"
  )
})

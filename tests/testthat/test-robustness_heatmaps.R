test_that("join_geneset_aggregates() attaches gs.aggregate and gs.label by gs.name", {
  robustness_df <- tibble::tibble(gs.name = c("gs1", "gs2", "gs3"), robustness = c(0.1, 0.5, 0.9))
  genesets <- list(
    geneset.names              = c("gs1", "gs2", "gs3"),
    geneset.aggregates          = factor(c("A", "B", "A"), levels = c("A", "B", "C")),
    geneset.names.descriptions = c("Gene set one", "Gene set two", "Gene set three")
  )

  out <- join_geneset_aggregates(robustness_df, genesets)

  expect_true(all(c("gs.aggregate", "gs.label") %in% colnames(out)))
  expect_equal(as.character(out$gs.aggregate), c("A", "B", "A"))
  expect_equal(levels(out$gs.aggregate), c("A", "B", "C"))
  expect_equal(out$gs.label, c("Gene set one", "Gene set two", "Gene set three"))
})

test_that("order_robustness_comparisons() orders condition/time and appends unrecognised conditions with a warning", {
  robustness_df <- tibble::tibble(
    condition = c("V2", "V1", "V3"),
    time      = c("7", "1", "3")
  )

  expect_warning(
    out <- order_robustness_comparisons(robustness_df, conditions_order = c("V1", "V2")),
    "V3"
  )

  expect_equal(levels(out$condition), c("V1", "V2", "V3"))
  expect_equal(levels(out$time), c("1", "3", "7"))
})

test_that("order_robustness_comparisons() does not warn when all conditions are recognised", {
  robustness_df <- tibble::tibble(condition = c("V2", "V1"), time = c(3, 1))
  expect_no_warning(order_robustness_comparisons(robustness_df, conditions_order = c("V1", "V2")))
})

test_that("summarise_robustness_by_aggregate() computes the unweighted mean per aggregate", {
  robustness_df <- tibble::tibble(
    gs.aggregate = factor(c("A", "A", "B"), levels = c("A", "B")),
    condition    = "V1",
    time         = 7,
    robustness   = c(0.2, 0.8, 0.5)
  )

  out <- summarise_robustness_by_aggregate(robustness_df)

  expect_equal(out$mean_robustness[out$gs.aggregate == "A"], 0.5)
  expect_equal(out$mean_robustness[out$gs.aggregate == "B"], 0.5)
  expect_equal(out$n_gene_sets[out$gs.aggregate == "A"], 2)
})

test_that("plot_robustness_heatmap_aggregate() returns a ggplot without erroring", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("ggh4x")

  set.seed(1)
  robustness_df <- tidyr::expand_grid(
    gs.name   = c("gs1", "gs2"),
    condition = c("Vaccine A", "Vaccine B"),
    time      = c(1, 7)
  ) |>
    dplyr::mutate(
      gs.aggregate = factor(ifelse(gs.name == "gs1", "A", "B"), levels = c("A", "B")),
      robustness    = stats::runif(dplyr::n())
    )

  p <- plot_robustness_heatmap_aggregate(robustness_df, conditions_order = c("Vaccine A", "Vaccine B"))

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Mean signal-robustness by gene-set aggregate")
})

test_that("plot_robustness_heatmap_genesets() returns a single-page list by default and drops always-null gene sets", {
  testthat::skip_if_not_installed("patchwork")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("ggh4x")

  robustness_df <- tidyr::expand_grid(
    gs.name   = c("gs1", "gs2", "gs3"),
    condition = c("Vaccine A", "Vaccine B"),
    time      = c(1, 7)
  ) |>
    dplyr::mutate(
      gs.aggregate = factor(dplyr::if_else(gs.name %in% c("gs1", "gs2"), "A", "B"), levels = c("A", "B")),
      gs.label      = gs.name,
      robustness     = ifelse(gs.name == "gs3", 0, 0.5)  # gs3 is always null
    )

  expect_message(
    pages <- plot_robustness_heatmap_genesets(
      robustness_df,
      conditions_order = c("Vaccine A", "Vaccine B"),
      aggregate_colors  = c("#111111", "#222222")
    ),
    "Dropping 1 gene set"
  )

  expect_type(pages, "list")
  expect_length(pages, 1)
  expect_s3_class(pages[[1]], "patchwork")
})

test_that("plot_robustness_heatmap_genesets() keeps all gene sets when drop_null_gene_sets = FALSE", {
  testthat::skip_if_not_installed("patchwork")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("ggh4x")

  robustness_df <- tidyr::expand_grid(
    gs.name   = c("gs1", "gs2"),
    condition = "Vaccine A",
    time      = 1
  ) |>
    dplyr::mutate(
      gs.aggregate = factor("A", levels = "A"),
      gs.label      = gs.name,
      robustness     = 0
    )

  pages <- plot_robustness_heatmap_genesets(
    robustness_df,
    conditions_order    = "Vaccine A",
    aggregate_colors     = "#111111",
    drop_null_gene_sets   = FALSE
  )

  expect_length(pages, 1)
  expect_s3_class(pages[[1]], "patchwork")
})

test_that("plot_robustness_heatmap_genesets() splits gene sets across pages when rows_per_page is exceeded", {
  testthat::skip_if_not_installed("patchwork")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("ggh4x")

  robustness_df <- tidyr::expand_grid(
    gs.name   = paste0("gs", 1:5),
    condition = "Vaccine A",
    time      = 1
  ) |>
    dplyr::mutate(
      gs.aggregate = factor("A", levels = "A"),
      gs.label      = gs.name,
      robustness     = 0.5
    )

  pages <- plot_robustness_heatmap_genesets(
    robustness_df,
    conditions_order = "Vaccine A",
    aggregate_colors  = "#111111",
    rows_per_page      = 2
  )

  expect_length(pages, 3)  # 5 gene sets, 2 per page -> 3 pages
  for (p in pages) expect_s3_class(p, "patchwork")
})

test_that("plot_robustness_heatmap_genesets() rows_per_page = Inf gives a single page", {
  testthat::skip_if_not_installed("patchwork")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("ggh4x")

  robustness_df <- tidyr::expand_grid(
    gs.name   = paste0("gs", 1:5),
    condition = "Vaccine A",
    time      = 1
  ) |>
    dplyr::mutate(
      gs.aggregate = factor("A", levels = "A"),
      gs.label      = gs.name,
      robustness     = 0.5
    )

  pages <- plot_robustness_heatmap_genesets(
    robustness_df,
    conditions_order = "Vaccine A",
    aggregate_colors  = "#111111",
    rows_per_page      = Inf
  )

  expect_length(pages, 1)
})

test_that("plot_robustness_heatmap_genesets() displays the full gs.label, not the short gs.name", {
  testthat::skip_if_not_installed("patchwork")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("ggh4x")

  robustness_df <- tidyr::expand_grid(
    gs.name   = c("gs1", "gs2"),
    condition = "Vaccine A",
    time      = 1
  ) |>
    dplyr::mutate(
      gs.aggregate = factor("A", levels = "A"),
      gs.label      = c("Full description one", "Full description two")[match(gs.name, c("gs1", "gs2"))],
      robustness     = 0.5
    )

  pages <- plot_robustness_heatmap_genesets(
    robustness_df, conditions_order = "Vaccine A", aggregate_colors = "#111111"
  )

  p_main <- pages[[1]][[2]]
  expect_setequal(as.character(p_main$data$gs.label), c("Full description one", "Full description two"))
  expect_false(any(c("gs1", "gs2") %in% as.character(p_main$data$gs.label)))
})

test_that("plot_robustness_heatmap_genesets() keeps the aggregate legend identical across pages", {
  testthat::skip_if_not_installed("patchwork")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("ggh4x")

  robustness_df <- tidyr::expand_grid(
    gs.name   = paste0("gs", 1:4),
    condition = "Vaccine A",
    time      = 1
  ) |>
    dplyr::mutate(
      gs.aggregate = factor(c("A", "A", "B", "B")[match(gs.name, paste0("gs", 1:4))], levels = c("A", "B")),
      gs.label      = gs.name,
      robustness     = 0.5
    )

  pages <- plot_robustness_heatmap_genesets(
    robustness_df,
    conditions_order = "Vaccine A", aggregate_colors = c("#111111", "#222222"),
    rows_per_page = 2
  )

  # Page 2 (gs3, gs4) only has aggregate "B" gene sets, but its legend
  # should still list both aggregates, matching page 1.
  strip_scale_page2 <- ggplot2::layer_scales(pages[[2]][[1]])$fill
  expect_equal(strip_scale_page2$get_limits(), c("A", "B"))
})

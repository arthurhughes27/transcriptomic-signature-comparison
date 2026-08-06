test_that("extract_dgsa_results() parses comparison names into conditions and times", {
  results_list <- list(
    "V1 vs Control - Day 1" = list(
      pvals = list(rawPval = c(0.01, 0.2)),
      score = list(activation.scores = list(0.5, -0.2), fc.scores = list(0.6, -0.3)),
      cor   = list(mean.corr = c(0.3, -0.1), corr.mean = c(0.25, -0.05))
    ),
    "V2 vs Control - Day 7" = list(
      pvals = list(rawPval = c(0.4, 0.02)),
      score = list(activation.scores = list(0.1, 0.2), fc.scores = list(0.15, 0.25)),
      cor   = list(mean.corr = c(0.1, 0.2), corr.mean = c(0.05, 0.15))
    )
  )

  out <- extract_dgsa_results(results_list)

  expect_equal(out$conditions, c("V1", "V2"))
  expect_equal(out$times, c(1, 7))
  expect_equal(out$raw_pvals[[1]], c(0.01, 0.2))
})

test_that("apply_pvalue_adjustments() produces correct adjusted p-values for every scope x method", {
  df <- tibble::tibble(
    comparison = c("V1 - Day 1", "V1 - Day 1", "V1 - Day 7", "V1 - Day 7"),
    time       = c(1, 1, 7, 7),
    rawPval    = c(0.01, 0.2, 0.03, 0.5)
  )

  out <- apply_pvalue_adjustments(df, methods = c("BH", "bonferroni"))

  expect_true(all(c(
    "global.adjPval_BH", "withinTime.adjPval_BH", "withinComparison.adjPval_BH",
    "global.adjPval_bonferroni", "withinTime.adjPval_bonferroni", "withinComparison.adjPval_bonferroni"
  ) %in% colnames(out)))

  expect_equal(out$global.adjPval_BH, stats::p.adjust(df$rawPval, method = "BH"))

  expected_within_time <- df |>
    dplyr::group_by(time) |>
    dplyr::mutate(adj = stats::p.adjust(rawPval, method = "BH")) |>
    dplyr::pull(adj)
  expect_equal(out$withinTime.adjPval_BH, expected_within_time)

  expected_within_comparison <- df |>
    dplyr::group_by(comparison) |>
    dplyr::mutate(adj = stats::p.adjust(rawPval, method = "bonferroni")) |>
    dplyr::pull(adj)
  expect_equal(out$withinComparison.adjPval_bonferroni, expected_within_comparison)
})

test_that("build_tidy_dgsa_results() tidies a results list, dropping NULL comparisons", {
  genesets <- list(
    geneset.names       = c("gs1", "gs2", "gs3"),
    geneset.descriptions = c("d1", "d2", "d3"),
    geneset.aggregates    = factor(c("A", "A", "B"))
  )

  results_list <- list(
    "V1 vs Control - Day 1" = list(
      pvals = list(rawPval = c(0.01, 0.2, 0.5)),
      score = list(activation.scores = list(0.5, -0.2, 0.1), fc.scores = list(0.6, -0.3, 0.05)),
      cor   = list(mean.corr = c(0.3, -0.1, 0.2), corr.mean = c(0.25, -0.05, 0.15))
    ),
    "V1 vs Control - Day 7" = list(
      pvals = list(rawPval = c(0.001, 0.3, 0.04)),
      score = list(activation.scores = list(0.9, 0.05, -0.4), fc.scores = list(1.0, 0.02, -0.5)),
      cor   = list(mean.corr = c(0.4, 0.0, -0.2), corr.mean = c(0.35, 0.01, -0.25))
    ),
    "V2 vs Control - Day 1" = NULL
  )

  out <- build_tidy_dgsa_results(results_list, genesets, method = "dearseq")

  # 2 completed comparisons x 3 gene sets = 6 rows
  expect_equal(nrow(out), 6)
  expect_setequal(unique(out$comparison), c("V1 - Day 1", "V1 - Day 7"))
  expect_true(all(c("global.adjPval_BH", "withinTime.adjPval_holm") %in% colnames(out)))
  expect_equal(unique(out$method), "dearseq")
})

test_that("build_tidy_dgsa_results() errors when geneset counts mismatch", {
  genesets <- list(
    geneset.names       = c("gs1", "gs2"),
    geneset.descriptions = c("d1", "d2"),
    geneset.aggregates    = factor(c("A", "B"))
  )
  results_list <- list(
    "V1 vs Control - Day 1" = list(
      pvals = list(rawPval = c(0.01, 0.2, 0.5)),  # 3 p-values, but genesets has 2
      score = list(activation.scores = list(0.5, -0.2, 0.1), fc.scores = list(0.6, -0.3, 0.05)),
      cor   = list(mean.corr = c(0.3, -0.1, 0.2), corr.mean = c(0.25, -0.05, 0.15))
    )
  )
  expect_error(build_tidy_dgsa_results(results_list, genesets))
})

test_that("build_tidy_dgsa_results() errors when all comparisons are NULL", {
  results_list <- list("V1 vs Control - Day 1" = NULL)
  genesets <- list(geneset.names = "gs1", geneset.descriptions = "d1", geneset.aggregates = factor("A"))
  expect_error(build_tidy_dgsa_results(results_list, genesets))
})

test_that("build_tidy_dgsa_results() orders conditions and warns about unlisted ones", {
  genesets <- list(
    geneset.names       = "gs1",
    geneset.descriptions = "d1",
    geneset.aggregates    = factor("A")
  )
  results_list <- list(
    "V2 vs Control - Day 1" = list(
      pvals = list(rawPval = 0.1),
      score = list(activation.scores = list(0.1), fc.scores = list(0.1)),
      cor   = list(mean.corr = 0.1, corr.mean = 0.1)
    ),
    "V1 vs Control - Day 1" = list(
      pvals = list(rawPval = 0.2),
      score = list(activation.scores = list(0.2), fc.scores = list(0.2)),
      cor   = list(mean.corr = 0.2, corr.mean = 0.2)
    )
  )

  expect_warning(
    out <- build_tidy_dgsa_results(results_list, genesets, conditions_order = c("V1")),
    "V2"
  )
  expect_equal(levels(out$condition), c("V1", "V2"))
})

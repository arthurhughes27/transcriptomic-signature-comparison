make_synthetic_raw_results_list <- function(seed) {
  set.seed(seed)
  list(
    "V1 vs Control - Day 1" = list(
      pvals = list(rawPval = c(0.01, 0.2)),
      score = list(activation.scores = c(1.1, 0.2), fc.scores = c(0.5, 0.1)),
      cor   = list(mean.corr = c(0.1, 0.2), corr.mean = c(0.3, 0.4))
    )
  )
}

test_that("load_baseline_results_from_raw() loads only is_baseline == TRUE rows, tidied and adjusted", {
  genesets <- list(
    geneset.names        = c("gs1", "gs2"),
    geneset.descriptions = c("desc1", "desc2"),
    geneset.aggregates    = factor(c("A", "B"))
  )

  raw_grid <- tibble::tibble(
    raw_spec_id = 1:3,
    method       = c("dearseq", "dearseq", "qusage"),
    spec_label   = c("dearseq_baseline", "dearseq_other", "qusage_baseline"),
    is_baseline  = c(TRUE, FALSE, TRUE)
  )

  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  saveRDS(make_synthetic_raw_results_list(1), fs::path(tmp_dir, "dearseq_baseline.rds"))
  saveRDS(make_synthetic_raw_results_list(2), fs::path(tmp_dir, "qusage_baseline.rds"))
  # dearseq_other.rds deliberately not written - should never be read.

  out <- load_baseline_results_from_raw(raw_grid, tmp_dir, genesets)

  expect_equal(sort(unique(out$method)), c("dearseq", "qusage"))
  expect_true("withinTime.adjPval_BH" %in% colnames(out))
  expect_equal(nrow(out), 4)  # 2 gene sets x 1 comparison x 2 baseline specs
})

test_that("load_baseline_results_from_raw() applies a shared conditions_order across specifications", {
  genesets <- list(
    geneset.names        = c("gs1", "gs2"),
    geneset.descriptions = c("desc1", "desc2"),
    geneset.aggregates    = factor(c("A", "B"))
  )

  raw_grid <- tibble::tibble(
    raw_spec_id = 1:2,
    method       = c("dearseq", "qusage"),
    spec_label   = c("dearseq_baseline", "qusage_baseline"),
    is_baseline  = c(TRUE, TRUE)
  )

  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  saveRDS(make_synthetic_raw_results_list(1), fs::path(tmp_dir, "dearseq_baseline.rds"))
  saveRDS(make_synthetic_raw_results_list(2), fs::path(tmp_dir, "qusage_baseline.rds"))

  out <- load_baseline_results_from_raw(raw_grid, tmp_dir, genesets, conditions_order = c("V1", "V0"))

  expect_equal(levels(out$condition), c("V1", "V0"))
})

test_that("load_baseline_results_from_raw() errors informatively when a baseline results file is missing", {
  raw_grid <- tibble::tibble(
    raw_spec_id = 1, method = "dearseq", spec_label = "dearseq_baseline", is_baseline = TRUE
  )
  genesets <- list(geneset.names = "gs1", geneset.descriptions = "d1", geneset.aggregates = factor("A"))

  expect_error(
    load_baseline_results_from_raw(raw_grid, tempfile(), genesets),
    "not found"
  )
})

test_that("load_baseline_results_from_raw() errors when raw_grid has no baseline rows", {
  raw_grid <- tibble::tibble(raw_spec_id = 1, method = "dearseq", spec_label = "x", is_baseline = FALSE)
  genesets <- list(geneset.names = "gs1", geneset.descriptions = "d1", geneset.aggregates = factor("A"))

  expect_error(load_baseline_results_from_raw(raw_grid, tempfile(), genesets), "is_baseline")
})

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

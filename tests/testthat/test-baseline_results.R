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

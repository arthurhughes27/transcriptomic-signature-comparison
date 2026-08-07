synthetic_genesets <- list(
  geneset.names             = c("gs1", "gs2", "gs3"),
  geneset.names.descriptions = c(
    "Gene set one description (gs1)",
    "Gene set two description (gs2)",
    "Gene set three description (gs3)"
  )
)

test_that("attach_geneset_full_names() joins the full label by gs.name", {
  df  <- tibble::tibble(gs.name = c("gs2", "gs1"))
  out <- attach_geneset_full_names(df, synthetic_genesets)

  expect_equal(out$gs.label, c("Gene set two description (gs2)", "Gene set one description (gs1)"))
})

test_that("build_top_robust_table() returns the n highest-robustness rows, most robust first", {
  robustness_df <- tibble::tibble(
    gs.name    = c("gs1", "gs2", "gs3"),
    condition  = "Vaccine A",
    time       = 1,
    robustness = c(0.2, 0.9, 0.5)
  )

  out <- build_top_robust_table(robustness_df, synthetic_genesets, n = 2)

  expect_equal(nrow(out), 2)
  expect_equal(out$Robustness, c(0.9, 0.5))
  expect_true(out$`Gene set`[1] == "Gene set two description (gs2)")
  expect_equal(colnames(out), c("Gene set", "Vaccine", "Timepoint", "Robustness"))
})

test_that("build_significant_nonrobust_table() filters to significant-but-non-robust rows only", {
  robustness_df <- tibble::tibble(
    gs.name    = c("gs1", "gs2", "gs3"),
    condition  = "Vaccine A",
    time       = 1,
    robustness = c(0.1, 0.9, 0.3)  # gs2 is robust, should be excluded regardless of significance
  )
  baseline_df <- tibble::tibble(
    gs.name           = c("gs1", "gs2", "gs3"),
    condition         = "Vaccine A",
    time              = 1,
    method                = "dearseq",
    withinTime.adjPval_BH = c(0.001, 0.001, 0.5)  # gs3 not significant, should be excluded
  )

  out <- build_significant_nonrobust_table(
    robustness_df, baseline_df, synthetic_genesets,
    robustness_threshold = 0.5
  )

  expect_equal(nrow(out), 1)
  expect_equal(out$`Gene set`, "Gene set one description (gs1)")
  expect_equal(out$Robustness, 0.1)
  expect_equal(out$`Baseline p-value`, formatC(0.001, format = "e", digits = 2))
  expect_equal(colnames(out), c("Gene set", "Vaccine", "Timepoint", "Method", "Baseline p-value", "Robustness"))
})

test_that("build_significant_nonrobust_table() orders least robust first and respects n", {
  robustness_df <- tibble::tibble(
    gs.name    = c("gs1", "gs2", "gs3"),
    condition  = "Vaccine A",
    time       = 1,
    robustness = c(0.4, 0.1, 0.3)
  )
  baseline_df <- tibble::tibble(
    gs.name           = c("gs1", "gs2", "gs3"),
    condition         = "Vaccine A",
    time              = 1,
    method                = "dearseq",
    withinTime.adjPval_BH = c(0.01, 0.01, 0.01)
  )

  out <- build_significant_nonrobust_table(
    robustness_df, baseline_df, synthetic_genesets,
    robustness_threshold = 0.5, n = 2
  )

  expect_equal(nrow(out), 2)
  expect_equal(out$Robustness, c(0.1, 0.3))
})

test_that("save_latex_table() writes a .tex file containing booktabs markup and an optional .csv", {
  testthat::skip_if_not_installed("knitr")

  df <- tibble::tibble(`Gene set` = "gs1", Robustness = 0.5)
  tmp_tex <- tempfile(fileext = ".tex")
  tmp_csv <- tempfile(fileext = ".csv")
  on.exit(unlink(c(tmp_tex, tmp_csv)), add = TRUE)

  save_latex_table(df, path_tex = tmp_tex, caption = "A caption", label = "tab:test", path_csv = tmp_csv)

  tex_lines <- readLines(tmp_tex)
  expect_true(any(grepl("\\\\toprule", tex_lines)))
  expect_true(any(grepl("tab:test", tex_lines)))
  expect_true(fs::file_exists(tmp_csv))
})

# Runs the full test suite for R/. Run from the project root:
#
#   Rscript tests/run_tests.R
#
# Individual tests skip gracefully if optional heavy dependencies (e.g.
# dearseq) are not installed.

if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("The 'testthat' package is required to run tests. Install it with install.packages('testthat').")
}
if (!requireNamespace("fs", quietly = TRUE)) {
  stop("The 'fs' package is required to run tests. Install it with install.packages('fs').")
}

testthat::test_dir(fs::path("tests", "testthat"), reporter = "summary")

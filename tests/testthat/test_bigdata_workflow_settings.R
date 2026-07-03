test_that("bigData workflow uses higher RUST-compatible read depth", {
  script_path <- testthat::test_path("..", "tests_and_demos", "bigData.R")

  expect_error(parse(script_path), NA)

  script <- readLines(script_path, warn = FALSE)
  expect_true(any(grepl("3000a_high_depth_different_gene_counts", script, fixed = TRUE)))
  expect_true(any(grepl("n_genes <- 3000", script, fixed = TRUE)))
  expect_true(any(grepl("interceptMean = 6", script, fixed = TRUE)))
  expect_false(any(grepl("interceptMean = 4", script, fixed = TRUE)))
  expect_true(any(grepl("interceptSD = 1.5", script, fixed = TRUE)))
  expect_true(any(grepl("betaSD = 0.2", script, fixed = TRUE)))
})

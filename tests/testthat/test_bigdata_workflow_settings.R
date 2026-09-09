test_that("bigData workflow uses higher RUST-compatible read depth", {
  script_path <- testthat::test_path("..", "tests_and_demos", "bigData.R")

  expect_error(parse(script_path), NA)

  script <- readLines(script_path, warn = FALSE)
  expect_true(any(grepl("COVSIM_BIGDATA_N_GENES", script, fixed = TRUE)))
  expect_true(any(grepl("COVSIM_BIGDATA_BIAS_PROFILE", script, fixed = TRUE)))
  expect_true(any(grepl(
    'allowed_bias_profiles <- c("start_codon", "stop_codon", "similar")',
    script, fixed = TRUE
  )))
  expect_true(any(grepl("bias = bias_profile", script, fixed = TRUE)))
  expect_false(any(grepl('bias = "all"', script, fixed = TRUE)))
  expect_true(any(grepl("simulation_settings.tsv", script, fixed = TRUE)))
  expect_true(any(grepl('unset = "3000"', script, fixed = TRUE)))
  expect_true(any(grepl("interceptMean = 6", script, fixed = TRUE)))
  expect_false(any(grepl("interceptMean = 4", script, fixed = TRUE)))
  expect_true(any(grepl("interceptSD = 1.5", script, fixed = TRUE)))
  expect_true(any(grepl("betaSD = 0.2", script, fixed = TRUE)))
  expect_true(any(grepl("trailer_length <- rep(120, n_genes)", script, fixed = TRUE)))
  expect_true(any(grepl("max_uorfs = 2", script, fixed = TRUE)))
  expect_true(any(grepl('site_reference = "a_site"', script, fixed = TRUE)))
  expect_true(any(grepl("ground_truth = TRUE", script, fixed = TRUE)))
  expect_true(any(grepl("simulated_region_counts.tsv", script, fixed = TRUE)))
})

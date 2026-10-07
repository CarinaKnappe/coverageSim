test_that("the dataset simulation workflow uses a RUST-compatible read depth", {
  script_path <- testthat::test_path("..", "..", "analysis", "workflows", "simulate_dataset.R")

  expect_error(parse(script_path), NA)

  script <- readLines(script_path, warn = FALSE)
  expect_true(any(grepl("COVSIM_BIGDATA_N_GENES", script, fixed = TRUE)))
  expect_true(any(grepl("COVSIM_BIGDATA_BIAS_PROFILE", script, fixed = TRUE)))
  expect_true(any(grepl("COVSIM_BIGDATA_DMN_ALPHA_SCALE", script, fixed = TRUE)))
  expect_true(any(grepl(
    'allowed_bias_profiles <- c("median", "start_codon", "stop_codon", "similar")',
    script, fixed = TRUE
  )))
  expect_true(any(grepl('unset = "median"', script, fixed = TRUE)))
  expect_true(any(grepl("bias = bias_profile", script, fixed = TRUE)))
  expect_true(any(grepl("dmn_alpha_scale = dmn_alpha_scale", script, fixed = TRUE)))
  expect_false(any(grepl('bias = "all"', script, fixed = TRUE)))
  expect_true(any(grepl("simulation_settings.tsv", script, fixed = TRUE)))
  expect_true(any(grepl('unset = "3000"', script, fixed = TRUE)))
  # Read depth and fragment lengths are settable per run, but their defaults are
  # part of what this workflow promises, so both the mechanism and the default
  # are pinned here. A depth of 4 was too sparse once and must not come back.
  expect_true(any(grepl("COVSIM_BIGDATA_INTERCEPT_MEAN", script, fixed = TRUE)))
  expect_true(any(grepl("interceptMean = intercept_mean", script, fixed = TRUE)))
  expect_true(any(grepl('unset = "6"', script, fixed = TRUE)))
  expect_false(any(grepl('unset = "4"', script, fixed = TRUE)))
  expect_true(any(grepl("COVSIM_BIGDATA_READ_LENGTHS", script, fixed = TRUE)))
  expect_true(any(grepl("read_lengths_per = list(RFP = read_lengths)", script, fixed = TRUE)))
  expect_true(any(grepl('unset = "26:33"', script, fixed = TRUE)))
  # Both settings have to reach the run's own record, or a stored run cannot be
  # told apart from one made with different values.
  expect_true(any(grepl('"read_lengths", "intercept_mean"', script, fixed = TRUE)))
  expect_true(any(grepl("interceptSD = 1.5", script, fixed = TRUE)))
  expect_true(any(grepl("betaSD = 0.2", script, fixed = TRUE)))
  expect_true(any(grepl("trailer_length <- rep(120, n_genes)", script, fixed = TRUE)))
  expect_true(any(grepl("max_uorfs = 2", script, fixed = TRUE)))
  expect_true(any(grepl('site_reference = "a_site"', script, fixed = TRUE)))
  expect_true(any(grepl("ground_truth = TRUE", script, fixed = TRUE)))
  expect_true(any(grepl("simulated_region_counts.tsv", script, fixed = TRUE)))
})

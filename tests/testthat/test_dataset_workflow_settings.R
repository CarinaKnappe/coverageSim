test_that("the dataset simulation workflow uses a RUST-compatible read depth", {
  script_path <- testthat::test_path("..", "..", "analysis", "workflows", "simulate_dataset.R")

  expect_error(parse(script_path), NA)

  script <- readLines(script_path, warn = FALSE)
  expect_true(any(grepl("COVSIM_DATASET_N_GENES", script, fixed = TRUE)))
  expect_true(any(grepl("COVSIM_DATASET_BIAS_PROFILE", script, fixed = TRUE)))
  expect_true(any(grepl("COVSIM_DATASET_DMN_ALPHA_SCALE", script, fixed = TRUE)))
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
  expect_true(any(grepl("COVSIM_DATASET_INTERCEPT_MEAN", script, fixed = TRUE)))
  expect_true(any(grepl("interceptMean = intercept_mean", script, fixed = TRUE)))
  expect_true(any(grepl('unset = "6"', script, fixed = TRUE)))
  # The end bias: read from a learned fit, scaled by one strength, and actually
  # handed to the simulation. Without the second half a renamed argument would
  # leave the setting without effect, and a biased run would differ from the
  # control in nothing at all.
  expect_true(any(grepl("COVSIM_DATASET_END_BIAS_FIT", script, fixed = TRUE)))
  expect_true(any(grepl("COVSIM_DATASET_END_BIAS_STRENGTH", script, fixed = TRUE)))
  expect_true(any(grepl("five_prime_bias = five_prime_bias", script, fixed = TRUE)))
  expect_true(any(grepl("three_prime_bias = three_prime_bias", script, fixed = TRUE)))
  expect_true(any(grepl('"end_bias_fit", "end_bias_strength"', script, fixed = TRUE)))
  expect_false(any(grepl('unset = "4"', script, fixed = TRUE)))
  expect_true(any(grepl("COVSIM_DATASET_READ_LENGTHS", script, fixed = TRUE)))
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

test_that("the dataset workflow hands the end-bias fit and strength to the simulation", {
  # Run the workflow's own end_bias_settings() rather than searching its text,
  # so a strength that is read but not passed on fails here. Only that one
  # definition is evaluated; sourcing the whole script would start a simulation.
  script_path <- test_path("..", "..", "analysis", "workflows", "simulate_dataset.R")
  exprs <- parse(script_path)
  is_definition <- vapply(exprs, function(e) {
    is.call(e) && identical(e[[1]], as.name("<-")) &&
      identical(e[[2]], as.name("end_bias_settings"))
  }, logical(1))
  expect_equal(sum(is_definition), 1L)
  env <- new.env()
  eval(exprs[[which(is_definition)]], env)
  end_bias_settings <- env$end_bias_settings

  # No fit: no end bias at either end, whatever the strength.
  none <- end_bias_settings("", "3")
  expect_equal(none$five_prime_bias, list(source = "none"))
  expect_equal(none$three_prime_bias, list(source = "none"))

  table_5 <- data.table::data.table(kmer = c("A", "C", "G", "T"),
                                    weight = c(1, 1, 0.5, 2))
  table_3 <- data.table::data.table(kmer = c("A", "C", "G", "T"),
                                    weight = c(0.5, 1, 2, 1))
  fit <- structure(list(
    five_prime_bias = list(source = "learned", table = table_5, strength = 1),
    three_prime_bias = list(source = "learned", table = table_3, strength = 1)
  ), class = c("covsim_end_bias_fit", "list"))
  fit_path <- tempfile(fileext = ".rds")
  saveRDS(fit, fit_path)

  for (strength in c(0, 1, 3)) {
    settings <- end_bias_settings(fit_path, as.character(strength))
    expect_equal(settings$five_prime_bias$table, table_5)
    expect_equal(settings$three_prime_bias$table, table_3)
    expect_equal(settings$five_prime_bias$strength, strength)
    expect_equal(settings$three_prime_bias$strength, strength)
    # And the simulator turns it into the expected weights: the measured
    # weight raised to the strength, at both ends.
    geometry <- normalize_fragment_geometry(list(
      five_prime_bias = settings$five_prime_bias,
      three_prime_bias = settings$three_prime_bias
    ))
    expect_equal(end_bias_weight(paste0("T", strrep("A", 27)), 28L,
                                 geometry$five_prime_bias, NULL), 2^strength)
    expect_equal(end_bias_weight(paste0(strrep("C", 27), "G"), 28L,
                                 NULL, geometry$three_prime_bias), 2^strength)
  }

  expect_error(end_bias_settings(fit_path, "-1"), "END_BIAS_STRENGTH")
  expect_error(end_bias_settings(fit_path, "strong"), "END_BIAS_STRENGTH")
  expect_error(end_bias_settings(tempfile(), "1"), "does not exist")
  wrong_class <- tempfile(fileext = ".rds")
  saveRDS(list(five_prime_bias = fit$five_prime_bias), wrong_class)
  expect_error(end_bias_settings(wrong_class, "1"), "learn_end_bias")
})

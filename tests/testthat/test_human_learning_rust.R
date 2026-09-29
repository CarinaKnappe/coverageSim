test_that("human-learning RUST workflow is syntactically valid and uses every scenario", {
  scripts <- c(
    testthat::test_path("..", "..", "analysis", "benchmarks", "rust", "run_human_learning_rust.py"),
    testthat::test_path("..", "..", "analysis", "benchmarks", "rust", "plot_human_learning_rust.py")
  )
  for (script in scripts) {
    status <- system2("python3", c("-m", "py_compile", shQuote(script)))
    expect_equal(status, 0L)
  }
  run_script <- readLines(scripts[1], warn = FALSE)
  settings <- data.table::fread(testthat::test_path(
    "..", "..", "..", "coverageSim_data", "runs", "2026-09-11_human-learning",
    "scenario_settings.tsv"
  ))
  for (scenario in settings$scenario) {
    expect_true(any(grepl(paste0('"', scenario, '"'), run_script, fixed = TRUE)))
  }
  expect_true(any(grepl('"28"', run_script, fixed = TRUE)))
  expect_true(any(grepl("offset15_only=True", run_script, fixed = TRUE)))
  plot_script <- readLines(scripts[2], warn = FALSE)
  expect_true(any(grepl("LEFT_YLIM = (-3.0, 2.0)", plot_script, fixed = TRUE)))
  expect_true(any(grepl("RIGHT_YLIM = (0.0, 1.0)", plot_script, fixed = TRUE)))
})

test_that("optional upstream RUST run covers all human-learning simulations", {
  path <- Sys.getenv("COVSIM_HUMAN_RUST_OUTPUT", "")
  skip_if(!nzchar(path), "Set COVSIM_HUMAN_RUST_OUTPUT for the upstream RUST run")
  settings <- data.table::fread(file.path(path, "scenario_settings.tsv"))
  summary <- data.table::fread(file.path(path, "rust_signal_summary.tsv"))
  expect_setequal(summary$scenario, settings$scenario)
  expect_setequal(summary$stratum, c("length28", "offset15_pooled"))
  expect_equal(nrow(summary), 2L * nrow(settings))
  expect_true(all(summary$bam_records >= summary$accepted_records))
  expect_true(all(summary$accepted_records >= summary$retained_records))
  expect_true(all(summary$retained_records > 0))
  plot_check <- jsonlite::fromJSON(file.path(path, "rust_plot_validation.json"))
  expect_equal(plot_check$profiles, 32)
  expect_true(plot_check$plotted_curves_match_upstream_numerical_output)
  expect_false(plot_check$upstream_numerical_files_modified)
  expect_equal(unlist(plot_check$axis_limits$log2_observed_expected), c(-3, 2))
  expect_equal(unlist(plot_check$axis_limits$divergence), c(0, 1))
  expect_true(plot_check$axis_limits_verified_on_all_panels)
  expect_equal(plot_check$primary_pdf_mode, "codon")
  expect_equal(plot_check$nucleotide_output_role, "separate diagnostic")
  pdf <- file.path(path, "human_learning_RUST_metafootprints.pdf")
  expect_true(file.exists(pdf))
  expect_gt(file.info(pdf)$size, 10000)
  nucleotide_pdf <- file.path(path, "human_learning_RUST_nucleotide_diagnostics.pdf")
  expect_true(file.exists(nucleotide_pdf))
  expect_gt(file.info(nucleotide_pdf)$size, 10000)
})

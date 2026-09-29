# The analysis code lives in analysis/, not in the package's own R/ directory, so
# that loading coverageSim cannot pull in analysis adapters and so that the built
# package stays limited to the simulator itself.

test_that("analysis helpers are kept outside the package R directory", {
  expect_false(file.exists(test_path("..", "..", "R", "Learn_from_input.R")))
  expect_false(file.exists(test_path("..", "..", "R", "G_end_dropout_helpers.R")))

  helper_files <- c(
    test_path("..", "..", "analysis", "R", "Learn_from_input.R"),
    test_path("..", "..", "analysis", "R", "G_end_dropout_helpers.R")
  )

  expect_true(all(file.exists(helper_files)))
  expect_error(parse(helper_files[1]), NA)
  expect_error(parse(helper_files[2]), NA)
})

test_that("the run-producing workflows source their helpers explicitly", {
  workflow_dir <- test_path("..", "..", "analysis", "workflows")

  real_input_workflow <- file.path(
    workflow_dir, "Workflow coverageSim learn from real BAM.R"
  )
  dropout_workflow <- file.path(workflow_dir, "create_G_end_dropout_datasets.R")
  genome_only_workflow <- file.path(
    workflow_dir, "Workflow coverageSim human genome only.R"
  )

  expect_true(file.exists(real_input_workflow))
  expect_true(file.exists(dropout_workflow))
  expect_true(file.exists(genome_only_workflow))

  expect_error(parse(real_input_workflow), NA)
  expect_error(parse(dropout_workflow), NA)
  expect_error(parse(genome_only_workflow), NA)

  real_input_script <- readLines(real_input_workflow)
  dropout_script <- readLines(dropout_workflow)
  genome_only_script <- readLines(genome_only_workflow)

  # Sourced by path from analysis/R, never loaded as part of the package.
  expect_true(any(grepl("analysis.*Learn_from_input\\.R", real_input_script)))
  expect_true(any(grepl("analysis.*G_end_dropout_helpers\\.R", dropout_script)))

  # load_all(repo_dir) rather than load_all("."), so a workflow does not depend
  # on which directory it happens to be started from.
  for (script in list(real_input_script, dropout_script, genome_only_script)) {
    expect_true(any(grepl("devtools::load_all(repo_dir)", script, fixed = TRUE)))
    expect_false(any(grepl('devtools::load_all(".")', script, fixed = TRUE)))
  }

  expect_true(any(grepl("simCountTables\\(", genome_only_script)))
  expect_true(any(grepl("simNGScoverage\\(", genome_only_script)))
  expect_true(any(grepl("loadRegion(txdb_file, \"cds\")", genome_only_script, fixed = TRUE)))
  expect_false(any(grepl("learn_", genome_only_script)))
  expect_false(any(grepl("fimport\\(", genome_only_script)))
})

test_that("the RUST benchmark scripts sit beside the vendored tool they drive", {
  rust_dir <- test_path("..", "..", "analysis", "benchmarks", "rust")
  expect_true(dir.exists(rust_dir))
  # natural_end_model.py is imported by name, so it has to be in this directory.
  expect_true(file.exists(file.path(rust_dir, "natural_end_model.py")))
  expect_true(file.exists(file.path(rust_dir, "rust_plot_compat.py")))
})

test_that("RUST helpers are kept outside the coverageSim package R directory", {
  expect_false(file.exists(test_path("..", "..", "R", "Learn_from_input.R")))
  expect_false(file.exists(test_path("..", "..", "R", "G_end_dropout_helpers.R")))

  helper_files <- c(
    test_path("..", "..", "rust_helpers", "R", "Learn_from_input.R"),
    test_path("..", "..", "rust_helpers", "R", "G_end_dropout_helpers.R")
  )

  expect_true(all(file.exists(helper_files)))
  expect_error(parse(helper_files[1]), NA)
  expect_error(parse(helper_files[2]), NA)
})

test_that("RUST helper workflows explicitly source their helper files", {
  workflow_dir <- test_path("..", "..", "rust_helpers", "workflows")

  real_input_workflow <- file.path(
    workflow_dir,
    "Workflow coverageSim learn from real BAM.R"
  )
  dropout_workflow <- file.path(workflow_dir, "create_G_end_dropout_datasets.R")
  genome_only_workflow <- file.path(
    workflow_dir,
    "Workflow coverageSim human genome only.R"
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

  expect_true(any(grepl("rust_helpers.*Learn_from_input\\.R", real_input_script)))
  expect_true(any(grepl("rust_helpers.*G_end_dropout_helpers\\.R", dropout_script)))

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

test_that("human genome-only workflow is parseable and does not require BAM learning", {
  repo_dir <- normalizePath(
    Sys.getenv("COVSIM_REPO", unset = testthat::test_path("..", "..")),
    mustWork = TRUE
  )
  workflow <- file.path(
    repo_dir,
    "rust_helpers",
    "workflows",
    "Workflow coverageSim human genome only.R"
  )

  expect_true(file.exists(workflow))
  expect_error(parse(workflow), NA)

  script <- readLines(workflow)

  expect_true(any(grepl("simCountTables\\(", script)))
  expect_true(any(grepl("simNGScoverage\\(", script)))
  expect_true(any(grepl("loadRegion(txdb_file, \"cds\")", script, fixed = TRUE)))

  expect_true(any(grepl("interceptMean = 4", script, fixed = TRUE)))
  expect_true(any(grepl("interceptSD = 1.5", script, fixed = TRUE)))
  expect_true(any(grepl("betaSD = 0.2", script, fixed = TRUE)))
  expect_true(any(grepl("read_lengths_per = list(RFP = 28:30)", script, fixed = TRUE)))
  expect_true(any(grepl("load_seq_bias(type = \"codon\", shift = \"a-site\", bias = \"all\")", script, fixed = TRUE)))

  expect_false(any(grepl("learn_", script)))
  expect_false(any(grepl("fimport\\(", script)))
  expect_false(any(grepl("read.experiment\\(", script)))
})

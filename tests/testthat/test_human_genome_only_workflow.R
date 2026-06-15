test_that("human genome-only workflow is parseable and does not require BAM learning", {
  workflow <- test_path(
    "..",
    "tests_and_demos",
    "real_human_riboseq",
    "sources",
    "Workflow coverageSim human genome only.R"
  )

  expect_true(file.exists(workflow))
  expect_error(parse(workflow), NA)

  script <- readLines(workflow)

  expect_true(any(grepl("simCountTables\\(", script)))
  expect_true(any(grepl("simNGScoverage\\(", script)))
  expect_true(any(grepl("loadRegion(txdb_file, \"cds\")", script, fixed = TRUE)))

  expect_false(any(grepl("learn_", script)))
  expect_false(any(grepl("fimport\\(", script)))
  expect_false(any(grepl("read.experiment\\(", script)))
})

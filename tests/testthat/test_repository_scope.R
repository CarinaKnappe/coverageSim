test_that("Git excludes local analysis material but retains simulator sources", {
  testthat::skip_if(Sys.which("git") == "", "Git is unavailable")
  rules <- testthat::test_path("..", "..", ".gitignore")
  repository <- tempfile("coverageSim-git-scope-")
  dir.create(repository)
  file.copy(rules, file.path(repository, ".gitignore"))
  system2("git", c("-C", shQuote(repository), "init", "--quiet"))

  excluded <- c(
    "choros/choros_utils.R",
    "rust_helpers/R/Learn_from_input.R",
    "tests/testthat/test_bias_overview_workbook.R",
    "tests/testthat/test_choros_scripts.R",
    "tests/testthat/test_rust_helpers_separation.R",
    "tests/testthat/test_deprecated_dataset_layout.R",
    "tests/tests_and_demos/Parameter_settings_per_run.ods",
    "tests/tests_and_demos/learn_from_HumanGenome.R"
  )
  included <- c(
    "R/Sim_reads_from_counts.R", "R/fragment_geometry.R",
    "tests/testthat/test_simulations.R",
    "tests/testthat/test_simulated_rpf_fragments.R",
    "tests/tests_and_demos/bigData.R", "man/simNGScoverage.Rd"
  )
  ignored <- system2(
    "git",
    c("-C", shQuote(repository), "check-ignore", "--no-index", "--stdin"),
    input = c(excluded, included), stdout = TRUE
  )
  expect_setequal(ignored, excluded)
})

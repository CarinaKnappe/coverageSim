test_that("Git excludes data and vendored code but retains every analysis source", {
  testthat::skip_if(Sys.which("git") == "", "Git is unavailable")
  rules <- testthat::test_path("..", "..", ".gitignore")
  repository <- tempfile("coverageSim-git-scope-")
  dir.create(repository)
  file.copy(rules, file.path(repository, ".gitignore"))
  system2("git", c("-C", shQuote(repository), "init", "--quiet"))

  # Excluded for one of three reasons: it is bulk data, it is somebody else's
  # code that we only reference, or it is personal material. Note that the
  # folders themselves are not excluded -- only what sits inside them.
  excluded <- c(
    "rust_helpers/runs/2026-09-09_300tx_codon-bias_run2/reads/RFP_benchmark.bam",
    "rust_helpers/upstream/RUST/RUST/codon.py",
    "rust_helpers/.venv-rust-original/lib/python3.12/site-packages/numpy/version.py",
    "rust_helpers/workflows/__pycache__/natural_end_model.cpython-312.pyc",
    "choros/results/summary.tsv",
    "tests/tests_and_demos/Parameter_settings_per_run.ods",
    "git_ignore/coverageSim_bias_overview.xlsx",
    ".agents/onboarding.md",
    "CoverageSim_manuscript.pdf",
    "tests/tests_and_demos/human_flavoured_riboseq/reads/RFP_WT_1.bam",
    "test_human.2bit",
    "test_human.2bit_Human.test_seed",
    "test_human.gtf",
    "BSgenome.Human.test.genc2522M/inst/extdata/single_sequences.2bit"
  )
  # Every source file is tracked, including the analysis code that lives beside
  # the package and the benchmark tests that skip without their external data.
  included <- c(
    "R/Sim_reads_from_counts.R", "R/fragment_geometry.R", "R/code_version.R",
    "man/simNGScoverage.Rd",
    "tests/testthat/test_simulations.R",
    "tests/testthat/test_simulated_rpf_fragments.R",
    "tests/testthat/test_choros_scripts.R",
    "tests/testthat/test_rust_helpers_separation.R",
    "tests/testthat/test_bias_overview_workbook.R",
    "tests/testthat/test_deprecated_dataset_layout.R",
    "tests/testthat/natural_end_model_checks.py",
    "tests/tests_and_demos/bigData.R",
    "tests/tests_and_demos/learn_from_HumanGenome.R",
    "tests/tests_and_demos/figure3c_alpha_calibration.R",
    "choros/choros_utils.R",
    "rust_helpers/R/Learn_from_input.R",
    "rust_helpers/workflows/natural_end_model.py",
    "rust_helpers/tests/test_learn_from_input.R"
  )
  ignored <- system2(
    "git",
    c("-C", shQuote(repository), "check-ignore", "--no-index", "--stdin"),
    input = c(excluded, included), stdout = TRUE
  )
  expect_setequal(ignored, excluded)
})

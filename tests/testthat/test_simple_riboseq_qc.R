test_that("simple Ribo-seq QC converts fragments to per-base coverage", {
  old <- Sys.getenv("COVSIM_SIMPLE_QC_AUTORUN", unset = NA_character_)
  on.exit({
    if (is.na(old)) Sys.unsetenv("COVSIM_SIMPLE_QC_AUTORUN") else
      Sys.setenv(COVSIM_SIMPLE_QC_AUTORUN = old)
  })
  Sys.setenv(COVSIM_SIMPLE_QC_AUTORUN = "false")
  script <- testthat::test_path("..", "tests_and_demos", "simple_riboseq_qc.R")
  env <- new.env(parent = globalenv())
  sys.source(script, envir = env)
  env$load_simulated_rpf_qc_helpers()

  expect_true(exists("map_genomic_to_transcript", envir = env, inherits = FALSE))

  truth <- data.table::data.table(
    tx_position = c(3L, 5L), site_offset = c(2L, 2L),
    fragment_length = c(4L, 3L), score = c(2L, 1L)
  )
  coverage <- env$fragment_coverage_vectors(truth, 8L)

  expect_equal(coverage$fragment, c(2, 2, 3, 3, 1, 0, 0, 0))
  expect_equal(coverage$a_site, c(0, 0, 2, 0, 1, 0, 0, 0))
})

test_that("simple Ribo-seq QC chooses covered examples on each strand", {
  old <- Sys.getenv("COVSIM_SIMPLE_QC_AUTORUN", unset = NA_character_)
  on.exit({
    if (is.na(old)) Sys.unsetenv("COVSIM_SIMPLE_QC_AUTORUN") else
      Sys.setenv(COVSIM_SIMPLE_QC_AUTORUN = old)
  })
  Sys.setenv(COVSIM_SIMPLE_QC_AUTORUN = "false")
  script <- testthat::test_path("..", "tests_and_demos", "simple_riboseq_qc.R")
  env <- new.env(parent = globalenv())
  sys.source(script, envir = env)

  truth <- data.table::data.table(
    transcript_id = c("plus_low", "plus_high", "minus"),
    strand = c("+", "+", "-"), score = c(2L, 10L, 7L),
    cigar = c("10M1N18M", "10M1N18M", "10M1N18M")
  )
  selected <- env$select_example_transcripts(truth)

  expect_setequal(selected$transcript_id, c("plus_high", "minus"))
})

test_that("simulated-RPF QC maps plus and minus transcript coordinates", {
  old <- Sys.getenv("COVSIM_QC_AUTORUN", unset = NA_character_)
  on.exit({
    if (is.na(old)) Sys.unsetenv("COVSIM_QC_AUTORUN") else
      Sys.setenv(COVSIM_QC_AUTORUN = old)
  })
  Sys.setenv(COVSIM_QC_AUTORUN = "false")
  script <- testthat::test_path("..", "tests_and_demos", "simulated_rpf_qc.R")
  env <- new.env(parent = globalenv())
  sys.source(script, envir = env)

  plus <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(c(101L, 201L), width = 10L), "+",
    exon_rank = 1:2
  )
  minus <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(c(401L, 301L), width = 10L), "-",
    exon_rank = 1:2
  )

  expect_equal(env$map_genomic_to_transcript(c(101L, 110L, 201L), plus, "+"),
               c(1L, 10L, 11L))
  expect_equal(env$map_genomic_to_transcript(c(410L, 401L, 310L), minus, "-"),
               c(1L, 10L, 11L))
})

test_that("simulated-RPF QC calculates frame and splicing summaries", {
  old <- Sys.getenv("COVSIM_QC_AUTORUN", unset = NA_character_)
  on.exit({
    if (is.na(old)) Sys.unsetenv("COVSIM_QC_AUTORUN") else
      Sys.setenv(COVSIM_QC_AUTORUN = old)
  })
  Sys.setenv(COVSIM_QC_AUTORUN = "false")
  script <- testthat::test_path("..", "tests_and_demos", "simulated_rpf_qc.R")
  env <- new.env(parent = globalenv())
  sys.source(script, envir = env)

  truth <- data.table::data.table(
    sample = "RFP_WT_1", transcript_id = "tx1", region = "cds",
    region_position = c(1L, 4L, 2L), cds_length = 9L,
    fragment_length = 28L, site_offset = 15L,
    strand = "+", cigar = c("28M", "10M90N18M", "28M"),
    sequence = c(paste(rep("A", 28), collapse = ""),
                 paste(rep("C", 28), collapse = ""),
                 paste(rep("G", 28), collapse = "")),
    region_length = 9L, score = c(3L, 2L, 1L)
  )
  dataset <- tempfile()
  dir.create(file.path(dataset, "reads"), recursive = TRUE)
  # Test summaries that do not require opening a BAM.
  frame <- truth[, .(reads = sum(score)),
                 by = .(frame = (region_position - 1L) %% 3L)]
  splicing <- truth[, .(reads = sum(score)),
                    by = .(spliced = grepl("N", cigar, fixed = TRUE))]

  expect_equal(frame[frame == 0L, reads], 5L)
  expect_equal(splicing[spliced == TRUE, reads], 2L)
  expect_equal(splicing[spliced == FALSE, reads], 4L)

  end_nt <- env$summarize_fragment_end_nt(truth)
  expect_equal(end_nt[end == "5-prime" & nucleotide == "A", reads], 3L)
  expect_equal(end_nt[end == "3-prime" & nucleotide == "C", reads], 2L)
  expect_equal(sum(end_nt[end == "5-prime", reads]), sum(truth$score))
})

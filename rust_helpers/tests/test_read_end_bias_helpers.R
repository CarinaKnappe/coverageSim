repo_dir <- normalizePath(
  Sys.getenv("COVSIM_REPO", unset = testthat::test_path("..", "..")),
  mustWork = TRUE
)

source(file.path(repo_dir, "rust_helpers", "R", "Read_end_bias_helpers.R"))

test_that("read-end nucleotide counts compare last base to whole-read composition", {
  seqs <- c("AAG", "AAA", "CCG", "TTT")

  observed <- read_end_nt_counts_from_sequences(seqs, run_id = "toy")

  expect_equal(observed$runID, rep("toy", 4))
  expect_equal(observed$nt, c("A", "C", "G", "T"))
  expect_equal(observed$last_nt_count, c(1L, 0L, 2L, 1L))
  expect_equal(sum(observed$all_read_nt_count), 12L)

  g_row <- observed[nt == "G"]
  expect_equal(g_row$last_nt_fraction, 0.5)
  expect_equal(g_row$all_read_nt_fraction, 2 / 12)
  expect_equal(g_row$last_vs_all_enrichment, 3)
})

test_that("read-end nucleotide helper handles empty input", {
  observed <- read_end_nt_counts_from_sequences(character(), run_id = "empty")

  expect_equal(nrow(observed), 4)
  expect_true(all(observed$last_nt_count == 0L))
  expect_true(all(observed$all_read_nt_count == 0L))
})


test_that("read sequences can be reconstructed from reference coordinates", {
  fasta_file <- tempfile(fileext = ".fa")
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(c(chr1 = "AACCGGTT")),
    fasta_file
  )
  Rsamtools::indexFa(fasta_file)

  chunk <- list(
    rname = factor(c("chr1", "chr1")),
    pos = c(2L, 2L),
    cigar = c("3M", "3M"),
    strand = factor(c("+", "-"))
  )

  observed <- read_sequences_from_reference(chunk, fasta_file)

  expect_equal(observed, c("ACC", "GGT"))
})

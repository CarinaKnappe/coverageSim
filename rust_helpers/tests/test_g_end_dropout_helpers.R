repo_dir <- normalizePath(
  Sys.getenv("COVSIM_REPO", unset = testthat::test_path("..", "..")),
  mustWork = TRUE
)

source(file.path(repo_dir, "rust_helpers", "R", "G_end_dropout_helpers.R"))

test_that("G-end dropout labels are stable", {
  expect_equal(g_end_dropout_label(0.30), "drop_G_30")
  expect_equal(g_end_dropout_label(0.60), "drop_G_60")
  expect_equal(g_end_dropout_label(1.00), "drop_G_100")
})

test_that("experiment paths are rewritten to the dropout dataset", {
  lines <- c(
    '"fasta","old/base/genome/genome.fa"',
    '"RFP","","1","WT","","old/base/reads/RFP_WT_1.bam"'
  )

  rewritten <- rewrite_experiment_base(lines, "old/base", "new/base")

  expect_equal(
    rewritten,
    c(
      '"fasta","new/base/genome/genome.fa"',
      '"RFP","","1","WT","","new/base/reads/RFP_WT_1.bam"'
    )
  )
})

test_that("terminal bases are reconstructed from reference and strand", {
  fasta_file <- tempfile(fileext = ".fa")
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(c(chr1 = "AACCGGTT")),
    fasta_file
  )
  Rsamtools::indexFa(fasta_file)

  chunk <- list(
    rname = factor(c("chr1", "chr1", "chr1")),
    pos = c(1L, 2L, 2L),
    cigar = c("3M", "3M", "3M"),
    strand = factor(c("+", "+", "-"))
  )

  expect_equal(
    reference_terminal_nt_from_alignment(chunk, fasta_file),
    c("C", "C", "T")
  )
})

test_that("reference-based G-end filter drops reconstructed G-ending reads", {
  fasta_file <- tempfile(fileext = ".fa")
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(c(chr1 = "AAAGCCCT")),
    fasta_file
  )
  Rsamtools::indexFa(fasta_file)

  chunk <- list(
    rname = factor(c("chr1", "chr1", "chr1")),
    pos = c(1L, 4L, 6L),
    cigar = c("3M", "3M", "3M"),
    strand = factor(c("+", "+", "-"))
  )

  filter <- make_reference_g_end_filter(
    fasta_file = fasta_file,
    drop_fraction = 1,
    seed = 1
  )

  keep <- filter$filter_fun(chunk)

  expect_equal(keep, c(TRUE, TRUE, FALSE))
  expect_equal(filter$stats$total, 3L)
  expect_equal(filter$stats$g_end, 1L)
  expect_equal(filter$stats$removed, 1L)
  expect_equal(filter$stats$kept, 2L)
})

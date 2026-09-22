test_that("secondary alignments always differ from their primary alignment", {
  # Two alignment positions of the same fragment length: each molecule has exactly
  # one alternative, the case where sample() used to draw from 1:value.
  fragments <- data.table::data.table(
    seqnames = "chr1", start = c(100L, 200L), end = c(127L, 227L), strand = "+",
    cigar = "28M", score = c(2L, 1L), fragment_id = c("f1", "f2"),
    fragment_length = 28L, sequence = "A", reference_sequence = "A"
  )
  set.seed(1)
  for (i in 1:60) {
    records <- simulate_alignment_artifacts(
      fragments, list(multimapping_rate = 1, secondary_alignments = 1)
    )
    secondary <- records[is_secondary == TRUE]
    primary <- records[is_secondary == FALSE & is_duplicate == FALSE]
    expect_gt(nrow(secondary), 0)
    partner <- primary[match(secondary$qname, primary$qname)]
    expect_true(all(secondary$start != partner$start))
  }
})

test_that("resample_values never treats one number as a range", {
  set.seed(2)
  expect_equal(resample_values(37L, 1L, replace = TRUE), 37L)
  expect_equal(sort(resample_values(c(5L, 9L), 2L)), c(5L, 9L))
  expect_length(resample_values(1:10, 4L), 4L)
})

test_that("residual autocorrelation works for transcripts with few sites", {
  set.seed(3)
  make_sites <- function(n) {
    data.table::data.table(
      transcript_id = "t1", site_tx = seq(1, by = 3, length.out = n),
      observed = stats::rpois(n, 3), expected = rep(3, n)
    )
  }
  for (n in c(1, 2, 3, 5, 9, 10, 11, 20)) {
    result <- estimate_residual_autocorrelation(make_sites(n), 9L)
    expect_equal(nrow(result$diagnostics), 9L)
    expect_equal(sum(result$kernel), 1)
    expect_false(anyNA(result$kernel))
  }
})

test_that("a BAM chunk without usable reads yields an empty read table", {
  raw <- list(
    rname = factor(c("chr1", "chr1")), pos = c(10L, 20L), cigar = c("28M", "28M"),
    strand = factor(c("+", "-"), levels = c("+", "-", "*")), flag = c(4L, 4L),
    mapq = c(30L, 30L), tag = list(NH = c(1L, 1L))
  )
  result <- filter_end_learning_reads(raw, 20)
  expect_equal(nrow(result$reads), 0L)
  expect_setequal(names(result$reads),
                  c("chromosome", "position", "five_position", "strand", "cigar",
                    "fragment_length", "count"))
  expect_equal(result$diagnostics[["bam_records"]], 2)
  expect_equal(result$diagnostics[["flag_or_mapq_or_NH_excluded"]], 2)
})

test_that("autocor_window keeps the requested amount of padding", {
  signal <- rep(0, 50)
  signal[25] <- 1
  expect_length(autocor_window(signal, 4, padding.rm = TRUE), 50)
  expect_length(autocor_window(signal, 4, padding.rm = 1), 52)
  expect_length(autocor_window(signal, 4, padding.rm = 3), 56)
  # At most max.lag positions of padding hold computed values; more is capped.
  for (keep in c(3, 4, 5, 8)) {
    expect_false(anyNA(autocor_window(signal, 4, padding.rm = keep)))
  }
  expect_length(autocor_window(signal, 4, padding.rm = 8), 50 + 2 * 4)
})

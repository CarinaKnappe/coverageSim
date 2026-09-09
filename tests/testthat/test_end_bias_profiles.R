test_that("synthetic end-bias profiles are valid and neutral defaults stay neutral", {
  profile <- make_synthetic_end_bias(
    k = 2L, enriched_kmer = "AC", enriched_weight = 4,
    depleted_kmer = "TT", depleted_weight = 0.25
  )
  expect_s3_class(profile, "data.table")
  expect_equal(nrow(profile), 16L)
  expect_equal(profile[kmer == "AC", weight], 4)
  expect_equal(profile[kmer == "TT", weight], 0.25)
  expect_no_error(normalize_fragment_geometry(list(
    five_prime_bias = list(source = "default"),
    three_prime_bias = list(source = "none")
  )))
  expect_error(make_synthetic_end_bias(k = 0), "positive integer")
  expect_error(make_synthetic_end_bias(k = 2, enriched_kmer = "N"), "invalid")
})

test_that("synthetic end bias changes fragment selection while preserving counts", {
  sequence <- paste(rep(c("A", "C", "G", "T"), 30L), collapse = "")
  exons <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(1L, nchar(sequence)), "+",
    exon_rank = 1L
  )
  model <- list(
    transcript_id = "tx1", exons = exons, cumulative_start = 1L,
    length = nchar(sequence), strand = "+", sequence = sequence
  )
  models <- list(tx1 = model)
  signal <- data.table::data.table(
    seqnames = "chr1", start = 50L, end = 50L, strand = "+",
    score = 10000L, transcript_id = "tx1", signal_position = 50L
  )
  candidates <- data.frame(
    fragment_length = c(8L, 8L), site_offset = c(2L, 3L),
    probability = c(0.5, 0.5)
  )
  site_tx <- genomic_site_to_transcript(model, signal$signal_position)
  starts <- site_tx - candidates$site_offset
  kmers <- vapply(starts, function(x) substr(model$sequence, x, x), character(1))
  enriched <- kmers[1]
  profile <- make_synthetic_end_bias(k = 1L, enriched_kmer = enriched,
                                     enriched_weight = 20)
  set.seed(17)
  fragments <- make_simulated_rpf_fragments(
    signal, models, 8L,
    list(source = "user", distribution = candidates,
         five_prime_bias = list(source = "user", table = profile),
         boundary_action = "error")
  )
  expect_equal(sum(fragments$score), 10000L)
  enriched_count <- sum(fragments$score[fragments$site_offset == candidates$site_offset[1]])
  expect_gt(enriched_count, 9000L)
  expect_true(all(fragments$end_bias_weight[fragments$site_offset == candidates$site_offset[1]] == 20))
})

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

end_selection_fixture <- function() {
  sequence <- paste(rep("A", 90), collapse = "")
  substr(sequence, 28, 28) <- "G"
  substr(sequence, 35, 35) <- "C"
  model <- list(transcript_id = "tx", exons = GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(1, 90), "+", exon_rank = 1L),
    cumulative_start = 1L, length = 90L, strand = "+", sequence = sequence)
  list(models = list(tx = model), signal = data.table::data.table(
    transcript_id = "tx", signal_position = c(30L, 60L), score = c(50000L, 50000L)))
}

run_end_selection <- function(fixture, five = 1, three = 0, distribution = NULL) {
  make_simulated_rpf_fragments(fixture$signal, fixture$models, 8L, list(
    source = "user", site_offset = if (is.null(distribution)) 2L else NULL,
    distribution = distribution,
    five_prime_bias = list(source = "user", strength = five,
      table = make_synthetic_end_bias(enriched_kmer = "G", enriched_weight = 4)),
    three_prime_bias = list(source = "user", strength = three,
      table = make_synthetic_end_bias(enriched_kmer = "C", enriched_weight = 4))))
}

test_that("end strength controls selection between codons at fixed geometry", {
  fixture <- end_selection_fixture()
  for (strength in c(0, 0.5, 1, 2)) {
    set.seed(103)
    result <- run_end_selection(fixture, five = strength)
    fraction <- sum(result[signal_position == 30, score]) / sum(result$score)
    expect_equal(fraction, 4^strength / (1 + 4^strength), tolerance = 0.01)
    expect_equal(sum(result$score), 100000)
    expect_true(all(result$signal_position %% 3 == 0))
    expect_true(all(result$site_offset == 2 & result$fragment_length == 8))
  }
  three <- run_end_selection(fixture, five = 0, three = 1)
  both <- run_end_selection(fixture, five = 1, three = 1)
  expect_equal(three[signal_position == 30, sum(score)] / 100000, 0.8, tolerance = 0.01)
  expect_equal(both[signal_position == 30, sum(score)] / 100000, 16/17, tolerance = 0.01)
  expect_equal(both[signal_position == 30, end_bias_weight], 16)
})

test_that("mixed geometry is weighted jointly and region budgets are retained", {
  fixture <- end_selection_fixture()
  fixture$signal <- data.table::rbindlist(list(fixture$signal, fixture$signal))
  fixture$signal[, sampling_group := rep(c("cds", "uorf"), each = 2)]
  fixture$signal[, read_count := rep(c(100000L, 20000L), each = 2)]
  fixture$signal[, score := rep(c(50000, 10000), each = 2)]
  distribution <- data.frame(fragment_length = c(8L, 9L), site_offset = c(2L, 3L),
                             probability = c(0.25, 0.75))
  set.seed(102)
  result <- run_end_selection(fixture, distribution = distribution)
  # Site 30 has total weight .25*4 + .75*1 = 1.75; site 60 has weight 1.
  expect_equal(result[signal_position == 30, sum(score)] / 120000,
               1.75/2.75, tolerance = 0.01)
  expect_equal(sum(result$score), 120000)
  expect_equal(result[sampling_group == "cds", sum(score)], 100000)
  expect_equal(result[sampling_group == "uorf", sum(score)], 20000)
  expect_equal(result[signal_position == 30 & fragment_length == 8, sum(score)] /
                 result[signal_position == 30, sum(score)], 1/1.75, tolerance = 0.01)
})

test_that("zero strength is identical to disabled end profiles", {
  fixture <- end_selection_fixture()
  set.seed(18)
  disabled <- make_simulated_rpf_fragments(fixture$signal, fixture$models, 8L,
    list(source = "user", site_offset = 2L))
  set.seed(18)
  zero <- run_end_selection(fixture, five = 0, three = 0)
  expect_equal(zero$score, disabled$score)
  expect_equal(zero$final_probability, disabled$final_probability)
  for (strength in list(-1, NA_real_, Inf, c(1, 2), "1")) {
    expect_error(run_end_selection(fixture, five = strength), "strength")
  }
})

test_that("latent coverage probabilities retain structural zeros and DMN variation", {
  set.seed(18)
  draws <- replicate(3000, draw_site_probabilities(c(2, 0, 6), TRUE))
  expect_true(all(draws[2, ] == 0))
  expect_equal(colSums(draws), rep(1, 3000), tolerance = 1e-12)
  expect_equal(mean(draws[1, ]), 0.25, tolerance = 0.015)
  expect_lt(abs(var(draws[1, ]) - 2*6/(8^2*9)), 0.003)
  expect_equal(sum(draw_site_probabilities(c(1e-24, 1e-24), TRUE)), 1)
})

test_that("joint end selection uses biological ends on the minus strand", {
  fixture <- end_selection_fixture()
  fixture$models$tx$strand <- "-"
  GenomicRanges::strand(fixture$models$tx$exons) <- "-"
  fixture$signal[, signal_position := 91L - signal_position]
  set.seed(33)
  result <- run_end_selection(fixture, five = 1, three = 1)
  expect_equal(result[signal_position == 61, sum(score)] / 100000,
               16/17, tolerance = 0.01)
  expect_true(all(result$five_prime_end > result$three_prime_end))
  expect_true(all(result$five_prime_end - result$signal_position == 2))
  expect_equal(result[signal_position == 61, sequence], "GAAAAAAC")
})

test_that("end profile normalization is repeatable and does not mutate input", {
  table <- data.table::data.table(kmer = "g", weight = 4)
  original <- data.table::copy(table)
  geometry <- normalize_fragment_geometry(list(five_prime_bias = list(
    source = "user", table = table, strength = 0.5)))
  expect_equal(table, original)
  expect_equal(normalize_fragment_geometry(geometry)$five_prime_bias,
               geometry$five_prime_bias)
})

test_that("length-specific codon and frame profiles select physical fragments", {
  fixture <- end_selection_fixture()
  substr(fixture$models$tx$sequence, 30L, 32L) <- "AAA"
  substr(fixture$models$tx$sequence, 60L, 62L) <- "CCC"
  fixture$signal[, `:=`(region_position = c(1L, 2L),
                        sampling_group = "cds", read_count = 100000L)]
  distribution <- data.frame(
    fragment_length = c(8L, 9L), site_offset = c(2L, 3L), probability = c(.5, .5)
  )
  geometry <- list(
    source = "user", distribution = distribution,
    codon_bias = list(source = "user", table = data.table::data.table(
      codon = rep(c("AAA", "CCC"), 2),
      fragment_length = rep(c(8L, 9L), each = 2),
      weight = c(4, 1, 1, 4)
    ))
  )
  set.seed(601)
  codon_result <- make_simulated_rpf_fragments(
    fixture$signal, fixture$models, 8:9, geometry
  )
  expect_equal(codon_result[fragment_length == 8 & codon == "AAA", sum(score)] /
                 codon_result[fragment_length == 8, sum(score)], 0.8,
               tolerance = 0.015)
  expect_equal(codon_result[fragment_length == 9 & codon == "AAA", sum(score)] /
                 codon_result[fragment_length == 9, sum(score)], 0.2,
               tolerance = 0.03)
  expect_equal(sum(codon_result$score), 100000)

  geometry$codon_bias <- list(source = "none")
  geometry$frame_bias <- list(source = "user", table = data.table::data.table(
    frame = rep(0:1, 2), fragment_length = rep(c(8L, 9L), each = 2),
    weight = c(4, 1, 1, 4)
  ))
  set.seed(602)
  frame_result <- make_simulated_rpf_fragments(
    fixture$signal, fixture$models, 8:9, geometry
  )
  expect_equal(frame_result[fragment_length == 8 & frame == 0, sum(score)] /
                 frame_result[fragment_length == 8, sum(score)], 0.8,
               tolerance = 0.015)
  expect_equal(frame_result[fragment_length == 9 & frame == 0, sum(score)] /
                 frame_result[fragment_length == 9, sum(score)], 0.2,
               tolerance = 0.03)
  expect_equal(sum(frame_result$score), 100000)
})

test_that("length-specific feature tables reject invalid rows", {
  expect_error(normalize_fragment_geometry(list(codon_bias = list(
    source = "user", table = data.frame(
      codon = "AXA", fragment_length = 28L, weight = 1
    )
  ))), "codons")
  expect_error(normalize_fragment_geometry(list(frame_bias = list(
    source = "user", table = data.frame(
      frame = 3L, fragment_length = 28L, weight = 1
    )
  ))), "frames")
})

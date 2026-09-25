test_that("BAM learning recovers both end biases with codon effects and both strands", {
  fixture <- make_end_learning_fixture()
  set.seed(212)
  bam <- make_learning_bam(fixture)
  fit <- fit_learning_fixture(fixture, bam)
  expect_s3_class(fit, "covsim_end_bias_fit")
  expect_equal(fit$five_prime_bias$source, "learned")
  expect_equal(relative_end_weight(fit$five_prime_bias, "G"), 4, tolerance = .1)
  expect_equal(relative_end_weight(fit$three_prime_bias, "C"), 3, tolerance = .1)
  expect_equal(fit$diagnostics$reads[["used"]], sum(fixture$signal$score))
  expect_equal(fit$diagnostics$reads[["unmatched_or_ambiguous"]], 0)
  expect_equal(fit$diagnostics$transcripts, 4)
  expect_equal(fit$diagnostics$convergence, 0L)
  expect_true(is.finite(fit$dmn_alpha_scale))
  expect_gt(fit$dmn_alpha_scale, 0)
  expect_equal(unique(fit$sequence_bias$dmn_alpha_scale), fit$dmn_alpha_scale)
  expect_equal(mean(fit$sequence_bias$alpha), 1)
  expect_equal(fit$diagnostics$dmn_alpha[usable == TRUE, .N], 4L)
  expect_s3_class(fit$auto_correlation, "covsim_autocorrelation")
  expect_equal(sum(fit$auto_correlation), 1)
  expect_equal(fit$diagnostics$auto_correlation$lag, 1:9)
  expect_setequal(fit$coverage_qc$summary$measure,
    c("reads", "sites", "zero_fraction", "variance_to_mean",
      "peak_fraction", "spike_fraction", "longest_zero_run"))
  # A biological codon effect must not simply be absorbed into the end profiles.
  codon_weights <- fit$diagnostics$codon_weights
  expect_gt(codon_weights[codon == "AAA", weight] / median(codon_weights$weight), 3.5)
  path <- tempfile(fileext = ".rds")
  saveRDS(fit, path)
  saved <- readRDS(path)
  other <- make_end_learning_fixture(seed = 209L)
  geometry <- other$geometry
  geometry$five_prime_bias <- saved$five_prime_bias
  geometry$three_prime_bias <- saved$three_prime_bias
  result <- make_simulated_rpf_fragments(other$signal, other$models, 28:29, geometry)
  expect_equal(sum(result$score), sum(other$signal$score))
  expect_gt(result[five_prime_kmer == "G", sum(score)] / sum(result$score), .45)
  expect_gt(result[three_prime_kmer == "C", sum(score)] / sum(result$score), .4)
  geometry$five_prime_bias$strength <- 0
  geometry$three_prime_bias$strength <- 0
  neutral <- make_simulated_rpf_fragments(other$signal, other$models, 28:29, geometry)
  expect_true(all(neutral$end_bias_weight == 1))
})

test_that("neutral data do not acquire strong end preferences", {
  fixture <- make_end_learning_fixture()
  set.seed(129)
  fit <- fit_learning_fixture(fixture, make_learning_bam(fixture, 1, 1))
  expect_true(all(abs(log(fit$five_prime_bias$table$weight)) < .15))
  expect_true(all(abs(log(fit$three_prime_bias$table$weight)) < .15))
})

test_that("length-specific preferences are learned separately", {
  fixture <- make_end_learning_fixture()
  geometry <- fixture$geometry
  geometry$five_prime_bias <- list(source = "user", table = data.table::data.table(
    kmer = "G", fragment_length = c(28L, 29L), weight = c(4, .5)))
  set.seed(151)
  fragments <- make_simulated_rpf_fragments(fixture$signal, fixture$models, 28:29, geometry)
  bam <- tempfile(fileext = ".bam")
  write_bam_library(simulated_rpf_alignments(fragments,
    GenomeInfoDb::seqinfo(fixture$transcripts)), bam)
  fit <- fit_learning_fixture(fixture, bam, by_length = TRUE)
  expect_equal(relative_end_weight(fit$five_prime_bias, "G"), c(4, .5), tolerance = .12)
  expect_setequal(fit$five_prime_bias$table$fragment_length, c(28L, 29L))
  expect_equal(fit$codon_bias$source, "learned")
  expect_setequal(fit$codon_bias$table$fragment_length, c(28L, 29L))
  expect_equal(fit$frame_bias$source, "learned")
  expect_setequal(fit$frame_bias$table$fragment_length, c(28L, 29L))
  expect_true(all(fit$diagnostics$frame_counts[frame == 0, weight] >
                  fit$diagnostics$frame_counts[frame != 0, weight]))
})

test_that("ambiguous annotations are excluded and geometry must be explicit", {
  fixture <- make_end_learning_fixture()
  bam <- make_learning_bam(fixture)
  expect_error(learn_end_bias(bam, fixture$fasta, fixture$transcripts, fixture$cds,
                             list(source = "default")), "explicit")
  geometry <- fixture$geometry
  geometry$distribution$fragment_length <- c(28L, 28L)
  expect_error(learn_end_bias(bam, fixture$fasta, fixture$transcripts, fixture$cds,
                             geometry), "one offset")
  expect_error(fit_learning_fixture(fixture, bam, dmn_min_reads = 0),
               "dmn_min_reads")
  expect_error(fit_learning_fixture(fixture, bam, dmn_min_sites = 1),
               "dmn_min_sites")
  expect_error(fit_learning_fixture(fixture, bam, acf_max_lag = 0),
               "acf_max_lag")
  opportunities <- end_learning_opportunities(fixture$models, fixture$cds,
    fixture$geometry$distribution, 1L)
  reads <- read_end_learning_bam(bam, 20)$reads
  duplicate <- data.table::copy(opportunities)
  duplicate[, transcript_id := paste0(transcript_id, "_isoform")]
  counted <- count_end_learning_reads(data.table::rbindlist(list(opportunities, duplicate)), reads)
  expect_equal(counted$diagnostics[["used"]], 0)
  expect_equal(counted$diagnostics[["unmatched_or_ambiguous"]], sum(reads$count))
  expect_equal(learning_alignment_cigar(c("10=2X16M", "3S25M", "12M60N16M")),
               c("28M", NA_character_, "12M60N16M"))
})

test_that("local residual structure yields an autocorrelation kernel and roughness QC", {
  sites <- data.table::rbindlist(lapply(paste0("tx", 1:20), function(id) {
    data.table::data.table(
      transcript_id = id,
      site_tx = seq.int(1L, 180L, by = 3L),
      observed = rep(c(rep(20, 5), rep(0, 5)), 6),
      expected = rep(10, 60)
    )
  }))
  learned <- learn_coverage_structure(sites, 4L)
  expect_gt(learned$acf[lag == 1, correlation], 0.5)
  expect_equal(names(learned$kernel), as.character(-4:4))
  expect_equal(learned$qc$summary[measure == "zero_fraction", median], 0.5)
  expect_equal(learned$qc$summary[measure == "longest_zero_run", median], 5)

  set.seed(92)
  shuffled <- sites[, observed := sample(observed), by = transcript_id]
  random_acf <- learn_coverage_structure(shuffled, 4L)$acf
  expect_gt(learned$acf[lag == 1, correlation],
            random_acf[lag == 1, correlation] + 0.3)
})

test_that("BAM filters report excluded reads and retain explicit duplicates", {
  sam <- tempfile(fileext = ".sam")
  flags <- c(0L, 1024L, 256L, 2048L, 512L, 1L, 0L, 0L, 0L)
  cigars <- c(rep("28M", 8L), "2S26M")
  mapq <- c(rep(30L, 6L), 0L, 30L, 30L)
  nh <- c(rep(1L, 7L), 2L, 1L)
  records <- vapply(seq_along(flags), function(i) paste(
    paste0("r", i), flags[i], "chr1", 100L, mapq[i], cigars[i], "*", 0, 0,
    paste(rep("A", 28L), collapse = ""), "*", paste0("NH:i:", nh[i]), sep = "\t"), character(1))
  writeLines(c("@SQ\tSN:chr1\tLN:1000", records), sam)
  bam <- Rsamtools::asBam(sam, destination = tempfile())
  result <- read_end_learning_bam(bam, 20)
  expect_equal(sum(result$reads$count), 2L)
  expect_equal(unname(result$diagnostics), c(9, 6, 1))
})

test_that("two-base motifs can be fit without inventing unsupported enrichment", {
  set.seed(451)
  kmers <- end_motif_levels(2)
  data <- data.table::as.data.table(expand.grid(five = kmers, three = kmers,
    stringsAsFactors = FALSE))
  data[, `:=`(transcript_id = "tx", fragment_length = 28L, codon = "AAA")]
  data[, count := as.integer(100 * ifelse(five == "GG", 3, 1) *
                              ifelse(three == "TC", 2, 1))]
  fit <- fit_end_preferences(data, 2, FALSE, 1, 500)
  expect_equal(relative_end_weight(fit$five_prime_bias, "GG", "AA"), 3, tolerance = .03)
  expect_equal(relative_end_weight(fit$three_prime_bias, "TC", "AA"), 2, tolerance = .03)
  expect_error(fit_end_preferences(data, 2, FALSE, 1, 1), "converge")
})

test_that("the dmn_alpha in-frame-fraction diagnostic ignores off-transcript reads and never rescales", {
  # An earlier version divided the raw dmn_alpha_scale estimate by this
  # fraction, reasoning that reads excluded for landing off the exact codon
  # boundary would otherwise bias the moment estimator low. That was
  # reverted after independent review showed MOM concentration estimates
  # are essentially unaffected by uniform random read thinning, so "divide
  # by the excluded fraction" has no general justification -- this is now
  # diagnostic only (dispersion$scale is never touched) and its denominator
  # is restricted to reads actually overlapping a modeled transcript.
  models <- list(tx1 = list(
    exons = GenomicRanges::GRanges("chr1", IRanges::IRanges(1000, 2000), strand = "+")
  ))
  dispersion <- list(scale = 0.065, diagnostics = NULL, sites = NULL)
  reads <- data.table::data.table(
    chromosome = c("chr1", "chr1", "chr1", "chr1", "chr1"),
    position = c(1100L, 1200L, 5000L, 1300L, 900L),
    strand = c("+", "+", "+", "+", "+"),
    # The last read's aligned blocks (900-913 and 2014-2027, skipping the
    # 1100 nt intron) never touch the exon (1000-2000) even though its
    # full reference span (900-2027) would appear to -- must not count.
    cigar = c("28M", "28M", "28M", "30M", "14M1100N14M"),
    fragment_length = c(28L, 28L, 28L, 30L, 28L),  # length 30 is unsupported
    count = c(10L, 5L, 100L, 999L, 100L)
  )
  distribution <- data.table::data.table(fragment_length = 28L)
  # Of the length-28 (supported) reads, only positions 1100 and 1200
  # overlap the modeled transcript (count 10 + 5 = 15); position 5000 and
  # the spliced read at 900 are excluded from the denominator despite
  # having a supported length. "used" = 10 (only the position-1100 read
  # was on the exact codon boundary and usable).
  count_diagnostics <- c(used = 10L, unmatched_or_ambiguous = 5L + 100L + 999L + 100L)

  result <- correct_dmn_alpha_scale_for_smearing(
    dispersion, reads, distribution, count_diagnostics, models
  )
  expect_equal(result$raw_scale, 0.065)
  expect_equal(result$scale, 0.065)  # never rescaled
  expect_equal(result$in_frame_fraction, 10 / 15)
})

test_that("learn_end_bias() does not correct dmn_alpha_scale when nothing was excluded", {
  # The main fixture (see the first test in this file) has zero
  # unmatched_or_ambiguous reads -- every usable read is exactly in-frame --
  # so the correction must be a no-op there.
  fixture <- make_end_learning_fixture()
  set.seed(212)
  bam <- make_learning_bam(fixture)
  fit <- fit_learning_fixture(fixture, bam)
  expect_equal(fit$diagnostics$dmn_alpha_in_frame_fraction, 1)
  expect_equal(fit$diagnostics$dmn_alpha_raw_scale, fit$dmn_alpha_scale)
})

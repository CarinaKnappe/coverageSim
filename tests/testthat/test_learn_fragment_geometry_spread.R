make_two_offset_bam <- function(fixture, offsets, probabilities, length = 28L, seed = 900L) {
  set.seed(seed)
  geometry <- list(source = "user", distribution = data.frame(
    fragment_length = rep(length, length(offsets)),
    site_offset = offsets, probability = probabilities
  ))
  fragments <- make_simulated_rpf_fragments(fixture$signal, fixture$models, length, geometry)
  bam <- tempfile(fileext = ".bam")
  write_bam_library(
    simulated_rpf_alignments(fragments, GenomeInfoDb::seqinfo(fixture$transcripts)), bam
  )
  list(bam = bam, fragments = fragments)
}

test_that("a low offset_relative_threshold recovers a genuine, non-aliased offset spread", {
  fixture <- make_end_learning_fixture(seed = 900L)
  # 15 and 17 differ by 2 nt (not a multiple of 3), so they occupy different
  # reading frames and are cleanly distinguishable from each other.
  sim <- make_two_offset_bam(fixture, offsets = c(15L, 17L), probabilities = c(0.7, 0.3))
  truth <- sim$fragments[, .(reads = sum(score)), by = site_offset]
  truth[, probability := reads / sum(reads)]
  data.table::setorder(truth, site_offset)

  learned <- learn_fragment_geometry(
    sim$bam, fixture$transcripts, fixture$cds,
    min_length = 28L, max_length = 28L, min_reads_per_length = 1L,
    offset_search = 2L, offset_relative_threshold = 0.2
  )
  distribution <- learned$distribution[order(site_offset)]
  expect_equal(distribution$site_offset, truth$site_offset)
  expect_equal(distribution$probability, truth$probability, tolerance = 0.02)
  expect_equal(sum(distribution$probability), 1)
})

test_that("a candidate offset a whole number of codons away from the best one is not reported as independent jitter", {
  # 14 = 17 - 3 aliases the true secondary offset 17 (same reading frame);
  # only the better-supported of the two should survive.
  fixture <- make_end_learning_fixture(seed = 900L)
  sim <- make_two_offset_bam(fixture, offsets = c(15L, 17L), probabilities = c(0.7, 0.3))
  learned <- learn_fragment_geometry(
    sim$bam, fixture$transcripts, fixture$cds,
    min_length = 28L, max_length = 28L, min_reads_per_length = 1L,
    offset_search = 2L, offset_relative_threshold = 0.2
  )
  expect_false(14L %in% learned$distribution$site_offset)
  expect_setequal(learned$distribution$site_offset, c(15L, 17L))
})

test_that("a single true offset is not split into spurious secondary candidates at a low threshold", {
  fixture <- make_end_learning_fixture(seed = 901L)
  set.seed(902)
  bam <- make_learning_bam(fixture, five = 1, three = 1)
  learned <- learn_fragment_geometry(
    bam, fixture$transcripts, fixture$cds,
    min_length = 28L, max_length = 29L, min_reads_per_length = 1L,
    offset_search = 2L, offset_relative_threshold = 0.2
  )
  expect_equal(nrow(learned$distribution), 2L)
  # Row order follows first-seen BAM length order, not a sorted guarantee.
  expect_setequal(learned$distribution$fragment_length, c(28L, 29L))
  expect_setequal(learned$distribution$site_offset, c(15L, 16L))
  expect_equal(sum(learned$distribution$probability), 1)
})

test_that("offset_relative_threshold = 1 (default) keeps exactly the historical single offset per length, even under ties", {
  # This fixture happens to tie two offsets at frame_fraction == 1 for one
  # length (low read counts make perfect in-frame scores common); the default
  # threshold must still collapse to a single, deterministically chosen row.
  fixture <- make_end_learning_fixture(seed = 900L)
  bam <- make_learning_bam(fixture, five = 1, three = 1)
  learned <- learn_fragment_geometry(
    bam, fixture$transcripts, fixture$cds,
    min_length = 28L, max_length = 28L, min_reads_per_length = 1L, offset_search = 2L
  )
  expect_equal(nrow(learned$distribution), 1L)
  expect_equal(learned$distribution$probability, 1)
})

test_that("offset_search_window returns no candidates when the search window is too short for the fragment", {
  # length=10, expected offset=14 (a_site), offset_search=2: lower clips to
  # 12, upper clips to 9 -- lower > upper, so no offset fits inside this
  # fragment at all. seq.int(12, 9) would otherwise silently count downward
  # through 12,11,10 (all >= length_value, invalid) before reaching 9.
  expect_equal(offset_search_window(10L, 14L, 2L), integer(0))
  # Ordinary, well-supported lengths are unaffected.
  expect_equal(offset_search_window(28L, 15L, 2L), 13:17)
  expect_equal(offset_search_window(28L, 15L, 0L), 15L)
  # A length that only barely fits the window's lower edge still works.
  expect_equal(offset_search_window(14L, 14L, 2L), 12:13)
})

test_that("a sub-1 offset_relative_threshold does not reorder rows relative to internal processing order", {
  # At the default threshold (1), the pre-existing tie-break sort
  # (order(-frame_fraction, ...) before grouping) intentionally reorders by
  # score -- unrelated to this diff and already pinned by
  # test_learn_fragment_geometry.R. Below 1, no such global sort happens
  # before merge() attaches each length's read share; merge(sort = FALSE) is
  # what keeps that step from introducing an unrelated reorder-by-join-key.
  fixture <- make_end_learning_fixture(seed = 707L)
  set.seed(708)
  bam <- make_learning_bam(fixture, five = 1, three = 1)
  learned <- learn_fragment_geometry(
    bam, fixture$transcripts, fixture$cds,
    min_length = 28L, max_length = 29L, min_reads_per_length = 1L,
    offset_search = 2L, offset_relative_threshold = 0.5
  )
  internal_order <- unique(attr(learned, "diagnostics")$fragment_length)
  expect_equal(unique(learned$distribution$fragment_length), internal_order)
})

test_that("a multi-offset fragment_geometry runs end to end through the public simNGScoverage()", {
  # Uses a directly user-supplied 2-offset distribution (not one round-tripped
  # through learn_fragment_geometry() first) to isolate exactly the link
  # code review flagged as untested: does simNGScoverage() itself, the public
  # entrypoint, forward and honor a multi-offset fragment_geometry$distribution?
  set.seed(910)
  genome <- suppressWarnings(suppressMessages(simGenome(
    n = 4, out_dir = tempfile("genome-"), cds_length = rep(300, 4),
    max_uorfs = 0, debug_on = FALSE
  )))
  cds <- ORFik::loadRegion(genome["txdb"], "cds")
  counts <- SummarizedExperiment::SummarizedExperiment(
    assays = list(gene = matrix(20000L, length(cds), 1L,
                               dimnames = list(names(cds), "RFP_x"))),
    rowRanges = cds,
    colData = S4Vectors::DataFrame(libtype = factor("RFP"), condition = factor("x"),
                                   replicate = "1", row.names = "RFP_x")
  )
  counts <- suppressMessages(simCountTablesRegions(
    counts, regionsToSample = "cds",
    region_proportion = list(cds = list(RFP = 1)), sampling = c(RFP = "MN")
  ))
  geometry <- list(
    source = "user", site_reference = "a_site",
    distribution = data.frame(
      fragment_length = c(28L, 28L), site_offset = c(15L, 17L), probability = c(0.7, 0.3)
    ),
    boundary_action = "renormalize"
  )
  out_dir <- tempfile("reads-"); exp_dir <- tempfile("exp-")
  dir.create(out_dir); dir.create(exp_dir)
  experiment <- suppressMessages(simNGScoverage(
    genome, counts, out_dir = out_dir, exp_name = "spread",
    exp_save_dir = exp_dir, sampling = list(cds = list(RFP = "DMN")),
    seq_bias = NULL, read_lengths_per = list(RFP = 28L),
    fragment_geometry = geometry, ground_truth = TRUE,
    libFormats = list(RFP = "ofst"), validate = FALSE
  ))
  truth <- data.table::fread(list.files(out_dir, pattern = "ground_truth", full.names = TRUE)[1])
  observed <- truth[, .(reads = sum(score)), by = site_offset]
  observed[, probability := reads / sum(reads)]
  data.table::setorder(observed, site_offset)
  expect_equal(observed$site_offset, c(15L, 17L))
  expect_equal(observed$probability, c(0.7, 0.3), tolerance = 0.03)
})

test_that("offset_relative_threshold is validated", {
  fixture <- make_end_learning_fixture(seed = 900L)
  bam <- make_learning_bam(fixture, five = 1, three = 1)
  for (invalid in list(0, -0.1, 1.1, NA_real_, Inf, "0.5")) {
    expect_error(
      learn_fragment_geometry(bam, fixture$transcripts, fixture$cds,
                              min_reads_per_length = 1L, offset_relative_threshold = invalid),
      "offset_relative_threshold"
    )
  }
})

test_that("a learned spread simulates fragments matching its own recovered proportions", {
  # Confirms the learned multi-offset distribution round-trips correctly
  # through simNGScoverage()'s existing simulated-RPF machinery.
  fixture <- make_end_learning_fixture(seed = 900L)
  sim <- make_two_offset_bam(fixture, offsets = c(15L, 17L), probabilities = c(0.7, 0.3))
  learned <- learn_fragment_geometry(
    sim$bam, fixture$transcripts, fixture$cds,
    min_length = 28L, max_length = 28L, min_reads_per_length = 1L,
    offset_search = 2L, offset_relative_threshold = 0.2
  )
  learned$boundary_action <- "renormalize"
  learned$five_prime_bias <- list(source = "none")
  learned$three_prime_bias <- list(source = "none")
  set.seed(903)
  fragments <- make_simulated_rpf_fragments(
    fixture$signal, fixture$models, 28L, learned
  )
  observed <- fragments[, .(reads = sum(score)), by = site_offset]
  observed[, probability := reads / sum(reads)]
  data.table::setorder(observed, site_offset)
  expect_equal(observed$site_offset, c(15L, 17L))
  expect_equal(observed$probability, c(0.7, 0.3), tolerance = 0.03)
})

test_that("a length whose search window contains no valid offset fails with a clear error, not an internal one", {
  # offset_search_window(10L, 14L, 2L) is empty (see the dedicated unit test
  # above); this confirms the *pipeline* handles that cleanly end to end,
  # for both an all-degenerate and a mixed valid/degenerate request, instead
  # of the rbindlist() column-mismatch crash a previous review round found.
  fixture <- make_end_learning_fixture(seed = 707L)
  set.seed(708)
  geometry <- list(source = "user", distribution = data.frame(
    fragment_length = 10L, site_offset = 5L, probability = 1
  ))
  fragments <- make_simulated_rpf_fragments(fixture$signal, fixture$models, 10L, geometry)
  short_bam <- tempfile(fileext = ".bam")
  write_bam_library(
    simulated_rpf_alignments(fragments, GenomeInfoDb::seqinfo(fixture$transcripts)), short_bam
  )

  expect_error(
    learn_fragment_geometry(short_bam, fixture$transcripts, fixture$cds,
                            min_length = 10L, max_length = 10L,
                            min_reads_per_length = 1L, offset_search = 2L),
    "No unambiguous CDS frame evidence"
  )

  long_bam <- make_learning_bam(fixture, five = 1, three = 1)
  mixed_bam <- tempfile(fileext = ".bam")
  Rsamtools::mergeBam(c(short_bam, long_bam), mixed_bam,
                      indexDestination = TRUE, overwrite = TRUE)
  expect_error(
    learn_fragment_geometry(mixed_bam, fixture$transcripts, fixture$cds,
                            min_length = 10L, max_length = 29L,
                            min_reads_per_length = 1L, offset_search = 2L),
    "No unambiguous CDS frame evidence"
  )
})

test_that("row order and per-length offset sets are both correct when several lengths each retain multiple frames", {
  # Two lengths, each with its own genuine (non-aliased) two-offset spread,
  # merged into one BAM: confirms .SD[order(...)][1L] neither reorders the
  # length groups nor drops a row when several lengths simultaneously
  # contribute more than one surviving offset.
  fixture <- make_end_learning_fixture(seed = 900L)
  sim28 <- make_two_offset_bam(fixture, offsets = c(15L, 17L),
                               probabilities = c(0.7, 0.3), length = 28L, seed = 900L)
  # default_ribosome_site_offset(29, "a_site") = 15, so its default search
  # window (offset_search = 2) is 13:17; 16 and 14 both fall inside it.
  sim29 <- make_two_offset_bam(fixture, offsets = c(16L, 14L),
                               probabilities = c(0.6, 0.4), length = 29L, seed = 901L)
  merged_bam <- tempfile(fileext = ".bam")
  Rsamtools::mergeBam(c(sim28$bam, sim29$bam), merged_bam,
                      indexDestination = TRUE, overwrite = TRUE)

  learned <- learn_fragment_geometry(
    merged_bam, fixture$transcripts, fixture$cds,
    min_length = 28L, max_length = 29L, min_reads_per_length = 1L,
    offset_search = 2L, offset_relative_threshold = 0.2
  )
  diagnostics <- attr(learned, "diagnostics")
  expect_equal(unique(learned$distribution$fragment_length), unique(diagnostics$fragment_length))

  # The exact, complete row order (not just the set of offsets per length):
  # within-group order follows the tie-break criteria, not probability
  # magnitude, so offset 17 legitimately sorts before 15 within length 28.
  expect_equal(learned$distribution$fragment_length, c(28L, 28L, 29L, 29L))
  expect_equal(learned$distribution$site_offset, c(17L, 15L, 16L, 14L))
  expect_equal(learned$distribution$probability,
              c(0.15, 0.35, 0.30, 0.20), tolerance = 0.03)
  expect_equal(sum(learned$distribution$probability), 1)
})

small_pipeline_fixture <- function(regions = "cds") {
  set.seed(9)
  genome <- suppressWarnings(suppressMessages(simGenome(
    n = 4, out_dir = tempfile("genome-"), cds_length = rep(300, 4),
    max_uorfs = 0, debug_on = FALSE
  )))
  cds <- ORFik::loadRegion(genome["txdb"], "cds")
  counts <- SummarizedExperiment::SummarizedExperiment(
    assays = list(gene = matrix(2000L, length(cds), 1L,
                                dimnames = list(names(cds), "RFP_x"))),
    rowRanges = cds,
    colData = S4Vectors::DataFrame(libtype = factor("RFP"), condition = factor("x"),
                                   replicate = "1", row.names = "RFP_x")
  )
  proportions <- list(leader = list(RFP = 0.2), cds = list(RFP = 0.8))
  counts <- suppressMessages(simCountTablesRegions(
    counts, regionsToSample = regions,
    region_proportion = if (length(regions) > 1) proportions else list(cds = list(RFP = 1)),
    sampling = c(RFP = "MN")
  ))
  list(genome = genome, counts = counts)
}

test_that("a region missing from `sampling` defaults to multinomial sampling", {
  fixture <- small_pipeline_fixture(c("leader", "cds"))
  out_dir <- tempfile("reads-"); exp_dir <- tempfile("exp-")
  dir.create(out_dir); dir.create(exp_dir)
  experiment <- suppressMessages(simNGScoverage(
    fixture$genome, fixture$counts, out_dir = out_dir, exp_name = "missing_sampling",
    exp_save_dir = exp_dir, fragment_mode = "legacy_point",
    sampling = list(cds = list(RFP = "DMN")), rnase_bias = list(RFP = 1),
    libFormats = list(RFP = "ofst"), validate = FALSE
  ))
  reads <- ORFik::fimport(ORFik::filepath(experiment, "default")[[1]])
  expect_equal(sum(S4Vectors::mcols(reads)$score), 4 * 2000)
})

test_that("MN sampling of a region ignores a real RNase kernel instead of crashing", {
  # A region-wide multi-element RNase kernel used to unconditionally extend
  # dt_range regardless of sampling mode, but only DMN sampling ever filled
  # those extra rows -- MN sampling never applies rnase_bias at all, so its
  # score vector stayed at the original (unextended) length and the later
  # dt_region[, score := sample] assignment crashed with a length mismatch.
  fixture <- small_pipeline_fixture("cds")
  out_dir <- tempfile("reads-"); exp_dir <- tempfile("exp-")
  dir.create(out_dir); dir.create(exp_dir)
  experiment <- suppressMessages(simNGScoverage(
    fixture$genome, fixture$counts, out_dir = out_dir, exp_name = "mn_with_rnase",
    exp_save_dir = exp_dir, fragment_mode = "legacy_point",
    sampling = list(cds = list(RFP = "MN")),
    rnase_bias = list(RFP = c(0.5, 2, 1, 10, 2, 1, 0.5)),
    libFormats = list(RFP = "ofst"), validate = FALSE
  ))
  reads <- ORFik::fimport(ORFik::filepath(experiment, "default")[[1]])
  expect_equal(sum(S4Vectors::mcols(reads)$score), 4 * 2000)
})

test_that("RNase extension follows the transcript direction on both strands", {
  # Two genes with six positions each: plus strand ascending, minus strand descending.
  dt_range <- data.table::data.table(
    genes = rep(1:2, each = 6), position = rep(1:6, 2),
    start = c(101:106, 206:201), strand = rep(c("+", "-"), each = 6)
  )
  dt_range[, end := start]
  result <- append_rnase_to_dt(dt_range, c(6, 6), list(RFP = c(0.5, 2, 1, 10, 2, 1, 0.5)))
  # Three extra positions at each end, moving away from the gene in transcript direction.
  expect_equal(result[genes == 1, start], 98:109)
  expect_equal(result[genes == 2, start], 209:198)
  expect_equal(nrow(result), 2 * (6 + 6))
})

test_that("a sequence bias table with missing motifs is rejected", {
  fixture <- small_pipeline_fixture()
  cds <- ORFik::loadRegion(fixture$genome["txdb"], "cds")
  lengths <- ORFik::widthPerGroup(cds, FALSE)
  dt_range <- data.table::data.table(genes = rep(seq_along(lengths), lengths))
  table <- load_seq_bias(type = "codon", shift = "p-site")
  table$variable <- NULL
  complete <- add_sequence_bias(fixture$genome, data.table::copy(dt_range), table, cds, lengths, "cds")
  expect_equal(unname(lengths(complete)), unname(lengths / 3))
  expect_error(
    add_sequence_bias(fixture$genome, data.table::copy(dt_range), table[seqs != "TGG"],
                      cds, lengths, "cds"),
    "does not contain every"
  )
  duplicated_rows <- rbind(table, table[seqs == "TGG"])
  expect_error(
    add_sequence_bias(fixture$genome, data.table::copy(dt_range), duplicated_rows,
                      cds, lengths, "cds"),
    "exactly one finite, positive alpha"
  )
  missing_motif_name <- data.table::copy(table)
  missing_motif_name[seqs == "TGG", seqs := NA_character_]
  expect_error(
    add_sequence_bias(fixture$genome, data.table::copy(dt_range), missing_motif_name,
                      cds, lengths, "cds"),
    "exactly one finite, positive alpha"
  )
  missing_alpha <- data.table::copy(table)
  missing_alpha[seqs == "TGG", alpha := NA_real_]
  expect_error(
    add_sequence_bias(fixture$genome, data.table::copy(dt_range), missing_alpha,
                      cds, lengths, "cds"),
    "exactly one finite, positive alpha"
  )
})

test_that("uORF reads can be split by uORF length", {
  uorfs <- GenomicRanges::GRangesList(
    tx1 = GenomicRanges::GRanges("chr1", IRanges::IRanges(100, width = 30), strand = "+"),
    tx1 = GenomicRanges::GRanges("chr1", IRanges::IRanges(200, width = 90), strand = "+")
  )
  assay <- matrix(40000L, 1, 1, dimnames = list("tx1", "RFP_x"))
  counts <- c(tx1 = 40000)
  set.seed(4)
  by_length <- distribute_reads_to_uORFs(counts, assay, uorfs, "character", "length")
  expect_equal(sum(by_length), 40000)
  expect_equal(unname(by_length) / 40000, c(0.25, 0.75), tolerance = 0.02)
  uniform <- distribute_reads_to_uORFs(counts, assay, uorfs, "character", "uniform")
  expect_equal(unname(uniform) / 40000, c(0.5, 0.5), tolerance = 0.02)
})

test_that("Dirichlet-multinomial region sampling handles regions with proportion zero", {
  set.seed(5)
  proportions <- c(leader = 1, cds = 0, trailer = 0, uorf = 0)
  expect_equal(sample_region_counts_dmn(50, proportions), c(50, 0, 0, 0))
  proportions <- c(leader = 0.3, cds = 0, trailer = 0.7, uorf = 0)
  for (i in 1:20) {
    draw <- sample_region_counts_dmn(200, proportions)
    expect_equal(sum(draw), 200)
    expect_true(all(draw[c(2, 4)] == 0))
    expect_false(anyNA(draw))
  }
  expect_error(sample_region_counts_dmn(10, c(1.1, -0.1, 0, 0)))
  counts <- suppressMessages(simCountTables(
    n = 5, libtypes = "CAGE", conditions = "WT", replicates = 2,
    betaLibSD = c(CAGE = 1), print_statistics = FALSE, plot_PCA = FALSE
  ))
  regions <- suppressMessages(simCountTablesRegions(
    counts, sampling = c(RFP = "MN", RNA = "MN", CAGE = "DMN", PAS = "MN")
  ))
  expect_false(anyNA(SummarizedExperiment::assay(regions, "leader")))
  expect_equal(
    Reduce(`+`, lapply(c("leader", "cds", "trailer", "uorf"),
                       function(x) SummarizedExperiment::assay(regions, x))),
    SummarizedExperiment::assay(regions, "gene")
  )
})

test_that("region_dmn_concentration dampens draw-to-draw region-proportion spread", {
  # Previously the region proportions (summing to 1) were used directly as
  # the Dirichlet alpha vector, i.e. always alpha0 = 1 -- extremely
  # overdispersed regardless of what the user wanted (SD ~0.31 on a mean-0.75
  # proportion; the default, concentration = 1, reproduces this unchanged).
  # A higher concentration scales alpha0 up, shrinking the spread towards
  # the plain-multinomial level.
  set.seed(1)
  proportions <- c(cds = 0.75, leader = 0.15, trailer = 0.1)
  draw_cds_share <- function(concentration, n = 300) {
    vapply(seq_len(n), function(i) {
      # sample_region_counts_dmn() returns an unnamed vector in the same
      # order as `proportions` (cds is position 1 here) -- not a names<-
      # oversight, matching this file's other existing test of it above.
      draw <- sample_region_counts_dmn(1000, proportions, concentration)
      draw[1] / sum(draw)
    }, numeric(1))
  }
  spread_default <- stats::sd(draw_cds_share(1))
  spread_high <- stats::sd(draw_cds_share(1000))
  expect_equal(mean(draw_cds_share(1)), 0.75, tolerance = 0.05)
  expect_gt(spread_default, 0.2)
  expect_lt(spread_high, 0.05)
  expect_gt(spread_default, spread_high * 5)

  expect_error(sample_region_counts_dmn(10, proportions, concentration = 0))
  expect_error(sample_region_counts_dmn(10, proportions, concentration = c(1, 2)))
  expect_error(sample_region_counts_dmn(10, proportions, concentration = NA_real_))

  # A concentration small enough to underflow rdirmnom()'s internal Gamma
  # draws to all-zero (0/0 = NaN) falls back to the mathematically correct
  # limit of a Dirichlet(alpha) draw as concentration -> 0: all reads land
  # on one randomly chosen region, chosen with probability equal to its own
  # share of `proportions` -- not silently invalid, and not a uniform
  # random pick among the regions either.
  set.seed(1)
  draws <- t(vapply(seq_len(400), function(i) {
    sample_region_counts_dmn(1000, proportions, concentration = 1e-6)
  }, numeric(3)))
  expect_true(all(rowSums(draws) == 1000))
  expect_true(all(apply(draws, 1, function(r) sum(r == 1000) == 1 && sum(r == 0) == 2)))
  winner_share <- colMeans(draws == 1000)
  # expect_equal()'s `tolerance` is a relative measure, not a per-category
  # absolute bound (e.g. it would accept c(.79, .13, .08) against these
  # proportions at tolerance = 0.07) -- assert the intended absolute
  # difference directly. The largest standard error at n = 400 is ~0.022,
  # so 0.07 is a generous, non-flaky bound.
  expect_true(all(abs(winner_share - unname(proportions)) < 0.07))

  # Threaded through simCountTablesRegions(): default unchanged, invalid rejected.
  counts <- suppressMessages(simCountTables(
    n = 5, libtypes = "CAGE", conditions = "WT", replicates = 2,
    betaLibSD = c(CAGE = 1), print_statistics = FALSE, plot_PCA = FALSE
  ))
  expect_error(suppressMessages(simCountTablesRegions(
    counts, sampling = c(RFP = "MN", RNA = "MN", CAGE = "DMN", PAS = "MN"),
    region_dmn_concentration = -1
  )))
  regions <- suppressMessages(simCountTablesRegions(
    counts, sampling = c(RFP = "MN", RNA = "MN", CAGE = "DMN", PAS = "MN"),
    region_dmn_concentration = 100
  ))
  expect_equal(
    Reduce(`+`, lapply(c("leader", "cds", "trailer", "uorf"),
                       function(x) SummarizedExperiment::assay(regions, x))),
    SummarizedExperiment::assay(regions, "gene")
  )
})

test_that("missing sequence lengths are filled from the FASTA index", {
  fasta <- tempfile(fileext = ".fa")
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(c(chrA = "ACGTACGTAC", chrB = "ACGT")), fasta)
  Rsamtools::indexFa(fasta)
  info <- GenomeInfoDb::Seqinfo(c("chrA", "chrB", "chrC"), c(NA, 4L, NA))
  filled <- fill_missing_seqlengths(info, fasta)
  expect_equal(unname(GenomeInfoDb::seqlengths(filled)), c(10L, 4L, NA))
  complete <- GenomeInfoDb::Seqinfo(c("chrA", "chrB"), c(10L, 4L))
  expect_identical(fill_missing_seqlengths(complete, fasta), complete)
})

test_that("BAM headers carry the real chromosome lengths", {
  fixture <- small_pipeline_fixture()
  for (mode in c("legacy_point", "simulated_rpf")) {
    out_dir <- tempfile("reads-"); exp_dir <- tempfile("exp-")
    dir.create(out_dir); dir.create(exp_dir)
    experiment <- suppressMessages(simNGScoverage(
      fixture$genome, fixture$counts, out_dir = out_dir, exp_name = paste0("header_", mode),
      exp_save_dir = exp_dir, fragment_mode = mode, sampling = list(cds = list(RFP = "DMN")),
      libFormats = list(RFP = "bam"), validate = FALSE
    ))
    header <- Rsamtools::scanBamHeader(ORFik::filepath(experiment, "default")[[1]])[[1]]$targets
    fasta_lengths <- GenomeInfoDb::seqlengths(Rsamtools::FaFile(fixture$genome["genome"]))
    expect_equal(header[names(fasta_lengths)], fasta_lengths)
  }
})

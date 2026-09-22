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

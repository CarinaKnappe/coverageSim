rnase_mode_fixture <- function() {
  set.seed(9)
  genome <- suppressWarnings(suppressMessages(simGenome(
    n = 4, out_dir = tempfile("genome-"), cds_length = rep(300, 4),
    max_uorfs = 0, debug_on = FALSE
  )))
  cds <- ORFik::loadRegion(genome["txdb"], "cds")
  gene <- matrix(2000L, length(cds), 2L,
                 dimnames = list(names(cds), c("RFP_x", "RNA_x")))
  counts <- SummarizedExperiment::SummarizedExperiment(
    assays = list(gene = gene, leader = gene * 0L, cds = gene),
    rowRanges = cds,
    colData = S4Vectors::DataFrame(
      libtype = factor(c("RFP", "RNA")), condition = factor(c("x", "x")),
      replicate = c("1", "1"), row.names = colnames(gene)
    )
  )
  list(genome = genome, counts = counts)
}

run_rnase_mode_fixture <- function(fixture, sampling, kernel = c(0.5, 2, 1, 10, 2, 1, 0.5),
                                   transcripts = rownames(fixture$counts)) {
  out_dir <- tempfile("reads-")
  exp_dir <- tempfile("exp-")
  dir.create(out_dir)
  dir.create(exp_dir)
  on.exit(unlink(c(out_dir, exp_dir), recursive = TRUE))
  experiment <- suppressMessages(simNGScoverage(
    fixture$genome, fixture$counts, out_dir = out_dir, exp_name = "rnase_modes",
    exp_save_dir = exp_dir, fragment_mode = "legacy_point",
    sampling = sampling, rnase_bias = list(RFP = kernel, RNA = 1),
    read_lengths_per = list(RFP = 27:29, RNA = 100),
    libFormats = list(RFP = "ofst", RNA = "ofst"),
    transcripts = transcripts, validate = FALSE
  ))
  vapply(ORFik::filepath(experiment, "default"), function(path) {
    sum(S4Vectors::mcols(ORFik::fimport(path))$score)
  }, numeric(1), USE.NAMES = FALSE)
}

test_that("an omitted active libtype defaults to MN and is zero-padded to match the DMN-extended region", {
  # RNA is left out of `sampling` for cds, so it defaults to MN while RFP
  # uses DMN on the same (RNase-extended) region -- this is also
  # simNGScoverage()'s own default sampling/rnase_bias combination for any
  # multi-libtype RFP+RNA experiment with real cds counts for both, so it
  # must keep working rather than require every libtype to opt into DMN.
  fixture <- rnase_mode_fixture()
  expect_equal(run_rnase_mode_fixture(fixture, list(cds = list(RFP = "DMN"))),
               c(8000, 8000))
})

test_that("sampling entries for absent libtypes do not affect RNase setup", {
  fixture <- rnase_mode_fixture()
  fixture$counts <- fixture$counts[, 1, drop = FALSE]
  expect_equal(run_rnase_mode_fixture(fixture, list(cds = list(RFP = "MN", RNA = "DMN"))),
               8000)
  expect_equal(run_rnase_mode_fixture(fixture, list(cds = list(RFP = "DMN", RNA = "MN"))),
               8000)
})

test_that("a libtype with zero CDS counts cannot trigger mixed-mode rejection", {
  fixture <- rnase_mode_fixture()
  SummarizedExperiment::assay(fixture$counts, "cds")[, "RNA_x"] <- 0L
  SummarizedExperiment::assay(fixture$counts, "leader")[, "RNA_x"] <- 2000L
  expect_equal(run_rnase_mode_fixture(fixture, list(cds = list(RFP = "DMN"))),
               c(8000, 8000))
})

test_that("a zero-count DMN libtype cannot extend an MN region", {
  fixture <- rnase_mode_fixture()
  SummarizedExperiment::assay(fixture$counts, "cds")[, "RFP_x"] <- 0L
  SummarizedExperiment::assay(fixture$counts, "leader")[, "RFP_x"] <- 2000L
  expect_equal(run_rnase_mode_fixture(fixture, list(cds = list(RFP = "DMN"))),
               c(8000, 8000))
})

test_that("RNase setup only considers counts in the selected transcripts", {
  fixture <- rnase_mode_fixture()
  SummarizedExperiment::assay(fixture$counts, "cds")[1:2, "RNA_x"] <- 0L
  SummarizedExperiment::assay(fixture$counts, "leader")[1:2, "RNA_x"] <- 2000L
  expect_equal(run_rnase_mode_fixture(fixture, list(cds = list(RFP = "DMN")),
                                     transcripts = rownames(fixture$counts)[1:2]),
               c(4000, 4000))
})

test_that("a scalar RNase kernel permits active mixed sampling modes", {
  fixture <- rnase_mode_fixture()
  expect_equal(run_rnase_mode_fixture(fixture, list(cds = list(RFP = "DMN")), kernel = 1),
               c(8000, 8000))
})

test_that("a zero-padded MN libtype places reads exactly at the true CDS start, not the RNase flank", {
  # Deterministic check that zero-padding lines up with the correct rows: a
  # spike ideal_coverage weight vector plus region_counts == 1 forces RNA's
  # single read per gene onto position 1 with no randomness left
  # (rmultinom(1, 1, weights) is fully determined once only one position has
  # nonzero weight). If the padding were misaligned (wrong reach, wrong
  # side, or applied to the wrong rows), this read would land inside the
  # RNase-only flank instead of the gene's true first CDS position.
  set.seed(9)
  genome <- suppressWarnings(suppressMessages(simGenome(
    n = 4, out_dir = tempfile("genome-"), cds_length = rep(300, 4),
    max_uorfs = 0, debug_on = FALSE
  )))
  cds <- ORFik::loadRegion(genome["txdb"], "cds")
  gene <- matrix(c(rep(2000L, 4), rep(1L, 4)), 4, 2,
                 dimnames = list(names(cds), c("RFP_x", "RNA_x")))
  counts <- SummarizedExperiment::SummarizedExperiment(
    assays = list(gene = gene, leader = gene * 0L, cds = gene),
    rowRanges = cds,
    colData = S4Vectors::DataFrame(
      libtype = factor(c("RFP", "RNA")), condition = factor(c("x", "x")),
      replicate = c("1", "1"), row.names = colnames(gene)
    )
  )
  out_dir <- tempfile("reads-"); exp_dir <- tempfile("exp-")
  dir.create(out_dir); dir.create(exp_dir)
  experiment <- suppressMessages(simNGScoverage(
    genome, counts, out_dir = out_dir, exp_name = "rnase_padding_position",
    exp_save_dir = exp_dir, fragment_mode = "legacy_point",
    sampling = list(cds = list(RFP = "DMN")),
    rnase_bias = list(RFP = c(0.5, 2, 1, 10, 2, 1, 0.5), RNA = 1),
    ideal_coverage = list(cds = list(
      RFP = quote(rep(c(1, 0, 0), length.out = x)),
      RNA = quote(c(1, rep.int(0, x - 1)))
    )),
    read_lengths_per = list(RFP = 27:29, RNA = 100),
    libFormats = list(RFP = "ofst", RNA = "ofst"), validate = FALSE
  ))
  paths <- ORFik::filepath(experiment, "default")
  rna_reads <- ORFik::fimport(paths[grepl("RNA_x", paths)])
  expect_equal(sum(S4Vectors::mcols(rna_reads)$score), 4)

  true_cds_start <- vapply(seq_along(cds), function(i) {
    exons <- cds[[i]]
    if (as.character(unique(GenomicRanges::strand(exons))) == "-") {
      max(GenomicRanges::end(exons))
    } else {
      min(GenomicRanges::start(exons))
    }
  }, numeric(1))
  observed_positions <- sort(GenomicRanges::start(rna_reads)[S4Vectors::mcols(rna_reads)$score > 0])
  expect_equal(observed_positions, sort(true_cds_start))
})

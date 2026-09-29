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

test_that("a region-nested rnase_bias is rejected instead of silently disabling the kernel", {
  # rnase_bias is keyed by library type, but ideal_coverage, auto_correlation,
  # sampling and simCountTablesRegions()'s region_proportion are all nested by
  # region, so the nested shape is an easy mistake. It used to pass silently:
  # rnase_bias[["RFP"]] was NULL, so no region was RNase-extended and
  # sim_sequence_bias() ran with rnase_acf = NULL. Measured on 4 genes with
  # cds-only DMN, that collapsed the reading-frame distribution from 49/13/13
  # to 68/0/0 and removed every read outside the CDS -- i.e. it quietly dropped
  # the manuscript's ~5:1:2 frame ratio, with no error and no warning.
  nested <- list(
    leader = list(RFP = 1),
    cds = list(RFP = c(0.5, 2, 1, 10, 2, 1, 0.5)),
    trailer = list(RFP = 1)
  )
  expect_error(validate_rnase_bias(nested), "not by region")
  expect_error(validate_rnase_bias(nested), "switches RNase smearing off everywhere")
  # The message has to name the offending entry and say what to write instead.
  expect_error(validate_rnase_bias(nested), "rnase_bias entries 'leader', 'cds', 'trailer' contain")
  expect_error(validate_rnase_bias(nested), "one kernel per library type")

  # Entries must be numeric: unlike ideal_coverage, a quoted expression is never
  # evaluated for this argument, so it would also switch smearing off.
  expect_error(
    validate_rnase_bias(list(RFP = quote(1 / seq.int(x)))),
    "is not a numeric vector"
  )
  expect_error(validate_rnase_bias(list(RFP = c(1, NA, 1))), "NA, NaN or Inf")
  expect_error(validate_rnase_bias(list(c(1, 2, 1))), "named list")

  # A duplicated library name is a third silent way to lose the kernel:
  # rnase_bias[["RFP"]] returns the first match, so this pair reads as "RFP has a
  # kernel" but resolves to NULL (measured: reach 0).
  expect_error(
    validate_rnase_bias(list(RFP = NULL, RFP = c(1, 2, 1))),
    "duplicates: RFP"
  )

  # An even-length kernel has no center position, so sim_sequence_bias() returns
  # one alpha value fewer than the row count that rnase_kernel_reach() creates
  # (measured: 10 alphas against 11 rows for a 9 nt region with c(1, 1)).
  expect_error(validate_rnase_bias(list(RFP = c(1, 1))), "which is an even number")
  expect_error(validate_rnase_bias(list(RFP = c(1, 1))), "shifted by one nucleotide")
  # The message must state the actual count, so the reader can see what to change.
  expect_error(validate_rnase_bias(list(RFP = c(1, 2, 1, 2))), "has 4 weights")

  # An all-zero kernel -- including the scalar 0, which reads like a natural way
  # to say "off" -- makes every smoothed weight 0, so the mean-preserving
  # rescale in sim_sequence_bias() divides by zero and the whole alpha vector
  # becomes NaN. That used to surface much later as a misleading
  # "dmn_alpha_scale produced invalid Dirichlet alpha values".
  expect_error(validate_rnase_bias(list(RFP = 0)), "no positive weight")
  # And it must point at the right way to switch smearing off.
  expect_error(validate_rnase_bias(list(RFP = 0)), "use 1\n?\\s*or NULL rather than 0")
  expect_error(validate_rnase_bias(list(RFP = c(0, 0, 0))), "no positive weight")
  expect_error(validate_rnase_bias(list(RFP = c(1, -1, 1))), "cannot be negative")

  # The supported flat shapes all stay valid, including the package default --
  # which is deliberately asymmetric, so symmetry must not be required.
  expect_silent(validate_rnase_bias(NULL))
  expect_silent(validate_rnase_bias(list(RFP = NULL)))
  expect_silent(validate_rnase_bias(list(RFP = 1, RNA = 1)))
  expect_silent(validate_rnase_bias(list(RFP = c(0, 1, 0))))
  expect_silent(validate_rnase_bias(
    list(RFP = c(0.5, 2, 1, 10, 2, 1, 0.5), RNA = 1, CAGE = 1, PAS = 1)
  ))
  expect_false(identical(c(0.5, 2, 1, 10, 2, 1, 0.5), rev(c(0.5, 2, 1, 10, 2, 1, 0.5))))

  # And simNGScoverage() rejects it up front: validate_rnase_bias() runs before
  # input_validation_controller(), so this fails on the rnase_bias shape rather
  # than on the deliberately bogus genome paths -- i.e. before any annotation is
  # loaded and before any read file is written.
  expect_error(
    simNGScoverage(
      simGenome = c(genome = "nope.fasta", gtf = "nope.gtf", txdb = "nope.db"),
      count_table = NULL, rnase_bias = nested
    ),
    "not by region"
  )
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

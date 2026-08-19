make_simulation_fixture <- function(max_uorfs = 0,
                                    regions = "cds",
                                    region_proportion = NULL,
                                    export_txdb = TRUE) {
  genome_dir <- tempfile("coverageSim-genome-")
  exp_dir <- tempfile("coverageSim-exp-")
  dir.create(genome_dir)
  dir.create(exp_dir)

  sim_genome <- suppressWarnings(simGenome(
    n = 6,
    out_dir = genome_dir,
    genome_name = paste0("artificial_uorf_", max_uorfs),
    max_uorfs = max_uorfs,
    cds_length = c(rep.int(330, 3), rep.int(300, 3)),
    export_txdb = export_txdb
  ))

  cds <- ORFik::loadRegion(sim_genome["txdb"], "cds")
  gene_count_table <- simCountTables(
    cds,
    libtypes = "RFP",
    print_statistics = FALSE,
    plot_PCA = FALSE,
    interceptMean = 10
  )
  region_args <- list(
    count_table = gene_count_table,
    regionsToSample = regions
  )
  if (!is.null(region_proportion)) {
    region_args$region_proportion <- region_proportion
  }
  region_count_table <- do.call(simCountTablesRegions, region_args)

  list(
    genome_dir = genome_dir,
    exp_dir = exp_dir,
    sim_genome = sim_genome,
    cds = cds,
    gene_count_table = gene_count_table,
    region_count_table = region_count_table
  )
}

run_simulated_experiment <- function(fixture,
                                     lib_formats = list(RFP = "ofst"),
                                     validate = FALSE, ...) {
  exp_name <- basename(tempfile("coverageSim-exp-"))
  simNGScoverage(
    fixture$sim_genome,
    fixture$region_count_table[, 1],
    exp_name = exp_name,
    exp_save_dir = fixture$exp_dir,
    libFormats = lib_formats,
    validate = validate,
    ...
  )
}

test_that("simNGScoverage retains the explicitly named legacy point mode", {
  fixture <- make_simulation_fixture(max_uorfs = 0, regions = "cds")
  experiment <- run_simulated_experiment(
    fixture,
    fragment_mode = "legacy_point"
  )
  imported <- ORFik::fimport(ORFik::filepath(experiment, "default")[1])

  expect_s4_class(imported, "GRanges")
  expect_false(is(imported, "GAlignments"))
  expect_true(all(GenomicRanges::width(imported) == 1L))
  expect_true(all(S4Vectors::mcols(imported)$size %in% 27:29))
})

test_that("simGenome exports the expected files for coding-only genomes", {
  fixture <- make_simulation_fixture(max_uorfs = 0)

  expect_setequal(names(fixture$sim_genome), c("genome", "gtf", "txdb"))
  expect_true(all(file.exists(unname(fixture$sim_genome))))
  expect_s4_class(ORFik::loadTxdb(fixture$sim_genome["txdb"]), "TxDb")
  expect_true(all(ORFik::widthPerGroup(fixture$cds, FALSE) %% 3 == 0))
})

test_that("simGenome supports multiple uORFs per transcript", {
  fixture <- make_simulation_fixture(max_uorfs = 2)

  expect_true("uorfs" %in% names(fixture$sim_genome))
  expect_true(file.exists(fixture$sim_genome[["uorfs"]]))

  uorf_ranges <- readRDS(fixture$sim_genome[["uorfs"]])
  expect_gt(length(uorf_ranges), 0)
  expect_true(all(ORFik::widthPerGroup(uorf_ranges, FALSE) %% 3 == 0))
})

test_that("simCountTables returns a summarized experiment with expected counts", {
  fixture <- make_simulation_fixture(max_uorfs = 0)
  counts <- fixture$gene_count_table

  expect_s4_class(counts, "RangedSummarizedExperiment")
  expect_equal(nrow(counts), length(fixture$cds))
  expect_equal(unique(as.character(SummarizedExperiment::colData(counts)$libtype)), "RFP")
  expect_true(all(SummarizedExperiment::assay(counts) >= 0))
})

test_that("simCountTablesRegions preserves total counts across sampled regions", {
  fixture <- make_simulation_fixture(
    max_uorfs = 1,
    regions = c("cds", "uorf"),
    region_proportion = list(
      cds = list(RFP = 0.8),
      uorf = list(RFP = 0.2)
    )
  )
  region_counts <- fixture$region_count_table

  expect_equal(SummarizedExperiment::assayNames(region_counts), c("gene", "cds", "uorf"))
  expect_equal(
    SummarizedExperiment::assay(region_counts, "gene"),
    SummarizedExperiment::assay(region_counts, "cds") +
      SummarizedExperiment::assay(region_counts, "uorf")
  )
})

test_that("simNGScoverage writes importable output for a small RFP simulation", {
  fixture <- make_simulation_fixture(max_uorfs = 1, regions = c("cds", "uorf"))
  experiment <- run_simulated_experiment(fixture)

  expect_s4_class(experiment, "experiment")

  default_path <- ORFik::filepath(experiment, "default")
  expect_true(all(file.exists(default_path)))

  imported <- ORFik::fimport(default_path)
  expect_s4_class(imported, "GAlignments")
  expect_gt(length(imported), 0)
})

test_that("simNGScoverage supports codon sequence bias across multiple CDS transcripts", {
  fixture <- make_simulation_fixture(max_uorfs = 0, regions = "cds")
  experiment <- simNGScoverage(
    fixture$sim_genome,
    fixture$region_count_table[, 1],
    exp_name = basename(tempfile("coverageSim-codon-bias-")),
    exp_save_dir = fixture$exp_dir,
    seq_bias = load_seq_bias(type = "codon", shift = "p-site", bias = "all"),
    libFormats = list(RFP = "ofst"),
    validate = FALSE
  )

  imported <- ORFik::fimport(ORFik::filepath(experiment, "default")[1])

  expect_s4_class(experiment, "experiment")
  expect_gt(length(imported), 0)
})

test_that("simNGScoverage supports default RFP frame coverage without sequence bias", {
  fixture <- make_simulation_fixture(max_uorfs = 0, regions = "cds")
  experiment <- simNGScoverage(
    fixture$sim_genome,
    fixture$region_count_table[, 1],
    exp_name = basename(tempfile("coverageSim-no-seq-bias-")),
    exp_save_dir = fixture$exp_dir,
    seq_bias = NULL,
    libFormats = list(RFP = "ofst"),
    validate = FALSE
  )

  imported <- ORFik::fimport(ORFik::filepath(experiment, "default")[1])

  expect_s4_class(experiment, "experiment")
  expect_gt(length(imported), 0)
})

test_that("simNGScoverage handles cds and uorf chromosome totals regardless of row order", {
  set.seed(303)
  fixture <- make_simulation_fixture(max_uorfs = 1, regions = c("cds", "uorf"))
  experiment <- run_simulated_experiment(fixture)

  expect_s4_class(experiment, "experiment")
  expect_true(all(file.exists(ORFik::filepath(experiment, "default"))))
})

test_that("simNGScoverage writes SAM output through the format-specific writer", {
  fixture <- make_simulation_fixture(max_uorfs = 1, regions = c("cds", "uorf"))
  experiment <- run_simulated_experiment(fixture, lib_formats = list(RFP = "sam"))

  sam_path <- ORFik::filepath(experiment, "default")
  expect_match(sam_path, "\\.sam$")
  expect_true(file.exists(sam_path))

  sam_lines <- readLines(sam_path)
  expect_true(any(grepl("^@SQ\\tSN:", sam_lines)))
  expect_gt(sum(!grepl("^@", sam_lines)), 0)

  converted_bam <- Rsamtools::asBam(sam_path)
  expect_true(file.exists(converted_bam))
  expect_gt(length(Rsamtools::scanBam(converted_bam)[[1]]$pos), 0)
})

test_that("simNGScoverage writes BAM output through the format-specific writer", {
  fixture <- make_simulation_fixture(max_uorfs = 1, regions = c("cds", "uorf"))
  experiment <- run_simulated_experiment(fixture, lib_formats = list(RFP = "bam"))

  bam_path <- ORFik::filepath(experiment, "default")
  expect_match(bam_path, "\\.bam$")
  expect_true(file.exists(bam_path))
  expect_gt(length(Rsamtools::scanBam(bam_path)[[1]]$pos), 0)
})

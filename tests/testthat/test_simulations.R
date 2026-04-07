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
  exp_name <- basename(tempfile("coverageSim-exp-"))

  experiment <- simNGScoverage(
    fixture$sim_genome,
    fixture$region_count_table[, 1],
    exp_name = exp_name,
    exp_save_dir = fixture$exp_dir,
    validate = FALSE
  )

  expect_s4_class(experiment, "experiment")

  default_path <- ORFik::filepath(experiment, "default")
  expect_true(all(file.exists(default_path)))

  imported <- ORFik::fimport(default_path)
  expect_s4_class(imported, "GRanges")
  expect_gt(length(imported), 0)
})

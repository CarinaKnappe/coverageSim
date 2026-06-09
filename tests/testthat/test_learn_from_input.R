make_learn_input_fixture <- function() {
  genome_dir <- tempfile("coverageSim-learn-genome-")
  exp_dir <- tempfile("coverageSim-learn-exp-")
  dir.create(genome_dir)
  dir.create(exp_dir)

  sim_genome <- suppressWarnings(simGenome(
    n = 6,
    out_dir = genome_dir,
    genome_name = "learn_input_fixture",
    max_uorfs = 0,
    cds_length = c(rep.int(330, 3), rep.int(300, 3)),
    export_txdb = TRUE
  ))

  cds <- ORFik::loadRegion(sim_genome["txdb"], "cds")
  gene_count_table <- simCountTables(
    cds,
    libtypes = "RFP",
    print_statistics = FALSE,
    plot_PCA = FALSE,
    interceptMean = 10
  )
  region_count_table <- simCountTablesRegions(
    count_table = gene_count_table,
    regionsToSample = "cds"
  )

  exp_name <- basename(tempfile("coverageSim-learn-exp-"))
  experiment <- simNGScoverage(
    sim_genome,
    region_count_table[, 1],
    exp_name = exp_name,
    exp_save_dir = exp_dir,
    libFormats = list(RFP = "ofst"),
    validate = FALSE
  )

  list(
    sim_genome = sim_genome,
    cds = cds,
    experiment = experiment
  )
}

test_that("learned real-input helpers create simulator-ready objects", {
  set.seed(404)
  fixture <- make_learn_input_fixture()
  reads <- ORFik::fimport(ORFik::filepath(fixture$experiment, "default")[1])

  read_lengths <- learn_read_lengths(
    reads,
    min_length = 1L,
    max_length = 100L,
    max_observations = 100L
  )

  expect_type(read_lengths, "integer")
  expect_gt(length(read_lengths), 0)

  count_table <- learn_cds_count_table(
    fixture$cds,
    reads,
    sample_name = "RFP_learned_1",
    condition = "learned",
    min_reads = 1L
  )

  expect_s4_class(count_table, "RangedSummarizedExperiment")
  expect_equal(SummarizedExperiment::assayNames(count_table), c("gene", "cds"))
  expect_equal(
    SummarizedExperiment::assay(count_table, "gene"),
    SummarizedExperiment::assay(count_table, "cds")
  )
  expect_gt(nrow(count_table), 0)

  seq_bias <- learn_codon_seq_bias(
    fixture$cds,
    reads,
    fa_file = fixture$sim_genome["genome"],
    min_tx_reads = 1L
  )

  expect_s3_class(seq_bias, "data.table")
  expect_setequal(names(seq_bias), c("variable", "seqs", "alpha"))
  expect_equal(nrow(seq_bias), length(Biostrings::GENETIC_CODE))
  expect_true(all(seq_bias$alpha > 0))
})

make_simulation_fixture <- function(max_uorfs = 0,
                                    regions = "cds",
                                    region_proportion = NULL,
                                    export_txdb = TRUE,
                                    seqnames = NULL) {
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
    seqnames = seqnames,
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

test_that("simGenome places multiple genes on compact chromosomes", {
  genome_dir <- tempfile("coverageSim-compact-genome-")
  dir.create(genome_dir)

  set.seed(101)
  simulated <- suppressWarnings(simGenome(
    n = 30,
    out_dir = genome_dir,
    genome_name = "compact_genome",
    max_uorfs = 0,
    export_txdb = FALSE,
    debug_on = FALSE
  ))

  genome <- Biostrings::readDNAStringSet(simulated[["genome"]])
  annotation <- rtracklayer::import(simulated[["gtf"]])
  genes <- annotation[annotation$type == "gene"]

  expect_length(genome, 24L)
  expect_equal(names(genome), paste0("chr", seq_len(24L)))
  expect_length(genes, 30L)
  expect_true(any(table(as.character(GenomicRanges::seqnames(genes))) > 1L))
  expect_true(all(vapply(
    split(genes, GenomicRanges::seqnames(genes)),
    function(x) GenomicRanges::isDisjoint(x, ignore.strand = TRUE),
    logical(1)
  )))
  genome_widths <- stats::setNames(Biostrings::width(genome), names(genome))
  maximum_gene_ends <- tapply(
    GenomicRanges::end(genes),
    as.character(GenomicRanges::seqnames(genes)),
    max
  )
  expect_true(all(
    maximum_gene_ends < genome_widths[names(maximum_gene_ends)]
  ))
})

test_that("simGenome supports explicit chromosomes and weighted gene assignment", {
  genome_dir <- tempfile("coverageSim-weighted-genome-")
  dir.create(genome_dir)

  set.seed(102)
  simulated <- suppressWarnings(simGenome(
    n = 12,
    out_dir = genome_dir,
    genome_name = "weighted_genome",
    seqnames = c("chrA", "chrB"),
    chromosome_weights = c(1, 0.01),
    max_uorfs = 0,
    export_txdb = FALSE,
    debug_on = FALSE
  ))

  genome <- Biostrings::readDNAStringSet(simulated[["genome"]])
  annotation <- rtracklayer::import(simulated[["gtf"]])
  genes <- annotation[annotation$type == "gene"]

  expect_equal(names(genome), c("chrA", "chrB"))
  expect_setequal(as.character(GenomicRanges::seqnames(genes)), names(genome))
  expect_gt(sum(as.character(GenomicRanges::seqnames(genes)) == "chrA"),
            sum(as.character(GenomicRanges::seqnames(genes)) == "chrB"))
})

test_that("simGenome retains the legacy one-gene-per-contig layout", {
  genome_dir <- tempfile("coverageSim-legacy-genome-")
  dir.create(genome_dir)

  simulated <- suppressWarnings(simGenome(
    n = 5,
    out_dir = genome_dir,
    genome_name = "legacy_genome",
    chromosome_layout = "legacy_one_gene_per_contig",
    max_uorfs = 0,
    export_txdb = FALSE,
    debug_on = FALSE
  ))

  genome <- Biostrings::readDNAStringSet(simulated[["genome"]])
  annotation <- rtracklayer::import(simulated[["gtf"]])
  genes <- annotation[annotation$type == "gene"]

  expect_equal(names(genome), paste0("chr", seq_len(5L)))
  gene_counts <- table(as.character(GenomicRanges::seqnames(genes)))
  expect_equal(names(gene_counts), paste0("chr", seq_len(5L)))
  expect_equal(as.integer(gene_counts), rep.int(1L, 5L))
})

test_that("multiple transcripts on one chromosome retain counts through BAM export", {
  set.seed(103)
  fixture <- make_simulation_fixture(
    max_uorfs = 1,
    regions = c("cds", "uorf"),
    region_proportion = list(
      cds = list(RFP = 0.9),
      uorf = list(RFP = 0.1)
    ),
    seqnames = "chr1"
  )
  experiment <- run_simulated_experiment(
    fixture,
    lib_formats = list(RFP = "bam"),
    ground_truth = TRUE
  )

  bam_path <- ORFik::filepath(experiment, "default")[1]
  truth_path <- sub("[.]bam$", "_ground_truth.tsv", bam_path)
  truth <- data.table::fread(truth_path)
  bam <- Rsamtools::scanBam(
    bam_path,
    param = Rsamtools::ScanBamParam(what = c("rname", "qname"))
  )[[1]]
  expected <- sum(SummarizedExperiment::assay(fixture$region_count_table[, 1], "cds")) +
    sum(SummarizedExperiment::assay(fixture$region_count_table[, 1], "uorf"))

  expect_equal(sum(truth$score), expected)
  expect_equal(length(bam$qname), expected)
  expect_equal(unique(as.character(bam$rname)), "chr1")
  expect_gt(data.table::uniqueN(truth$transcript_id), 1L)
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
    seq_bias = load_seq_bias(type = "codon", shift = "p-site", bias = "R2"),
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

test_that("active end selection runs through MN and DMN with conserved transcript counts", {
  set.seed(71)
  fixture <- make_simulation_fixture(max_uorfs = 0, regions = "cds")
  for (mode in c("MN", "DMN", "DMN_smearing")) {
    truth_dir <- tempfile("end-selection-truth-")
    experiment <- run_simulated_experiment(fixture,
      seq_bias = if (mode == "DMN_smearing") load_seq_bias(bias = "R10") else NULL,
      rnase_bias = list(RFP = if (mode == "DMN_smearing") c(0.5, 2, 1, 10, 2, 1, 0.5) else NULL),
      auto_correlation = NULL,
      sampling = list(cds = list(RFP = if (mode == "MN") "MN" else "DMN")),
      read_lengths_per = list(RFP = 28L), ground_truth = truth_dir,
      fragment_geometry = list(source = "user", site_offset = 15L,
        five_prime_bias = list(source = "user", strength = 1,
          table = make_synthetic_end_bias(enriched_weight = 4))))
    truth <- data.table::fread(list.files(truth_dir, full.names = TRUE)[1])
    expected <- SummarizedExperiment::assay(fixture$region_count_table, "cds")[, 1]
    observed <- truth[, .(count = sum(score)), by = transcript_id]
    expect_equal(observed$count[match(names(expected), observed$transcript_id)],
                 unname(expected))
    expect_true(all(truth$fragment_length == 28))
    imported <- ORFik::fimport(ORFik::filepath(experiment, "default")[1])
    expect_equal(sum(S4Vectors::mcols(imported)$score), sum(expected))
  }
})

test_that("dmn_alpha_scale controls positional roughness without changing read totals", {
  set.seed(720)
  fixture <- make_simulation_fixture(max_uorfs = 0, regions = "cds")
  for (assay_name in c("gene", "cds")) {
    SummarizedExperiment::assay(
      fixture$region_count_table, assay_name
    )[, 1] <- 10000L
  }
  run_scale <- function(scale, seq_bias = NULL) {
    set.seed(721)
    experiment <- run_simulated_experiment(
      fixture,
      fragment_mode = "legacy_point",
      seq_bias = seq_bias,
      ideal_coverage = list(cds = list(RFP = quote(rep(1, x)))),
      rnase_bias = list(RFP = NULL),
      auto_correlation = NULL,
      sampling = list(cds = list(RFP = "DMN")),
      dmn_alpha_scale = scale
    )
    reads <- ORFik::fimport(ORFik::filepath(experiment, "default")[1])
    data.table::data.table(
      seqnames = as.character(GenomicRanges::seqnames(reads)),
      score = S4Vectors::mcols(reads)$score
    )[, .(
      total = sum(score),
      peak_fraction = max(score) / sum(score),
      occupied_positions = .N
    ), by = seqnames]
  }

  rough <- run_scale(0.001)
  smooth <- run_scale(100)
  expect_equal(sum(rough$total), 60000)
  expect_equal(sum(smooth$total), 60000)
  expect_gt(stats::median(rough$peak_fraction),
            5 * stats::median(smooth$peak_fraction))
  expect_lt(stats::median(rough$occupied_positions),
            stats::median(smooth$occupied_positions))

  learned_profile <- load_seq_bias()
  learned_profile[, dmn_alpha_scale := 0.001]
  automatic <- run_scale(NULL, learned_profile)
  explicit <- run_scale(0.001, learned_profile)
  expect_equal(automatic, explicit)
})

test_that("shuffled input counts remain attached to their transcript identifiers", {
  set.seed(911)
  fixture <- make_simulation_fixture(max_uorfs = 0, regions = "cds", seqnames = "chr1")
  table <- fixture$region_count_table[c(6, 2, 4, 1, 5, 3), 1]
  values <- c(101L, 203L, 307L, 409L, 503L, 607L)
  for (assay in c("gene", "cds")) SummarizedExperiment::assay(table, assay)[, 1] <- values
  fixture$region_count_table <- table
  for (strength in c(0, 1)) {
    dir <- tempfile("shuffled-counts-truth-")
    run_simulated_experiment(fixture, seq_bias = NULL, auto_correlation = NULL,
      rnase_bias = list(RFP = NULL), ground_truth = dir,
      fragment_geometry = list(source = "user", site_offset = 15L,
        five_prime_bias = list(source = "user", strength = strength,
          table = make_synthetic_end_bias())))
    truth <- data.table::fread(list.files(dir, full.names = TRUE)[1])
    actual <- truth[, .(reads = sum(score)), by = transcript_id]
    expect_equal(actual$reads[match(rownames(table), actual$transcript_id)], values)
  }
})

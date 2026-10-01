test_that("CHOROS helper builds a zero-filled model universe", {
  helper_file <- testthat::test_path(
    "..", "..", "analysis", "benchmarks", "choros", "choros_utils.R"
  )
  helper_env <- new.env(parent = globalenv())
  sys.source(helper_file, envir = helper_env)

  tx_seqs <- c(tx1 = paste(rep(c("A", "C", "G", "T"), 100), collapse = ""))
  observed <- data.table::data.table(
    transcript_id = "tx1",
    cod_idx = 21L,
    d5 = 15L,
    d3 = 10L,
    utr5_len = 20L,
    cds_len = 300L,
    count = 7L
  )
  offsets <- data.table::data.table(
    qwidth = 28L,
    pos_frame = 0L,
    d5 = 15L
  )

  model <- helper_env$build_model_dt(
    observed, tx_seqs, offsets, excl_codons = 20L
  )

  expect_equal(nrow(model), 60L)
  expect_equal(model[cod_idx == 21L, count], 7L)
  expect_true(all(model[cod_idx != 21L, count] == 0L))
  expect_true(all(nchar(model$f5) == 3L))
  expect_true(all(nchar(model$f3) == 3L))
})

test_that("CHOROS resolves each coverageSim dataset and sample independently", {
  helper_env <- new.env(parent = globalenv())
  sys.source(testthat::test_path("..", "..", "analysis", "benchmarks", "choros", "choros_utils.R"),
             envir = helper_env)
  dataset <- tempfile("choros-covsim-")
  dir.create(dataset)
  config <- helper_env$resolve_covsim_config(
    dataset, file.path(dataset, "custom-output"), "RFP_KO_2"
  )

  expect_identical(config$base_dir, normalizePath(dataset))
  expect_identical(config$run_id, "RFP_KO_2")
  expect_identical(config$output_dir, file.path(dataset, "custom-output"))
  expect_identical(config$bam_file, file.path(dataset, "reads", "RFP_KO_2.bam"))
  expect_identical(
    config$fasta_file,
    file.path(dataset, "genome", "human_flavoured_sim.fasta")
  )
  expect_identical(
    config$txdb_file,
    file.path(dataset, "genome", "human_flavoured_sim.gtf.db")
  )
  expect_error(helper_env$resolve_covsim_config(""), "CHOROS_COVSIM_BASE")
})

test_that("CHOROS end profiles show complete symmetric flanks", {
  helper_file <- testthat::test_path(
    "..", "..", "analysis", "benchmarks", "choros", "choros_utils.R"
  )
  helper_env <- new.env(parent = globalenv())
  sys.source(helper_file, envir = helper_env)

  ends <- data.table::data.table(
    distance = c(-60L, -12L, -12L, 18L, 60L),
    count = c(1L, 2L, 3L, 4L, 5L)
  )
  profile <- helper_env$summarize_end_profile(ends, "distance")

  expect_equal(profile$position, -60:60)
  expect_equal(profile[position == -12L, count], 5L)
  expect_equal(profile[position == 0L, count], 0)
  expect_equal(profile[position == 60L, count], 5L)

  plot <- helper_env$plot_end_profile(ends, "distance", "Fragment end")
  expect_equal(plot$scales$get_scales("x")$limits, c(-60L, 60L))
  expect_equal(ggplot2::layer_data(plot, 1L)$x, -60:60)
  expect_true(any(vapply(plot$layers, function(layer) {
    identical(layer$geom$objname, "vline") ||
      inherits(layer$geom, "GeomVline")
  }, logical(1))))
})


test_that("real-human CHOROS separates preparation and modeling", {
  driver_file <- testthat::test_path("..", "..", "analysis", "benchmarks", "choros", "choros_real_human.R")
  expect_no_error(parse(driver_file))
  driver <- readLines(driver_file, warn = FALSE)
  expect_true(any(grepl("readRDS(prepared_file)", driver, fixed = TRUE)))
  expect_false(any(grepl("detect_ribosome_shifts_compat", driver, fixed = TRUE)))
  prepare_file <- testthat::test_path("..", "..", "analysis", "benchmarks", "choros", "prepare_choros_real_human.R")
  expect_no_error(parse(prepare_file))
  prepare <- readLines(prepare_file, warn = FALSE)
  expect_true(any(grepl("select_periodic_lengths", prepare, fixed = TRUE)))
  expect_true(any(grepl("infer_tis_offsets", prepare, fixed = TRUE)))
  expect_true(any(grepl(".sorted.unique_nh1.ofst", prepare, fixed = TRUE)))
})

test_that("CHOROS resolves unknown circularity without changing coordinates", {
  helper_file <- testthat::test_path(
    "..", "..", "analysis", "benchmarks", "choros", "choros_utils.R"
  )
  helper_env <- new.env(parent = globalenv())
  sys.source(helper_file, envir = helper_env)

  ranges <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chrM"),
    ranges = IRanges::IRanges(c(10L, 20L), width = 5L),
    seqinfo = GenomeInfoDb::Seqinfo(
      c("chr1", "chrM"),
      seqlengths = c(1000L, 100L),
      isCircular = c(NA, TRUE)
    )
  )
  normalized <- helper_env$set_unknown_circular_to_false(ranges)

  expect_identical(GenomicRanges::ranges(normalized), GenomicRanges::ranges(ranges))
  expect_equal(GenomeInfoDb::seqlengths(normalized), c(chr1 = 1000L, chrM = 100L))
  expect_identical(
    GenomeInfoDb::isCircular(normalized), c(chr1 = FALSE, chrM = FALSE)
  )
})


test_that("CHOROS maps read ends without duplicating transcript ranges", {
  helper_file <- testthat::test_path(
    "..", "..", "analysis", "benchmarks", "choros", "choros_utils.R"
  )
  helper_env <- new.env(parent = globalenv())
  sys.source(helper_file, envir = helper_env)

  reads <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges::IRanges(c(5L, 25L), width = 1L),
    strand = "+",
    seqinfo = GenomeInfoDb::Seqinfo(c("chr1", "unused_read_contig"))
  )
  transcripts <- GenomicRanges::GRangesList(
    tx1 = GenomicRanges::GRanges(
      seqnames = c("chr1", "chr1"),
      ranges = IRanges::IRanges(c(1L, 21L), c(10L, 30L)),
      strand = "+"
    )
  )

  result <- helper_env$map_read_ends_to_transcripts(reads, transcripts)

  expect_equal(start(result$mapped), c(5L, 15L))
  expect_equal(result$transcript_index, c(1L, 1L))
  expect_equal(mcols(result$mapped)$xHits, c(1L, 2L))
  expect_length(result$reads, 2L)
})

test_that("CHOROS aligns CDS and transcript regions by transcript name", {
  helper_file <- testthat::test_path(
    "..", "..", "analysis", "benchmarks", "choros", "choros_utils.R"
  )
  helper_env <- new.env(parent = globalenv())
  sys.source(helper_file, envir = helper_env)

  cds <- GenomicRanges::GRangesList(
    tx2 = GenomicRanges::GRanges("chr1", IRanges::IRanges(21L, 30L)),
    tx1 = GenomicRanges::GRanges("chr1", IRanges::IRanges(1L, 10L))
  )
  transcripts <- GenomicRanges::GRangesList(
    tx1 = GenomicRanges::GRanges("chr1", IRanges::IRanges(1L, 100L)),
    tx2 = GenomicRanges::GRanges("chr1", IRanges::IRanges(201L, 300L))
  )

  aligned <- helper_env$align_regions_by_name(
    cds, transcripts, c("tx1", "tx2")
  )

  expect_identical(names(aligned$regions), c("tx1", "tx2"))
  expect_identical(names(aligned$transcripts), c("tx1", "tx2"))
  expect_equal(start(aligned$regions[[1]]), 1L)
  expect_equal(start(aligned$transcripts[[1]]), 1L)
})

test_that("CHOROS derives plausible offsets from observed TIS peaks", {
  helper_env <- new.env(parent = globalenv())
  sys.source(testthat::test_path("..", "..", "analysis", "benchmarks", "choros", "choros_utils.R"),
             envir = helper_env)
  reads <- data.table::rbindlist(list(
    data.table::data.table(L = 28L, dist_start = -12L, score = rep(1, 120L)),
    data.table::data.table(L = 28L, dist_start = -9L, score = rep(1, 20L)),
    data.table::data.table(L = 29L, dist_start = -12L, score = rep(1, 100L))
  ))
  inferred <- helper_env$infer_tis_offsets(reads, c(28L, 29L))
  offsets <- helper_env$build_frame_offset_map(inferred$offsets)
  expect_equal(inferred$offsets$dist_start, c(-12L, -12L))
  expect_equal(nrow(offsets), 6L)
  expect_true(all(offsets$d5 >= 0L & offsets$d5 <= offsets$qwidth - 3L))
})

test_that("CHOROS periodic-length QC rejects weak frame enrichment", {
  helper_env <- new.env(parent = globalenv())
  sys.source(testthat::test_path("..", "..", "analysis", "benchmarks", "choros", "choros_utils.R"),
             envir = helper_env)
  reads <- data.table::data.table(
    L = c(rep(28L, 12L), rep(29L, 12L)),
    dist_start = c(rep(c(-12L, -11L, -10L), c(10L, 1L, 1L)),
                   rep(c(-12L, -11L, -10L), 4L)),
    score = 1
  )
  selected <- helper_env$select_periodic_lengths(
    reads, c(28L, 29L), min_reads = 10L, min_frame_prop = 0.6
  )
  expect_equal(selected$lengths, 28L)
})

test_that("coverageSim CHOROS scripts are configurable and parseable", {
  prepare_file <- testthat::test_path("..", "..", "analysis", "benchmarks", "choros",
                                     "prepare_choros_covsim_human_genome_only.R")
  model_file <- testthat::test_path("..", "..", "analysis", "benchmarks", "choros",
                                   "choros_covsim_human_genome_only.R")
  expect_no_error(parse(prepare_file))
  expect_no_error(parse(model_file))
  prepare <- readLines(prepare_file, warn = FALSE)
  model <- readLines(model_file, warn = FALSE)
  helper <- readLines(
    testthat::test_path("..", "..", "analysis", "benchmarks", "choros", "choros_utils.R"),
    warn = FALSE
  )
  expect_true(any(grepl("CHOROS_COVSIM_BASE", prepare, fixed = TRUE)))
  expect_true(any(grepl("CHOROS_COVSIM_BASE", model, fixed = TRUE)))
  expect_true(any(grepl("CHOROS_OUTPUT_DIR", prepare, fixed = TRUE)))
  expect_true(any(grepl("CHOROS_OUTPUT_DIR", model, fixed = TRUE)))
  expect_true(any(grepl("CHOROS_RUN_ID", prepare, fixed = TRUE)))
  expect_true(any(grepl("CHOROS_RUN_ID", model, fixed = TRUE)))
  expect_true(any(grepl("human_flavoured_sim.fasta", helper, fixed = TRUE)))
  expect_true(any(grepl("human_flavoured_sim.gtf.db", helper, fixed = TRUE)))
  expect_true(any(grepl("coverageSim_human_genome_only", model, fixed = TRUE)))
  expect_false(any(grepl("SRR32491292", c(prepare, model), fixed = TRUE)))
})

test_that("coverageSim-like reads produce valid CHOROS geometry", {
  helper_env <- new.env(parent = globalenv())
  sys.source(testthat::test_path("..", "..", "analysis", "benchmarks", "choros", "choros_utils.R"),
             envir = helper_env)
  reads <- data.table::data.table(
    L = rep(28L, 120L),
    dist_start = c(rep(-12L, 100L), rep(-9L, 20L)),
    score = 1
  )
  periodic <- helper_env$select_periodic_lengths(
    reads, 25:34, min_reads = 100L, min_frame_prop = 0.8
  )
  tis <- helper_env$infer_tis_offsets(
    reads, periodic$lengths, min_peak_reads = 50L
  )
  offsets <- helper_env$build_frame_offset_map(tis$offsets)
  geometry <- offsets[, .(d5, d3 = qwidth - d5 - 3L)]
  expect_equal(periodic$lengths, 28L)
  expect_equal(tis$offsets$dist_start, -12L)
  expect_true(all(geometry$d5 >= 0L & geometry$d3 >= 0L))
})

test_that("coverageSim preprocessing restores zero BAM header lengths", {
  prepare_file <- testthat::test_path(
    "..", "..", "analysis", "benchmarks", "choros", "prepare_choros_covsim_human_genome_only.R"
  )
  expressions <- parse(prepare_file)
  is_restore <- vapply(expressions, function(expr) {
    is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("restore_sequence_lengths"))
  }, logical(1))
  expect_equal(sum(is_restore), 1L)
  helper_env <- new.env(parent = globalenv())
  eval(expressions[[which(is_restore)]], envir = helper_env)

  reads <- suppressWarnings(GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(90L, width = 5L),
    seqinfo = GenomeInfoDb::Seqinfo("chr1", seqlengths = 0L)
  ))
  repaired <- helper_env$restore_sequence_lengths(reads, c(chr1 = 100L))
  expect_equal(GenomeInfoDb::seqlengths(repaired), c(chr1 = 100L))
  expect_identical(GenomeInfoDb::isCircular(repaired), c(chr1 = FALSE))
  expect_error(
    helper_env$restore_sequence_lengths(reads, c(chr1 = 92L)),
    "alignments remain outside"
  )
})

test_that("the CHOROS README names both workflows and the model selection rule", {
  readme_file <- testthat::test_path("..", "..", "analysis", "benchmarks", "choros", "README.md")
  expect_true(file.exists(readme_file))
  readme <- paste(readLines(readme_file, warn = FALSE), collapse = "\n")
  expect_match(readme, "SRR32491292.sorted.unique_nh1.ofst", fixed = TRUE)
  expect_match(readme, "RFP_WT_1.bam", fixed = TRUE)
  expect_match(readme, "prepare_choros_real_human.R", fixed = TRUE)
  expect_match(readme, "prepare_choros_covsim_human_genome_only.R", fixed = TRUE)
  expect_match(readme, "BIC_full < BIC_base", fixed = TRUE)
  # The two absolute /home/rstudio paths this used to require were removed from
  # the README on purpose: the scripts resolve their locations from environment
  # variables now, so a path from one machine no longer belongs in the docs.
  expect_false(grepl("/home/rstudio", readme, fixed = TRUE))
})

test_that("coverageSim keeps all abundant configured read lengths", {
  helper_env <- new.env(parent = globalenv())
  sys.source(testthat::test_path("..", "..", "analysis", "benchmarks", "choros", "choros_utils.R"),
             envir = helper_env)
  reads <- data.table::data.table(
    L = rep(28:30, each = 12L),
    dist_start = rep(c(-12L, -11L, -10L), 12L),
    score = 1
  )
  selected <- helper_env$select_periodic_lengths(
    reads, candidate_lengths = 28:30,
    min_reads = 10L, min_frame_prop = 0
  )
  expect_equal(selected$lengths, 28:30)
  prepare <- readLines(testthat::test_path(
    "..", "..", "analysis", "benchmarks", "choros", "prepare_choros_covsim_human_genome_only.R"
  ), warn = FALSE)
  expect_true(any(grepl("configured_lengths <- 28:30", prepare, fixed = TRUE)))
  expect_true(any(grepl("min_frame_prop = 0", prepare, fixed = TRUE)))
})

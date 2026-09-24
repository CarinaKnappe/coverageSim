test_that("assay_by_chromo() keeps the sample column name for a single-column DelayedMatrix", {
  # as.data.table() on a single-column DelayedMatrix collapses it via
  # as.array(x, drop = TRUE) to a bare vector, losing its column name -- only
  # reproducible with exactly one column, since a wider DelayedMatrix keeps
  # its real dimensions. This silently broke simNGScoverage()'s later
  # column lookup by sample name ("column not found: [RFP_x]").
  m <- matrix(2000L, 4L, 1L, dimnames = list(paste0("g", 1:4), "RFP_x"))
  seqnames_per <- rep("chr1", 4L)
  base_result <- assay_by_chromo(m, seqnames_per)
  delayed_result <- assay_by_chromo(DelayedArray::DelayedArray(m), seqnames_per)
  expect_true("RFP_x" %in% names(delayed_result))
  expect_equal(delayed_result, base_result)
})

test_that("assay_by_chromo() names its grouping column seqnamesPer for a single-gene assay", {
  # The nrow(assay_by_chromosome) == 1 branch (only one gene total in the
  # whole simulated genome) used to name this column seqnamesPerGroup, while
  # its only caller (simNGScoverage()) reads $seqnamesPer -- previously this
  # only resolved by data.table's $ partial-name matching rather than an
  # exact match, which is fragile (e.g. it would break if a second column
  # starting with "seqnamesPer" were ever added).
  m <- matrix(2000L, 1L, 1L, dimnames = list("g1", "RFP_x"))
  result <- assay_by_chromo(m, "chr1")
  expect_identical(names(result)[1], "seqnamesPer")
  expect_equal(result$seqnamesPer, "chr1")
})

test_that("assay_by_chromo() is unaffected for multi-column assays, delayed or not", {
  m <- matrix(2000L, 4L, 2L, dimnames = list(paste0("g", 1:4), c("RFP_x", "RNA_x")))
  seqnames_per <- rep("chr1", 4L)
  base_result <- assay_by_chromo(m, seqnames_per)
  delayed_result <- assay_by_chromo(DelayedArray::DelayedArray(m), seqnames_per)
  expect_equal(delayed_result, base_result)
})

test_that("a single-column DelayedMatrix experiment completes end to end", {
  set.seed(9)
  genome <- suppressWarnings(suppressMessages(simGenome(
    n = 4, out_dir = tempfile("genome-"), cds_length = rep(300, 4),
    max_uorfs = 0, debug_on = FALSE
  )))
  cds <- ORFik::loadRegion(genome["txdb"], "cds")
  gene <- matrix(2000L, length(cds), 1L, dimnames = list(names(cds), "RFP_x"))
  counts <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      gene = DelayedArray::DelayedArray(gene),
      leader = DelayedArray::DelayedArray(gene * 0L),
      cds = DelayedArray::DelayedArray(gene)
    ),
    rowRanges = cds,
    colData = S4Vectors::DataFrame(
      libtype = factor("RFP"), condition = factor("x"), replicate = "1",
      row.names = "RFP_x"
    )
  )
  out_dir <- tempfile("reads-"); exp_dir <- tempfile("exp-")
  dir.create(out_dir); dir.create(exp_dir)
  experiment <- suppressMessages(simNGScoverage(
    genome, counts, out_dir = out_dir, exp_name = "single_col_delayed",
    exp_save_dir = exp_dir, fragment_mode = "legacy_point",
    sampling = list(cds = list(RFP = "MN")), rnase_bias = list(RFP = 1),
    libFormats = list(RFP = "ofst"), validate = FALSE
  ))
  reads <- ORFik::fimport(ORFik::filepath(experiment, "default")[[1]])
  expect_equal(sum(S4Vectors::mcols(reads)$score), 4 * 2000)
})

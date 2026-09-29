source(testthat::test_path("..", "tests_and_demos", "real_input_validation.R"))

test_that("held-out scoring rewards known end preferences independently of depth", {
  data <- data.table::data.table(transcript_id = rep(c("train", "test"), each = 2),
    fragment_length = 28L, five = rep(c("G", "A"), 2), three = "A", codon = "AAA",
    count = c(400L, 100L, 40L, 10L))
  fit <- list(five_prime_bias = list(table = data.table::data.table(kmer = c("G", "A"), weight = c(4, 1))),
    three_prime_bias = list(table = data.table::data.table(kmer = "A", weight = 1)),
    diagnostics = list(codon_weights = data.table::data.table(codon = "AAA", weight = 1)))
  result <- end_validation_scores(data[transcript_id == "test"], fit)
  expect_gt(result[model == "all", log_likelihood], result[model == "neutral", log_likelihood])
  expect_equal(result[model == "all", log_likelihood], 40*log(.8) + 10*log(.2))
  expect_true(all(result$reads == 50))
  data[, partition := transcript_id]
  partitioned <- partitioned_end_validation_scores(data, fit)
  expect_true(all(partitioned[partition == "train", reads] == 500))
  expect_true(all(partitioned[partition == "test", reads] == 50))
})

test_that("optional human learning run holds out transcripts and preserves budgets", {
  path <- Sys.getenv("COVSIM_HUMAN_LEARNING_OUTPUT", "")
  skip_if(!nzchar(path), "Set COVSIM_HUMAN_LEARNING_OUTPUT for the real-data validation")
  split <- data.table::fread(file.path(path, "transcript_split.tsv"))
  expect_false(anyDuplicated(split$transcript_id) > 0)
  scores <- data.table::fread(file.path(path, "prediction_scores.tsv"))
  expect_setequal(scores$partition, c("train", "test"))
  expect_true(all(is.finite(scores$log_likelihood)))
  for (part in c("train", "test")) {
    scored_ids <- sub(":[0-9]+$", "", scores[partition == part, group])
    expect_true(all(scored_ids %in% split[partition == part, transcript_id]))
  }
  counts <- data.table::fread(file.path(path, "simulation_counts.tsv"))
  expect_true(all(counts$observed == counts$requested))
  expect_gte(data.table::uniqueN(counts$scenario), 6)
  for (name in unique(counts$scenario)) {
    bam <- list.files(file.path(path, "simulations", name), pattern = "\\.bam$", full.names = TRUE)
    expect_length(bam, 1)
    expect_equal(Rsamtools::countBam(bam)$records, counts[scenario == name, sum(requested)])
  }
  fit <- readRDS(file.path(path, "learned_end_bias.rds"))
  expect_equal(fit$diagnostics$convergence, 0L)
  expect_gt(fit$diagnostics$reads[["used"]], 0)
})

test_that("real-read projection handles splice junctions and the minus strand", {
  models <- list(
    plus = list(transcript_id = "plus", strand = "+", cumulative_start = c(1L, 11L),
      length = 20L, sequence = paste(rep("A", 20), collapse = ""),
      exons = GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 21), width = 10), "+")),
    minus = list(transcript_id = "minus", strand = "-", cumulative_start = c(1L, 11L),
      length = 20L, sequence = paste(rep("T", 20), collapse = ""),
      exons = GenomicRanges::GRanges("chr1", IRanges::IRanges(c(121, 101), width = 10), "-")))
  reads <- data.table::data.table(chromosome = "chr1", position = c(7L, 107L, 7L),
    five_position = c(7L, 124L, 7L), strand = c("+", "-", "+"),
    cigar = c("4M10N4M", "4M10N4M", "8M"), fragment_length = 8L, count = c(3, 5, 7))
  projected <- project_validation_reads(reads, models)
  expect_equal(projected$tx_start, c(7L, 7L))
  expect_equal(projected$count, c(3, 5))
  expect_setequal(projected$transcript_id, c("plus", "minus"))
})

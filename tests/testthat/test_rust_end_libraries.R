local_end_helpers <- function() {
  path <- test_path("..", "..", "analysis", "R", "end_bias_library_helpers.R")
  skip_if_not(file.exists(path), "Local RUST benchmark helpers are not distributed")
  env <- new.env(parent = globalenv())
  sys.source(path, envir = env)
  env
}

test_that("library selection applies independent biological-end weights across sites", {
  helper <- local_end_helpers()
  sequences <- c("GAC", "GAT", "AAC", "AAT")
  expect_equal(helper$end_library_weights(sequences, 4, 4), c(16,4,4,1))
  expect_equal(helper$end_library_weights(sequences, 4, 1), c(4,4,1,1))
  expect_equal(helper$end_library_weights(sequences, 1, 4), c(4,1,4,1))
  expect_error(helper$end_library_weights("NAC"), "DNA sequences")
  expect_error(helper$end_library_weights("GAC", 0), "positive")
  fragments <- data.table::data.table(
    transcript_id = rep(c("tx1", "tx2"), each = 4L),
    sequence = rep(sequences, 2L), tx_start = rep(1:4, 2L), score = 1L)
  before <- data.table::copy(fragments)
  sampled <- helper$resample_end_library(fragments, 4, 4, 100000L, seed = 42)
  expect_equal(fragments, before)
  expect_equal(sampled[, sum(score), by = transcript_id]$V1, c(100000,100000))
  expect_equal(sampled$score / 100000, rep(c(16,4,4,1)/25, 2), tolerance = 0.006)
  expect_equal(sampled, helper$resample_end_library(fragments, 4, 4, 100000L, seed = 42))
  neutral <- helper$resample_end_library(fragments, 1, 1, 100000L, seed = 42)
  expect_equal(neutral$score / 100000, rep(.25,8), tolerance = .006)
})

test_that("transcript projection preserves read ends across splicing on both strands", {
  helper <- local_end_helpers()
  sequence <- paste(rep("ACGT", 5), collapse = "")
  make_model <- function(starts, strand) list(
    exons = GenomicRanges::GRanges("chr1", IRanges::IRanges(starts, width = 10L), strand),
    strand = strand, sequence = sequence)
  models <- list(plus = make_model(c(101L,201L), "+"),
                 minus = make_model(c(401L,301L), "-"))
  truth <- data.table::data.table(transcript_id = c("plus", "minus"),
    five_prime_end = c(109L,402L), signal_position = c(201L,310L),
    fragment_length = 8L, site_offset = 2L,
    sequence = substr(sequence,9L,16L))
  projected <- helper$project_fragment_truth(truth, models)
  expect_equal(projected$tx_start, c(9L,9L))
  expect_equal(projected$tx_site, c(11L,11L))
  bad <- data.table::copy(truth)
  bad$five_prime_end[1] <- 150L
  expect_error(helper$project_fragment_truth(bad, models), "Non-exonic")
  bad <- data.table::copy(truth)
  bad$sequence[1] <- "AAAAAAAA"
  expect_error(helper$project_fragment_truth(bad, models), "sequence mismatch")
})

test_that("original RUST detects each injected read-end signal", {
  root <- Sys.getenv("COVSIM_RUST_END_OUTPUT")
  skip_if(!nzchar(root), "Set COVSIM_RUST_END_OUTPUT for the end-to-end benchmark")
  scenarios <- c("control", "five_prime", "three_prime", "both_ends")
  composition <- data.table::fread(file.path(root, "end_composition.tsv"))
  expect_true(all(composition$total_reads == 1500000))
  expect_true(all(composition$transcripts == 300L))
  for (scenario in scenarios) {
    qc <- jsonlite::fromJSON(file.path(root, scenario, "bam_validation.json"))
    expect_true(qc$exact_truth_match)
    expect_equal(qc$records, 1500000)
    expect_equal(qc$retained_transcripts, 300)
    expect_equal(qc$frame0_fraction, 1)
    for (mode in c("codon", "nucleotide")) {
      folder <- file.path(root, scenario, paste0("rust_", mode))
      raw <- list.files(folder, pattern = paste0("^RUST_", mode, "_file_"), full.names = TRUE)
      expect_length(raw, 1L)
      profile <- data.table::fread(raw, nrows = if (mode == "codon") 61L else 4L)
      expect_equal(nrow(profile), if (mode == "codon") 61L else 4L)
      expect_true(all(is.finite(as.matrix(profile[, -1, with = FALSE]))))
      plot <- list.files(folder, pattern = "metafootprint.*png$", full.names = TRUE)
      expect_length(plot, 1L)
      expect_gt(file.info(plot)$size, 1000)
    }
  }
  result <- data.table::fread(file.path(root, "rust_signal_summary.tsv"))
  for (scenario_name in scenarios[-1]) {
    current <- result[result$scenario == scenario_name, ]
    control <- result[result$scenario == "control", ]
    if (scenario_name %in% c("five_prime", "both_ends")) {
      expect_gt(current$five_G_ratio, 1.4)
      expect_gt(current$five_divergence, 5 * control$five_divergence)
    }
    if (scenario_name %in% c("three_prime", "both_ends")) {
      expect_gt(current$three_C_ratio, 1.4)
      expect_gt(current$three_divergence, 5 * control$three_divergence)
    }
    expect_lt(current$a_site_divergence, max(current$five_divergence, current$three_divergence) / 3)
  }
  five <- result[result$scenario == "five_prime", ]
  three <- result[result$scenario == "three_prime", ]
  expect_equal(five$max_nt_divergence_position, -15)
  expect_equal(three$max_nt_divergence_position, 12)
  expect_lt(five$three_divergence, five$five_divergence / 10)
  expect_lt(three$five_divergence, three$three_divergence / 10)
  plot_qc <- jsonlite::fromJSON(file.path(root, "plot_validation.json"))
  expect_true(plot_qc$all_plotted_feature_curves_match_upstream)
  expect_true(plot_qc$all_plotted_divergence_curves_match_upstream)
  expect_gt(file.info(file.path(root, "original_RUST_metafootprints.pdf"))$size, 1000)
})

# Optional end-to-end checks for the locally generated RUST benchmark.
test_that("large RUST benchmark retains transcripts and recovers the injected bias", {
  root <- Sys.getenv("COVSIM_RUST_LARGE_OUTPUT")
  skip_if(!nzchar(root), "Set COVSIM_RUST_LARGE_OUTPUT to validate the local benchmark")
  for (scenario in c("control_28", "strong_28", "strong_mixed")) {
    reads <- file.path(root, scenario, "reads")
    truth <- data.table::fread(file.path(reads, "RFP_benchmark_ground_truth.tsv"))
    expect_equal(sum(truth$score), 1500000)
    expect_equal(length(unique(truth$transcript_id)), 300L)
    expect_equal(nchar(truth$sequence), truth$fragment_length)
    expect_equal(Rsamtools::countBam(file.path(reads, "RFP_benchmark.bam"))$records, 1500000)
  }
  targets <- c("AAA", "GAA", "CCA", "CGT", "TTC", "GGT")
  analyses <- c("control_28/rust_correct", "strong_28/rust_correct",
                "strong_mixed/rust_correct", "strong_mixed/rust_common15")
  profiles <- list()
  for (analysis in analyses) {
    folder <- file.path(root, analysis)
    stats <- data.table::fread(file.path(folder, "RUST_transcript_stats.csv"))
    profile <- data.table::fread(file.path(folder, "RUST_A_site_profiles.csv"))
    expect_equal(nrow(stats), 300L)
    expect_true(all(stats$rust_included))
    expect_gt(sum(stats$reads), 1e6)
    expect_equal(nrow(profile), 61L)
    expect_true(all(is.finite(as.matrix(profile[, as.character(-40:19), with = FALSE]))))
    expect_true(all(file.info(file.path(folder, c("RUST_metafootprints.pdf",
                      "RUST_A_site_codon_ratios.pdf")))$size > 1000))
    if (grepl("common15", analysis)) {
      expect_true(all(stats$frame_0_ratio > 0.25 & stats$frame_0_ratio < 0.42))
    } else expect_true(all(stats$frame_0_ratio > 0.99))
    profiles[[analysis]] <- profile
  }
  expect_true(file.info(file.path(root, "RUST_comparison.pdf"))$size > 1000)
  expect_equal(nrow(data.table::fread(file.path(root, "comparison_summary.tsv"))), 4L)
  control <- profiles[[1]]
  expect_lt(max(abs(control[["0"]] / control$expected - 1)), 0.15)
  for (i in 2:3) {
    profile <- profiles[[i]]
    ratios <- profile[["0"]] / profile$expected
    near_background <- as.matrix(profile[, as.character(c(-40:-2, 2:19)), with = FALSE]) / profile$expected
    expect_lt(max(abs(near_background - 1)), 0.4)
    expect_gt(min(ratios[profile$codon %in% targets]), 3)
    expect_gt(mean(ratios[profile$codon %in% targets]),
              3 * mean(ratios[!profile$codon %in% targets]))
  }
})

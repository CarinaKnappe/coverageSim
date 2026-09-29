test_that("heterogeneous fragment generation preserves geometry and library weights", {
  python <- test_path("..", "..", "analysis", "benchmarks", "rust", ".venv-rust-original", "bin", "python")
  skip_if_not(file.exists(python), "Local RUST Python environment is unavailable")
  result <- system2(python, shQuote(test_path("natural_end_model_checks.py")), stdout = TRUE, stderr = TRUE)
  status <- attr(result, "status")
  if (is.null(status)) status <- 0L
  expect_equal(status, 0L, info = paste(result, collapse = "\n"))
})

test_that("heterogeneous libraries and upstream RUST show the intended end effects", {
  root <- Sys.getenv("COVSIM_NATURAL_END_OUTPUT")
  skip_if(!nzchar(root), "Set COVSIM_NATURAL_END_OUTPUT to validate this local benchmark")
  counts <- data.table::fread(file.path(root, "transcript_counts.tsv"))
  result <- data.table::fread(file.path(root, "signal_summary.tsv"))
  expect_equal(nrow(result), 16L)
  for (replicate_name in c("rep1", "rep2")) {
    expect_equal(sum(counts[[replicate_name]]), 1500000)
    expect_gt(stats::sd(counts[[replicate_name]]) / mean(counts[[replicate_name]]), 1)
    for (scenario_name in c("control", "five_prime", "three_prime", "both_ends")) {
      folder <- file.path(root, replicate_name, scenario_name)
      qc <- jsonlite::fromJSON(file.path(folder, "validation.json"))
      expect_equal(qc$records, 1500000)
      expect_true(qc$exact_truth_match)
      expect_setequal(names(qc$lengths), as.character(27:31))
      expect_gt(length(qc$offsets), 5)
      observed_counts <- unlist(qc$transcripts)[counts$transcript_id]
      expect_equal(as.numeric(observed_counts), as.numeric(counts[[replicate_name]]))
      frame0 <- qc$frame_common15[["0"]]/qc$records
      expect_gt(frame0, .4)
      expect_lt(frame0, .7)
      for (stratum_name in c("pooled", "length28")) {
        for (mode in c("codon", "nucleotide")) {
          path <- file.path(folder, paste0("rust_", stratum_name, "_", mode))
          raw <- list.files(path, pattern = paste0("^RUST_",mode,"_file_"), full.names = TRUE)
          expect_length(raw, 1L)
          profile <- data.table::fread(raw, nrows = if (mode == "codon") 61L else 4L)
          expect_equal(nrow(profile), if (mode == "codon") 61L else 4L)
          expect_true(all(is.finite(as.matrix(profile[,-1,with=FALSE]))))
        }
        row <- result[result$replicate == replicate_name & result$scenario == scenario_name & result$stratum == stratum_name, ]
        expect_gt(row$retained_reads, if (stratum_name == "pooled") 1e6 else 250000)
        if (scenario_name %in% c("five_prime", "both_ends")) expect_gt(row$five_excess, .1)
        if (scenario_name %in% c("three_prime", "both_ends")) expect_gt(row$three_excess, .1)
        if (scenario_name == "five_prime") expect_lt(abs(row$three_excess), .1)
        if (scenario_name == "three_prime") expect_lt(abs(row$five_excess), .1)
      }
    }
  }
  plots <- jsonlite::fromJSON(file.path(root, "plot_validation.json"))
  expect_equal(plots$profiles, 32)
  expect_true(plots$plotted_curves_match_upstream)
  expect_gt(file.info(file.path(root, "natural_RUST_metafootprints.pdf"))$size, 1000)
})

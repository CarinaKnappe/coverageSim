test_that("Figure 3C alpha calibration preserves reads and selects a valid profile", {
  out_dir <- Sys.getenv("COVSIM_FIGURE3C_OUTPUT", unset = "")
  testthat::skip_if(!nzchar(out_dir), "Set COVSIM_FIGURE3C_OUTPUT to validate this local analysis")

  required <- c(
    "figure3c_peak_positions.tsv", "figure3c_alpha_summary.tsv",
    "simulation_read_budgets.tsv", "figure3c_real_vs_selected.png",
    "figure3c_real_vs_selected.pdf", "figure3c_all_alpha_scales.png",
    "figure3c_all_alpha_scales.pdf", "manifest.tsv"
  )
  expect_true(all(file.exists(file.path(out_dir, required))))

  peaks <- data.table::fread(file.path(out_dir, "figure3c_peak_positions.tsv"))
  summary <- data.table::fread(file.path(out_dir, "figure3c_alpha_summary.tsv"))
  budgets <- data.table::fread(file.path(out_dir, "simulation_read_budgets.tsv"))
  expect_true(all(is.finite(peaks$relative_position)))
  expect_true(all(peaks$relative_position >= 0 & peaks$relative_position <= 100))
  expect_equal(length(unique(peaks$transcript_id[peaks$profile == "Real human reads"])),
               unique(summary[profile == "Real human reads", transcripts]))
  expect_equal(sum(summary$selected), 1L)
  expect_true(is.finite(summary[selected == TRUE, total_variation_to_real]))
  expect_lt(
    summary[selected == TRUE, mean_absolute_quantile_distance],
    summary[dmn_alpha_scale == 1, mean_absolute_quantile_distance]
  )
  expect_true(all(summary[!is.na(dmn_alpha_scale), simulation_replicates] >= 2L))
  expect_equal(budgets$observed_reads, budgets$expected_reads)
  expect_true(all(file.info(file.path(out_dir, required[4:7]))$size > 0))
})

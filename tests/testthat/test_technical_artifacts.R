make_artifact_records <- function() {
  fixture <- end_selection_fixture()
  fixture$signal[, score := c(10L, 10L)]
  fragments <- make_simulated_rpf_fragments(
    fixture$signal, fixture$models, 8L,
    list(source = "user", site_offset = 2L)
  )
  set.seed(901)
  simulate_alignment_artifacts(fragments, list(
    duplication_rate = 1, duplicate_copies = 2L,
    multimapping_rate = 1, secondary_alignments = 1L
  ))
}

test_that("PCR and multimapping artifacts create explicit SAM records", {
  records <- make_artifact_records()
  diagnostics <- attr(records, "diagnostics")
  expect_equal(diagnostics[["biological_molecules"]], 20)
  expect_equal(diagnostics[["pcr_duplicate_records"]], 40)
  expect_equal(diagnostics[["multimapping_molecules"]], 20)
  expect_equal(diagnostics[["secondary_alignment_records"]], 20)
  expect_equal(nrow(records), 80)
  expect_equal(sum(bitwAnd(records$flag, 1024L) != 0L), 40)
  expect_equal(sum(bitwAnd(records$flag, 256L) != 0L), 20)
  expect_equal(sum(records$is_duplicate), 40)
  expect_equal(sum(records$is_secondary), 20)
  expect_true(all(records[is_secondary == TRUE, nh] == 2L))
  expect_true(all(records[is_secondary == TRUE, qname] %in%
                  records[is_secondary == FALSE & is_duplicate == FALSE, qname]))
})

test_that("artifact BAM preserves duplicate flags, secondary flags, and NH tags", {
  records <- make_artifact_records()
  model <- end_selection_fixture()$models$tx
  seqinfo <- GenomeInfoDb::seqinfo(model$exons)
  GenomeInfoDb::seqlengths(seqinfo) <- 90L
  base <- tempfile("artifact-bam-")
  paths <- write_artifact_library(records, base, "bam", seqinfo)
  raw <- Rsamtools::scanBam(
    paths[["default"]],
    param = Rsamtools::ScanBamParam(what = c("qname", "flag"), tag = "NH")
  )[[1]]
  expect_length(raw$flag, 80)
  expect_equal(sum(bitwAnd(raw$flag, 1024L) != 0L), 40)
  expect_equal(sum(bitwAnd(raw$flag, 256L) != 0L), 20)
  expect_equal(sum(raw$tag$NH > 1L), 40)
})

test_that("technical artifact defaults are clean and settings are validated", {
  defaults <- normalize_technical_artifacts(NULL)
  expect_false(has_active_technical_artifacts(defaults))
  expect_equal(defaults$duplication_rate, 0)
  expect_equal(defaults$multimapping_rate, 0)
  expect_error(normalize_technical_artifacts(list(duplication_rate = 1.1)),
               "duplication_rate")
  expect_error(normalize_technical_artifacts(list(secondary_alignments = 0)),
               "secondary_alignments")
  expect_error(normalize_technical_artifacts(list(unknown = 1)), "Unknown")
})

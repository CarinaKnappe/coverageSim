region_point_annotation <- function(models, position) {
  ranges <- GenomicRanges::GRangesList(lapply(models, function(model) {
    genomic <- transcript_position_to_genomic(model, position)
    GenomicRanges::GRanges(
      as.character(GenomicRanges::seqnames(model$exons)[1]),
      IRanges::IRanges(genomic, width = 1L), strand = model$strand
    )
  }))
  names(ranges) <- names(models)
  ranges
}

make_region_learning_fixture <- function() {
  fixture <- make_end_learning_fixture(seed = 307L)
  sites <- c(leader = 20L, cds = 300L, trailer = 580L)
  scores <- c(leader = 100L, cds = 600L, trailer = 300L)
  signal <- data.table::rbindlist(lapply(fixture$models, function(model) {
    data.table::data.table(
      transcript_id = model$transcript_id,
      signal_position = transcript_position_to_genomic(model, sites),
      score = unname(scores)
    )
  }))
  set.seed(311)
  fragments <- make_simulated_rpf_fragments(
    signal, fixture$models, 28:29, fixture$geometry
  )
  bam <- tempfile(fileext = ".bam")
  write_bam_library(simulated_rpf_alignments(
    fragments, GenomeInfoDb::seqinfo(fixture$transcripts)
  ), bam)
  regions <- list(
    leader = region_point_annotation(fixture$models, sites[["leader"]]),
    cds = region_point_annotation(fixture$models, sites[["cds"]]),
    trailer = region_point_annotation(fixture$models, sites[["trailer"]])
  )
  list(fixture = fixture, bam = bam, regions = regions)
}

test_that("region proportions use unique A-sites across strands and splice junctions", {
  input <- make_region_learning_fixture()
  fit <- learn_region_proportions(
    input$bam, input$fixture$transcripts, input$regions,
    input$fixture$geometry, pseudocount = 0
  )
  expect_s3_class(fit, "covsim_region_fit")
  expect_equal(fit$global$count, c(400, 2400, 1200))
  expect_equal(fit$global$proportion, c(0.1, 0.6, 0.3))
  expect_equal(sum(fit$diagnostics$reads[c(
    "uniquely_projected", "unmatched_or_ambiguous_transcript"
  )]), 4000)
  expect_equal(fit$diagnostics$reads[["assigned"]], 4000)
  expect_equal(vapply(fit$region_proportion, `[[`, numeric(1), "RFP"),
               c(leader = 0.1, cds = 0.6, trailer = 0.3))
  expect_true(all(fit$per_transcript[, abs(sum(proportion) - 1),
                                     by = transcript_id]$V1 < 1e-12))
})

test_that("overlapping regions are handled by priority, fractions, or exclusion", {
  input <- make_region_learning_fixture()
  input$regions$uorf <- input$regions$leader
  learn <- function(action) learn_region_proportions(
    input$bam, input$fixture$transcripts, input$regions,
    input$fixture$geometry, overlap_action = action, pseudocount = 0
  )
  priority <- learn("priority")
  expect_equal(priority$global[region == "uorf", count], 400)
  expect_equal(priority$global[region == "leader", count], 0)
  expect_equal(priority$diagnostics$reads[["assigned"]], 4000)

  fractional <- learn("fractional")
  expect_equal(fractional$global[region == "uorf", count], 200)
  expect_equal(fractional$global[region == "leader", count], 200)
  expect_equal(fractional$diagnostics$reads[["assigned"]], 4000)

  excluded <- learn("exclude")
  expect_equal(excluded$global[region %in% c("uorf", "leader"), sum(count)], 0)
  expect_equal(excluded$diagnostics$reads[["no_region_or_overlap_excluded"]], 400)
  expect_equal(excluded$diagnostics$reads[["assigned"]], 3600)
})

test_that("region learner validates geometry, annotations, and smoothing", {
  input <- make_region_learning_fixture()
  duplicate <- input$fixture$geometry
  duplicate$distribution$fragment_length <- c(28L, 28L)
  expect_error(learn_region_proportions(
    input$bam, input$fixture$transcripts, input$regions, duplicate
  ), "one offset")
  expect_error(learn_region_proportions(
    input$bam, input$fixture$transcripts, input$regions,
    input$fixture$geometry, pseudocount = -1
  ), "pseudocount")
  bad_regions <- input$regions
  names(bad_regions)[1] <- "intron"
  expect_error(learn_region_proportions(
    input$bam, input$fixture$transcripts, bad_regions, input$fixture$geometry
  ), "regions")
})

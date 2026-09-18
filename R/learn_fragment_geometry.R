#' Learn fragment lengths and reference-site offsets from a Ribo-seq BAM
#'
#' For each sufficiently supported fragment length, candidate offsets around a
#' biologically plausible default are scored by the fraction of uniquely
#' projected CDS sites in frame zero. The selected length frequencies and
#' offsets form a simulator-ready fragment geometry.
#'
#' @param bam Path to a genomic BAM containing complete single-end fragments.
#' @param transcripts Named GRangesList of transcript exons in transcript order.
#' @param cds Named GRangesList of complete CDS exons using the same identifiers.
#' @param site_reference Either a_site (default) or p_site.
#' @param min_length,max_length Inclusive fragment-length limits.
#' @param min_reads_per_length Minimum BAM records required for a length.
#' @param offset_search Number of nucleotides searched on each side of the
#'   default offset.
#' @param min_mapq Minimum mapping quality. Default 20.
#' @return A simulator-ready fragment_geometry list. The diagnostics attribute
#'   contains every tested length/offset combination.
#' @export
learn_fragment_geometry <- function(
    bam, transcripts, cds, site_reference = c("a_site", "p_site"),
    min_length = 25L, max_length = 34L, min_reads_per_length = 100L,
    offset_search = 2L, min_mapq = 20L) {
  site_reference <- match.arg(site_reference)
  validate_geometry_learning_input(
    bam, transcripts, cds, min_length, max_length,
    min_reads_per_length, offset_search, min_mapq
  )
  models <- transcript_geometry_models(transcripts)
  observed <- read_end_learning_bam(bam, min_mapq)
  reads <- observed$reads[
    fragment_length >= min_length & fragment_length <= max_length
  ]
  support <- reads[, .(read_count = sum(count)), by = fragment_length]
  support <- support[read_count >= min_reads_per_length]
  if (!nrow(support)) stop("No fragment length reached min_reads_per_length")
  bounds <- geometry_learning_cds_bounds(models, cds)
  diagnostics <- data.table::rbindlist(lapply(
    support$fragment_length, function(length_value) {
      expected <- default_site_offset(length_value, site_reference)
      offsets <- seq.int(max(0L, expected - offset_search),
                         min(length_value - 1L, expected + offset_search))
      data.table::rbindlist(lapply(offsets, function(offset) {
        distribution <- data.table::data.table(
          fragment_length = length_value, site_offset = offset, probability = 1
        )
        projected <- project_region_learning_reads(
          reads[fragment_length == length_value], models, distribution
        )$reads
        projected <- merge(projected, bounds, by = "transcript_id")
        projected <- projected[site_tx >= cds_start & site_tx <= cds_end]
        total <- sum(projected$count)
        in_frame <- projected[(site_tx - cds_start) %% 3L == 0L, sum(count)]
        data.table::data.table(
          fragment_length = length_value, site_offset = offset,
          expected_offset = expected,
          bam_reads = support[fragment_length == length_value, read_count],
          cds_reads = total, in_frame_reads = in_frame,
          frame_fraction = if (total > 0) in_frame / total else NA_real_
        )
      }))
    }
  ))
  selected <- diagnostics[is.finite(frame_fraction)][
    order(-frame_fraction, -cds_reads, abs(site_offset - expected_offset)),
    .SD[1L], by = fragment_length
  ]
  if (nrow(selected) != nrow(support)) {
    stop("No unambiguous CDS frame evidence was available for one or more lengths")
  }
  selected[, probability := bam_reads / sum(bam_reads)]
  geometry <- list(
    source = "learned", site_reference = site_reference,
    distribution = selected[, .(fragment_length, site_offset, probability)],
    boundary_action = "renormalize"
  )
  attr(geometry, "diagnostics") <- diagnostics
  attr(geometry, "read_diagnostics") <- observed$diagnostics
  class(geometry) <- c("covsim_fragment_geometry", "list")
  geometry
}

default_site_offset <- function(fragment_length, site_reference) {
  p_offset <- ifelse(fragment_length <= 27L, 11L,
    ifelse(fragment_length <= 30L, 12L, 13L))
  as.integer(p_offset + ifelse(site_reference == "a_site", 3L, 0L))
}

geometry_learning_cds_bounds <- function(models, cds) {
  data.table::rbindlist(lapply(models, function(model) {
    sites <- learning_cds_positions(model, cds[[model$transcript_id]])
    data.table::data.table(
      transcript_id = model$transcript_id,
      cds_start = sites[1L], cds_end = sites[length(sites)] + 2L
    )
  }))
}

validate_geometry_learning_input <- function(
    bam, transcripts, cds, min_length, max_length,
    min_reads_per_length, offset_search, min_mapq) {
  if (!is.character(bam) || length(bam) != 1L || !file.exists(bam)) {
    stop("bam must name an existing file")
  }
  for (ranges in list(transcripts, cds)) {
    if (!methods::is(ranges, "GRangesList") || !length(ranges) ||
        is.null(names(ranges)) || anyDuplicated(names(ranges))) {
      stop("transcripts and cds must be named, nonempty GRangesLists")
    }
  }
  if (!setequal(names(transcripts), names(cds))) {
    stop("Transcript and CDS names must match")
  }
  integer_scalar <- function(value, lower, upper = Inf) {
    is.numeric(value) && length(value) == 1L && is.finite(value) &&
      value >= lower && value <= upper && value == as.integer(value)
  }
  if (!integer_scalar(min_length, 1) || !integer_scalar(max_length, min_length)) {
    stop("min_length and max_length must be ordered positive integers")
  }
  if (!integer_scalar(min_reads_per_length, 1)) {
    stop("min_reads_per_length must be a positive integer")
  }
  if (!integer_scalar(offset_search, 0, 20)) {
    stop("offset_search must be an integer from zero to 20")
  }
  if (!integer_scalar(min_mapq, 0, 255)) {
    stop("min_mapq must be an integer from zero to 255")
  }
}

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
#' @param offset_relative_threshold Numeric in (0, 1], default 1. For each
#'   fragment length, every candidate offset whose in-frame fraction is at
#'   least this fraction of that length's best candidate is kept, instead of
#'   only the single best one. At the default (1), this keeps only ties for
#'   the best -- i.e. exactly one offset per length, as before. A value below
#'   1 (e.g. 0.5) recovers real footprint-to-footprint jitter in where the
#'   ribosome site sits relative to the fragment's 5' end, as a
#'   probability-weighted spread of offsets per length instead of collapsing
#'   it to one fixed value; \code{fragment_geometry$distribution} and
#'   simulated-RPF fragment placement already support this directly (see
#'   \code{simNGScoverage()}). Only \code{learn_end_bias()} and
#'   \code{learn_region_proportions()} currently require one offset per
#'   length; feed them a single-offset geometry if a learned spread from here
#'   is not directly reusable there. In-frame fraction can only identify the
#'   correct reading frame (offset modulo 3), not an exact absolute offset;
#'   two candidates a whole number of codons apart are collapsed to the
#'   better-supported one instead of being reported as independent jitter.
#' @return A simulator-ready fragment_geometry list. The diagnostics attribute
#'   contains every tested length/offset combination.
#' @export
learn_fragment_geometry <- function(
    bam, transcripts, cds, site_reference = c("a_site", "p_site"),
    min_length = 25L, max_length = 34L, min_reads_per_length = 100L,
    offset_search = 2L, min_mapq = 20L, offset_relative_threshold = 1) {
  site_reference <- match.arg(site_reference)
  validate_geometry_learning_input(
    bam, transcripts, cds, min_length, max_length,
    min_reads_per_length, offset_search, min_mapq, offset_relative_threshold
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
      expected <- default_ribosome_site_offset(length_value, site_reference)
      offsets <- offset_search_window(length_value, expected, offset_search)
      if (!length(offsets)) {
        # Keep the same columns as the populated branch below (zero rows),
        # so rbindlist() across lengths never sees a column mismatch even
        # when every (or just this) length has no valid candidate offset --
        # letting the existing is.finite(frame_fraction)/setequal() checks
        # further down report their normal, clear error instead of an
        # internal "column not found" failure.
        return(data.table::data.table(
          fragment_length = integer(0), site_offset = integer(0),
          expected_offset = integer(0), bam_reads = integer(0),
          cds_reads = integer(0), in_frame_reads = integer(0),
          frame_fraction = numeric(0)
        ))
      }
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
  candidates <- diagnostics[is.finite(frame_fraction)]
  if (offset_relative_threshold >= 1) {
    # Exactly the historical behavior: one deterministically tie-broken
    # offset per length, even when several candidates are genuinely tied
    # (which small BAMs can produce by chance).
    selected <- candidates[
      order(-frame_fraction, -cds_reads, abs(site_offset - expected_offset)),
      .SD[1L], by = fragment_length
    ]
  } else {
    candidates[, best_fraction := max(frame_fraction), by = fragment_length]
    selected <- candidates[frame_fraction >= offset_relative_threshold * best_fraction]
  }
  # Two candidate offsets differing by an exact multiple of 3 (whole codons)
  # show statistically indistinguishable in-frame fractions: "in-frame" is
  # inherently a period-3 property, blind to which specific codon is
  # correct. This method can reliably pin down the frame (offset mod 3), not
  # an exact absolute offset a further 3n nt away; a true secondary/jitter
  # position exactly 3n nt from the primary one cannot be recovered this
  # way. Collapse any such same-frame collision within a length down to its
  # single better-supported candidate, rather than reporting the weaker
  # alias as if it were independent evidence of genuine offset jitter.
  # .SD[order(...)][1L] (sorting *within* each group) rather than
  # order(...); .SD[1L] (sorting the whole table first) keeps groups in
  # their original appearance order instead of reordering them by score.
  selected[, offset_frame := site_offset %% 3L]
  selected <- selected[, .SD[
    order(-frame_fraction, -in_frame_reads, abs(site_offset - expected_offset))
  ][1L], by = .(fragment_length, offset_frame)]
  selected[, offset_frame := NULL]
  if (!setequal(selected$fragment_length, support$fragment_length)) {
    stop("No unambiguous CDS frame evidence was available for one or more lengths")
  }
  # probability = P(length), from each length's overall BAM read share, times
  # P(offset | length), from each kept candidate's share of in-frame reads
  # among the candidates kept for that same length (falling back to an equal
  # split if none of them have any in-frame reads, which is only possible for
  # kept candidates whose best_fraction is itself 0).
  selected[, offset_share := if (sum(in_frame_reads) > 0) {
    in_frame_reads / sum(in_frame_reads)
  } else rep(1 / .N, .N), by = fragment_length]
  length_totals <- unique(selected[, .(fragment_length, bam_reads)])
  length_totals[, length_share := bam_reads / sum(bam_reads)]
  # sort = FALSE: merge() otherwise reorders rows by the join key, changing
  # the historical first-seen-length row order for no statistical reason.
  selected <- merge(selected, length_totals[, .(fragment_length, length_share)],
                    by = "fragment_length", sort = FALSE)
  selected[, probability := length_share * offset_share]
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

# The set of candidate offsets to test for one fragment length: every
# integer within offset_search of the expected default, clipped to the valid
# range [0, length_value). Empty if the window falls entirely outside that
# range -- only possible for a length shorter than roughly (expected offset
# - offset_search); seq.int() would otherwise silently count downward
# through invalid (>= length_value) offsets instead of yielding none.
offset_search_window <- function(length_value, expected, offset_search) {
  lower <- max(0L, expected - offset_search)
  upper <- min(length_value - 1L, expected + offset_search)
  if (lower > upper) integer(0) else seq.int(lower, upper)
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
    min_reads_per_length, offset_search, min_mapq, offset_relative_threshold) {
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
  if (!is.numeric(offset_relative_threshold) || length(offset_relative_threshold) != 1L ||
      !is.finite(offset_relative_threshold) || offset_relative_threshold <= 0 ||
      offset_relative_threshold > 1) {
    stop("offset_relative_threshold must be one finite number in (0, 1]")
  }
}

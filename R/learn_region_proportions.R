#' Learn read allocation among transcript regions
#'
#' Assigns each supported single-end fragment to one transcript and then to a
#' leader, CDS, trailer or uORF using its length-specific reference-site offset.
#' Ambiguous transcript assignments are excluded. Region overlaps are handled
#' explicitly rather than counted more than once.
#'
#' @param bam Path to a genomic BAM containing complete single-end fragments.
#' @param transcripts Named GRangesList of transcript exons in transcript order.
#' @param regions Named list of GRangesList annotations. Supported names are
#'   leader, cds, trailer and uorf; entries may omit transcripts with no region.
#' @param fragment_geometry Geometry list with exactly one site_offset for each
#'   supported fragment_length.
#' @param min_mapq Minimum mapping quality. Default 20.
#' @param overlap_action How to handle a site in several supplied regions:
#'   priority (default), fractional, or exclude.
#' @param region_priority Region order used by overlap_action = "priority".
#' @param pseudocount Non-negative count added to each region before calculating
#'   global and per-transcript proportions. Default 0.5.
#' @return A covsim_region_fit list with simulator-ready region_proportion,
#'   global and per-transcript estimates, diagnostics and provenance.
#' @export
learn_region_proportions <- function(
    bam, transcripts, regions, fragment_geometry, min_mapq = 20,
    overlap_action = c("priority", "fractional", "exclude"),
    region_priority = c("uorf", "cds", "leader", "trailer"),
    pseudocount = 0.5) {
  overlap_action <- match.arg(overlap_action)
  validate_region_learning_input(
    bam, transcripts, regions, min_mapq, region_priority, pseudocount
  )
  geometry <- normalize_fragment_geometry(fragment_geometry)
  if (is.null(geometry$distribution)) {
    stop("Region learning requires an explicit length/offset distribution")
  }
  distribution <- validate_fragment_distribution(geometry$distribution)
  distribution <- distribution[probability > 0]
  if (anyDuplicated(distribution$fragment_length)) {
    stop("Region learning requires one offset per fragment length")
  }
  models <- transcript_geometry_models(transcripts)
  observed <- read_end_learning_bam(bam, min_mapq)
  projected <- project_region_learning_reads(
    observed$reads, models, distribution
  )
  intervals <- transcript_region_intervals(regions, models)
  assigned <- assign_projected_regions(
    projected$reads, intervals, names(regions), overlap_action, region_priority
  )
  estimates <- summarize_region_assignments(
    assigned$reads, names(models), names(regions), pseudocount
  )
  fit <- list(
    region_proportion = stats::setNames(lapply(
      estimates$global$proportion, function(value) list(RFP = value)
    ), estimates$global$region),
    global = estimates$global,
    per_transcript = estimates$per_transcript,
    diagnostics = list(
      reads = c(observed$diagnostics, projected$diagnostics,
                assigned$diagnostics),
      overlaps = assigned$overlaps
    ),
    provenance = list(
      bam = normalizePath(bam), transcript_ids = names(transcripts),
      distribution = distribution, overlap_action = overlap_action,
      region_priority = region_priority, min_mapq = min_mapq,
      pseudocount = pseudocount, created = as.character(Sys.time()),
      package_version = as.character(utils::packageVersion("coverageSim"))
    )
  )
  class(fit) <- c("covsim_region_fit", "list")
  fit
}

validate_region_learning_input <- function(
    bam, transcripts, regions, min_mapq, region_priority, pseudocount) {
  if (!is.character(bam) || length(bam) != 1L || !file.exists(bam)) {
    stop("bam must name an existing file")
  }
  ids <- names(transcripts)
  if (!methods::is(transcripts, "GRangesList") || !length(transcripts) ||
      is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("transcripts must be a nonempty GRangesList with unique names")
  }
  allowed <- c("leader", "cds", "trailer", "uorf")
  if (!is.list(regions) || !length(regions) || is.null(names(regions)) ||
      anyDuplicated(names(regions)) || !all(names(regions) %in% allowed)) {
    stop("regions must be a named list of leader, cds, trailer or uorf annotations")
  }
  for (name in names(regions)) {
    annotation <- regions[[name]]
    if (!methods::is(annotation, "GRangesList") || is.null(names(annotation)) ||
        anyDuplicated(names(annotation)) || !all(names(annotation) %in% ids)) {
      stop("Every region must be a named GRangesList using transcript identifiers")
    }
  }
  if (!is.numeric(min_mapq) || length(min_mapq) != 1L ||
      !is.finite(min_mapq) || min_mapq < 0 || min_mapq > 255 ||
      min_mapq != as.integer(min_mapq)) {
    stop("min_mapq must be an integer from zero to 255")
  }
  if (!is.character(region_priority) || anyDuplicated(region_priority) ||
      !all(names(regions) %in% region_priority)) {
    stop("region_priority must contain every supplied region exactly once")
  }
  if (!is.numeric(pseudocount) || length(pseudocount) != 1L ||
      !is.finite(pseudocount) || pseudocount < 0) {
    stop("pseudocount must be one finite non-negative number")
  }
}

# See build_transcript_models() in fragment_geometry.R for the shared
# implementation; this wrapper additionally rejects transcripts whose exons
# span more than one chromosome, which the geometry-learning callers need.
transcript_geometry_models <- function(transcripts) {
  build_transcript_models(transcripts, check_chromosome = TRUE)
}

project_region_learning_reads <- function(reads, models, distribution) {
  reads <- data.table::copy(reads[fragment_length %in% distribution$fragment_length])
  if (!nrow(reads)) stop("No BAM reads have a supported fragment length")
  reads[, read_id := .I]
  points <- GenomicRanges::GRanges(
    reads$chromosome, IRanges::IRanges(reads$five_position, width = 1L),
    strand = reads$strand
  )
  exons <- GenomicRanges::GRangesList(lapply(models, `[[`, "exons"))
  names(exons) <- names(models)
  hits <- suppressWarnings(GenomicRanges::findOverlaps(points, exons))
  candidates <- unique(data.table::data.table(
    read_id = S4Vectors::queryHits(hits),
    model_id = S4Vectors::subjectHits(hits)
  ))
  mapped <- data.table::rbindlist(lapply(
    split(candidates, candidates$model_id), function(group) {
      model <- models[[group$model_id[1L]]]
      result <- data.table::copy(reads[group$read_id])
      result[, tx_start := vapply(
        five_position, function(position) genomic_site_to_transcript(model, position),
        integer(1)
      )]
      result <- result[!is.na(tx_start) &
        tx_start + fragment_length - 1L <= model$length]
      expected <- lapply(seq_len(nrow(result)), function(i) {
        transcript_fragment_alignment(model, result$tx_start[i],
                                      result$fragment_length[i])
      })
      compatible <- result$position == vapply(expected, `[[`, integer(1), "pos") &
        result$cigar == vapply(expected, `[[`, character(1), "cigar")
      result <- result[compatible]
      result[, transcript_id := model$transcript_id]
      result[, site_offset := distribution$site_offset[
        match(fragment_length, distribution$fragment_length)
      ]]
      result[, site_tx := tx_start + site_offset]
      result
    }
  ), use.names = TRUE, fill = TRUE)
  unique_mapping <- mapped[, .N, by = read_id][N == 1L, read_id]
  retained <- mapped[read_id %in% unique_mapping]
  list(
    reads = retained,
    diagnostics = c(
      supported_length = sum(reads$count),
      uniquely_projected = sum(retained$count),
      unmatched_or_ambiguous_transcript = sum(reads$count) - sum(retained$count)
    )
  )
}

transcript_region_intervals <- function(regions, models) {
  data.table::rbindlist(lapply(names(regions), function(region) {
    annotation <- regions[[region]]
    data.table::rbindlist(lapply(names(annotation), function(id) {
      model <- models[[id]]
      ranges <- annotation[[id]]
      if (!length(ranges)) return(NULL)
      starts <- vapply(seq_along(ranges), function(i) {
        genomic_site_to_transcript(model, GenomicRanges::start(ranges)[i])
      }, integer(1))
      ends <- vapply(seq_along(ranges), function(i) {
        genomic_site_to_transcript(model, GenomicRanges::end(ranges)[i])
      }, integer(1))
      keep <- !is.na(starts) & !is.na(ends)
      data.table::data.table(
        transcript_id = id, region = region,
        start_tx = pmin(starts[keep], ends[keep]),
        end_tx = pmax(starts[keep], ends[keep])
      )
    }))
  }))
}

assign_projected_regions <- function(
    reads, intervals, region_names, overlap_action, region_priority) {
  candidates <- merge(reads, intervals, by = "transcript_id", allow.cartesian = TRUE)
  candidates <- unique(candidates[site_tx >= start_tx & site_tx <= end_tx],
                       by = c("read_id", "region"))
  candidates[, overlapping_regions := .N, by = read_id]
  overlap_reads <- unique(
    candidates[, .(read_id, overlapping_regions, count)], by = "read_id"
  )
  overlaps <- overlap_reads[, .(
    alignments = .N, reads = sum(count)
  ), by = overlapping_regions][order(overlapping_regions)]
  if (overlap_action == "priority") {
    candidates[, priority := match(region, region_priority)]
    data.table::setorder(candidates, read_id, priority)
    assigned <- candidates[, .SD[1L], by = read_id]
    assigned[, assigned_count := count]
  } else if (overlap_action == "fractional") {
    assigned <- candidates
    assigned[, assigned_count := count / overlapping_regions]
  } else {
    assigned <- candidates[overlapping_regions == 1L]
    assigned[, assigned_count := count]
  }
  assigned_ids <- unique(assigned$read_id)
  list(
    reads = assigned,
    overlaps = overlaps,
    diagnostics = c(
      projected = sum(reads$count),
      assigned = sum(assigned$assigned_count),
      no_region_or_overlap_excluded = sum(reads[!read_id %in% assigned_ids, count])
    )
  )
}

summarize_region_assignments <- function(
    assigned, transcript_ids, region_names, pseudocount) {
  counts <- assigned[, .(count = sum(assigned_count)),
                     by = .(transcript_id, region)]
  grid <- data.table::CJ(
    transcript_id = transcript_ids, region = region_names, unique = TRUE
  )
  per_transcript <- merge(grid, counts, by = c("transcript_id", "region"),
                          all.x = TRUE, sort = FALSE)
  per_transcript[is.na(count), count := 0]
  per_transcript[, proportion := (count + pseudocount) /
    sum(count + pseudocount), by = transcript_id]
  global <- per_transcript[, .(count = sum(count)), by = region]
  global[, proportion := (count + pseudocount) / sum(count + pseudocount)]
  global[, region_order := match(region, region_names)]
  data.table::setorder(global, region_order)
  global[, region_order := NULL]
  list(global = global, per_transcript = per_transcript)
}

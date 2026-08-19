#' Convert alignments to strand-aware focal-site points
#'
#' @param reads a GRanges or alignment-like object with genomic read ranges.
#' @param focal_offset integer offset from read start on plus-strand reads and
#'   from read end on minus-strand reads.
#' @return a GRanges with width 1.
point_reads <- function(reads, focal_offset = 0L) {
  point <- GenomicRanges::granges(reads)
  plus <- as.character(GenomicRanges::strand(point)) != "-"

  # Ribo-seq footprints are imported as full genomic ranges. For codon-bias
  # learning we count one focal position per read, using the strand-aware 5' end
  # as the anchor and then applying an optional A-/P-site offset.
  plus_points <- GenomicRanges::resize(point[plus], width = 1, fix = "start")
  minus_points <- GenomicRanges::resize(point[!plus], width = 1, fix = "end")

  if (length(plus_points) > 0L && focal_offset != 0L) {
    plus_points <- suppressWarnings(
      GenomicRanges::shift(plus_points, shift = focal_offset)
    )
  }
  if (length(minus_points) > 0L && focal_offset != 0L) {
    minus_points <- suppressWarnings(
      GenomicRanges::shift(minus_points, shift = -focal_offset)
    )
  }

  keep_in_bounds <- function(x) {
    sequence_length <- GenomeInfoDb::seqlengths(x)[
      as.character(GenomicRanges::seqnames(x))
    ]
    keep <- GenomicRanges::start(x) >= 1L &
      (is.na(sequence_length) | GenomicRanges::end(x) <= sequence_length)
    x[keep]
  }

  c(keep_in_bounds(plus_points), keep_in_bounds(minus_points))
}

#' Learn empirical read lengths from imported reads
#'
#' @param reads a GRanges or alignment-like object with read widths.
#' @param min_length minimum read length to keep.
#' @param max_length maximum read length to keep.
#' @param max_observations optional cap for the returned empirical vector.
#' @return an integer vector. Repeated values preserve the empirical frequency
#'   distribution when passed to \code{simNGScoverage(read_lengths_per=)}.
learn_read_lengths <- function(reads, min_length = 25L, max_length = 34L,
                               max_observations = 1e6) {
  widths <- suppressWarnings(as.integer(ORFik::readWidths(reads)))

  # Keep only the footprint-size range we want the simulator to sample from.
  # Returning repeated lengths preserves the empirical frequency distribution.
  widths <- widths[widths >= min_length & widths <= max_length]
  if (!length(widths)) {
    stop("No read lengths remained after filtering.", call. = FALSE)
  }

  if (!is.null(max_observations) && length(widths) > max_observations) {
    widths <- sample(widths, size = max_observations, replace = FALSE)
  }

  widths
}

expected_site_offset <- function(fragment_length,
                                 site_reference = c("a_site", "p_site")) {
  site_reference <- match.arg(site_reference)
  fragment_length <- as.integer(fragment_length)
  p_offset <- ifelse(
    fragment_length <= 27L, 11L,
    ifelse(fragment_length <= 30L, 12L, 13L)
  )
  as.integer(p_offset + ifelse(site_reference == "a_site", 3L, 0L))
}

point_reads_from_geometry <- function(reads, geometry_distribution) {
  geometry <- data.table::as.data.table(geometry_distribution)
  required <- c("fragment_length", "site_offset")
  if (!all(required %in% names(geometry))) {
    stop("geometry_distribution requires fragment_length and site_offset")
  }
  if (anyDuplicated(geometry$fragment_length)) {
    stop("Coverage learning requires one learned offset per fragment length")
  }
  widths <- suppressWarnings(as.integer(ORFik::readWidths(reads)))
  points <- lapply(seq_len(nrow(geometry)), function(i) {
    keep <- widths == geometry$fragment_length[i]
    if (!any(keep)) return(NULL)
    point_reads(reads[keep], focal_offset = geometry$site_offset[i])
  })
  points <- points[!vapply(points, is.null, logical(1))]
  if (!length(points)) stop("No reads matched the learned fragment geometry")
  do.call(c, unname(points))
}

site_frame_score <- function(cds, points) {
  mapped <- GenomicFeatures::mapToTranscripts(points, cds)
  total <- length(mapped)
  in_frame <- sum((GenomicRanges::start(mapped) - 1L) %% 3L == 0L)
  c(in_frame = in_frame, total = total,
    frame_fraction = if (total > 0) in_frame / total else NA_real_)
}

learn_fragment_geometry <- function(reads, cds,
                                    site_reference = c("a_site", "p_site"),
                                    min_length = 25L, max_length = 34L,
                                    min_reads_per_length = 100L,
                                    max_observations = 1e6) {
  site_reference <- match.arg(site_reference)
  overlaps <- GenomicRanges::findOverlaps(
    GenomicRanges::granges(reads), unlist(cds, use.names = FALSE),
    ignore.strand = FALSE
  )
  keep <- unique(S4Vectors::queryHits(overlaps))
  if (!length(keep)) stop("No reads overlapped the CDS set for geometry learning")
  reads <- reads[keep]
  widths <- suppressWarnings(as.integer(ORFik::readWidths(reads)))
  eligible <- which(widths >= min_length & widths <= max_length)
  if (!length(eligible)) stop("No reads remained for fragment-geometry learning")
  if (!is.null(max_observations) && length(eligible) > max_observations) {
    eligible <- sample(eligible, as.integer(max_observations), replace = FALSE)
  }
  reads <- reads[eligible]
  widths <- widths[eligible]
  length_counts <- table(widths)
  usable_lengths <- as.integer(names(length_counts)[
    length_counts >= min_reads_per_length
  ])
  if (!length(usable_lengths)) {
    stop("No fragment length reached min_reads_per_length")
  }

  diagnostics <- data.table::rbindlist(lapply(usable_lengths, function(fragment_length) {
    expected <- expected_site_offset(fragment_length, site_reference)
    offsets <- pmax(0L, expected + c(-1L, 0L, 1L))
    length_reads <- reads[widths == fragment_length]
    data.table::rbindlist(lapply(offsets, function(offset) {
      score <- site_frame_score(cds, point_reads(length_reads, offset))
      data.table::data.table(
        fragment_length = fragment_length,
        site_offset = offset,
        expected_offset = expected,
        read_count = length(length_reads),
        in_frame_count = as.integer(score[["in_frame"]]),
        cds_count = as.integer(score[["total"]]),
        frame_fraction = score[["frame_fraction"]]
      )
    }))
  }))
  diagnostics[, selection_score := data.table::fifelse(
    is.na(frame_fraction), -Inf, frame_fraction
  )]
  selected <- diagnostics[
    order(-selection_score, abs(site_offset - expected_offset)),
    .SD[1L], by = fragment_length
  ]
  if (any(!is.finite(selected$selection_score))) {
    stop("No CDS-position evidence was available for one or more fragment lengths")
  }
  selected[, selection_score := NULL]
  diagnostics[, selection_score := NULL]
  selected[, probability := read_count / sum(read_count)]
  selected[, site_reference := site_reference]
  attr(selected, "offset_diagnostics") <- diagnostics
  selected[]
}

#' Learn CDS count table from real reads
#'
#' @param cds CDS GRanges/GRangesList to count over.
#' @param reads imported reads.
#' @param sample_name column name for the generated count table.
#' @param condition condition label.
#' @param replicate replicate label.
#' @param min_reads minimum CDS counts to keep.
#' @return a RangedSummarizedExperiment with \code{gene} and \code{cds} assays.
learn_cds_count_table <- function(cds, reads, sample_name = "RFP_real_1",
                                  condition = "real", replicate = "1",
                                  min_reads = 1L) {
  # coveragePerTiling returns one coverage vector per CDS/transcript. Summing
  # each vector gives a simulator-ready CDS count table for the real sample.
  cov <- ORFik::coveragePerTiling(
    cds,
    reads,
    is.sorted = TRUE,
    as.data.table = FALSE
  )
  counts <- vapply(cov, function(x) sum(as.numeric(x)), numeric(1))
  keep <- counts >= min_reads
  if (!any(keep)) {
    stop("No CDS regions passed min_reads.", call. = FALSE)
  }

  cds <- cds[keep]
  counts <- as.integer(round(counts[keep]))
  count_mat <- matrix(counts, ncol = 1)
  rownames(count_mat) <- names(cds)
  colnames(count_mat) <- sample_name

  col_data <- S4Vectors::DataFrame(
    libtype = factor("RFP"),
    condition = factor(condition),
    replicate = factor(replicate)
  )
  rownames(col_data) <- sample_name

  SummarizedExperiment::SummarizedExperiment(
    assays = list(gene = count_mat, cds = count_mat),
    rowRanges = cds,
    colData = col_data
  )
}

#' Learn codon-bias alpha weights from real focal-site coverage
#'
#' @param cds CDS GRanges/GRangesList.
#' @param reads imported reads.
#' @param fa_file genome FASTA matching \code{cds}.
#' @param focal_offset integer offset for focal-site conversion.
#' @param geometry_distribution optional learned table with one site offset per
#'   fragment length. When supplied it supersedes `focal_offset`.
#' @param min_tx_reads minimum transcript focal-site counts to include.
#' @param alpha_scale multiplier for learned probabilities.
#' @param pseudocount small count added to each codon class.
#' @return a data.table with columns \code{variable}, \code{seqs}, and
#'   \code{alpha}, compatible with \code{simNGScoverage(seq_bias=)}.
learn_codon_seq_bias <- function(cds, reads, fa_file, focal_offset = 0L,
                                 geometry_distribution = NULL,
                                 min_tx_reads = 20L, alpha_scale = 100,
                                 pseudocount = 1e-3) {
  # First reduce reads to focal-site points, then count those points over CDS
  # codons. This estimates codon-specific enrichment while avoiding full-read
  # footprint width as an extra signal.
  point <- if (is.null(geometry_distribution)) {
    point_reads(reads, focal_offset = focal_offset)
  } else {
    point_reads_from_geometry(reads, geometry_distribution)
  }
  GenomeInfoDb::seqlevels(point, pruning.mode = "coarse") <-
    GenomeInfoDb::seqlevels(cds)
  GenomeInfoDb::seqlengths(point) <-
    GenomeInfoDb::seqlengths(cds)[GenomeInfoDb::seqlevels(point)]

  cov <- ORFik::coveragePerTiling(
    cds,
    point,
    is.sorted = TRUE,
    as.data.table = FALSE
  )
  cds_seq <- ORFik::txSeqsFromFa(cds, fa_file)
  codons_by_tx <- stringi::stri_extract_all_regex(as.character(cds_seq), ".{3}")
  names(codons_by_tx) <- names(cds_seq)

  all_codons <- names(Biostrings::GENETIC_CODE)
  codon_scores <- stats::setNames(rep(pseudocount, length(all_codons)), all_codons)

  used_tx <- 0L
  for (tx in names(cds)) {
    cov_tx <- as.numeric(cov[[tx]])
    codons <- codons_by_tx[[tx]]
    usable_width <- min(length(cov_tx), length(codons) * 3L)
    if (usable_width < 3L) next

    cov_tx <- cov_tx[seq_len(usable_width)]
    codons <- codons[seq_len(floor(usable_width / 3L))]
    codon_cov <- colSums(matrix(cov_tx[seq_len(length(codons) * 3L)], nrow = 3))
    tx_sum <- sum(codon_cov)
    if (tx_sum < min_tx_reads) next

    tab <- data.table::data.table(codon = codons, count = codon_cov)
    tab <- tab[!grepl("N", codon)]
    if (!nrow(tab)) next

    # Normalize within each transcript before pooling so highly expressed CDSs
    # do not dominate the learned codon bias by count depth alone.
    tab[, occurrences := .N, by = codon]
    tab[, tx_norm := (count / tx_sum) / occurrences]
    tab <- tab[, .(tx_norm = sum(tx_norm)), by = codon]
    codon_scores[tab$codon] <- codon_scores[tab$codon] + tab$tx_norm
    used_tx <- used_tx + 1L
  }

  if (used_tx == 0L) {
    stop("No transcripts had enough focal-site coverage for codon-bias learning.",
         call. = FALSE)
  }

  probs <- codon_scores / sum(codon_scores)
  data.table::data.table(
    variable = "real_input",
    seqs = names(probs),
    alpha = pmax(as.numeric(probs) * alpha_scale, 1e-6)
  )
}

#' Keep the most highly covered CDSs in a learned count table
#'
#' @param count_table a RangedSummarizedExperiment produced by
#'   \code{learn_cds_count_table()}.
#' @param top_n maximum number of CDSs to keep. Use NULL to keep all.
#' @return a subset of \code{count_table}, ordered by decreasing CDS counts.
subset_top_cds_counts <- function(count_table, top_n = NULL) {
  if (is.null(top_n)) return(count_table)
  stopifnot(is.numeric(top_n), length(top_n) == 1L, top_n > 0)

  counts <- as.numeric(SummarizedExperiment::assay(count_table, "cds")[, 1])
  keep_n <- min(as.integer(top_n), length(counts))
  keep <- order(counts, decreasing = TRUE)[seq_len(keep_n)]

  count_table[keep, ]
}

#' Downscale learned counts for memory-safe simulations
#'
#' @param count_table a RangedSummarizedExperiment with \code{gene} and
#'   \code{cds} assays.
#' @param scale numeric multiplier in (0, 1]. Use 1 to keep counts unchanged.
#' @param min_count minimum non-zero count retained after scaling.
#' @return \code{count_table} with scaled integer assays.
scale_count_table <- function(count_table, scale = 1, min_count = 1L) {
  stopifnot(is.numeric(scale), length(scale) == 1L, scale > 0, scale <= 1)
  if (scale == 1) return(count_table)

  for (assay_name in c("gene", "cds")) {
    mat <- SummarizedExperiment::assay(count_table, assay_name)
    scaled <- round(mat * scale)
    scaled[mat > 0 & scaled < min_count] <- min_count
    SummarizedExperiment::assay(count_table, assay_name) <- scaled
  }

  count_table
}

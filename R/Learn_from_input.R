#' Convert alignments to strand-aware focal-site points
#'
#' @param reads a GRanges or alignment-like object with genomic read ranges.
#' @param focal_offset integer offset from read start on plus-strand reads and
#'   from read end on minus-strand reads.
#' @return a GRanges with width 1.
point_reads <- function(reads, focal_offset = 0L) {
  point <- GenomicRanges::granges(reads)
  plus <- as.character(GenomicRanges::strand(point)) != "-"

  plus_points <- GenomicRanges::resize(point[plus], width = 1, fix = "start")
  minus_points <- GenomicRanges::resize(point[!plus], width = 1, fix = "end")

  if (length(plus_points) > 0L && focal_offset != 0L) {
    plus_points <- GenomicRanges::shift(plus_points, shift = focal_offset)
  }
  if (length(minus_points) > 0L && focal_offset != 0L) {
    minus_points <- GenomicRanges::shift(minus_points, shift = -focal_offset)
  }

  c(plus_points, minus_points)
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
  widths <- as.integer(ORFik::readWidths(reads))
  widths <- widths[widths >= min_length & widths <= max_length]
  if (!length(widths)) {
    stop("No read lengths remained after filtering.", call. = FALSE)
  }

  if (!is.null(max_observations) && length(widths) > max_observations) {
    widths <- sample(widths, size = max_observations, replace = FALSE)
  }

  widths
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
#' @param min_tx_reads minimum transcript focal-site counts to include.
#' @param alpha_scale multiplier for learned probabilities.
#' @param pseudocount small count added to each codon class.
#' @return a data.table with columns \code{variable}, \code{seqs}, and
#'   \code{alpha}, compatible with \code{simNGScoverage(seq_bias=)}.
learn_codon_seq_bias <- function(cds, reads, fa_file, focal_offset = 0L,
                                 min_tx_reads = 20L, alpha_scale = 100,
                                 pseudocount = 1e-3) {
  point <- point_reads(reads, focal_offset = focal_offset)
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

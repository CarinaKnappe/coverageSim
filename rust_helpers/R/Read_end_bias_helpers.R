read_end_nt_counts_from_sequences <- function(seqs, run_id = "sample") {
  seqs <- as.character(seqs)
  seqs <- seqs[!is.na(seqs) & nzchar(seqs)]

  nts <- c("A", "C", "G", "T")
  if (!length(seqs)) {
    return(data.table::data.table(
      runID = run_id,
      nt = nts,
      last_nt_count = 0L,
      last_nt_fraction = 0,
      all_read_nt_count = 0L,
      all_read_nt_fraction = 0,
      last_vs_all_enrichment = NA_real_
    ))
  }

  last_nt <- substr(seqs, nchar(seqs), nchar(seqs))
  last_counts <- table(factor(last_nt, levels = nts))

  # Whole-read nucleotide composition is the baseline: how often each base occurs
  # anywhere in the read sequences, independent of read-end position.
  base_counts <- colSums(
    Biostrings::alphabetFrequency(
      Biostrings::DNAStringSet(seqs),
      baseOnly = TRUE
    )[, nts, drop = FALSE]
  )

  last_total <- sum(last_counts)
  base_total <- sum(base_counts)

  out <- data.table::data.table(
    runID = run_id,
    nt = nts,
    last_nt_count = as.integer(last_counts[nts]),
    last_nt_fraction = if (last_total > 0) as.numeric(last_counts[nts]) / last_total else 0,
    all_read_nt_count = as.integer(base_counts[nts]),
    all_read_nt_fraction = if (base_total > 0) as.numeric(base_counts[nts]) / base_total else 0
  )

  out[, last_vs_all_enrichment := data.table::fifelse(
    all_read_nt_fraction > 0,
    last_nt_fraction / all_read_nt_fraction,
    NA_real_
  )]

  out
}


read_sequences_from_reference <- function(chunk, fasta_file) {
  if (is.null(fasta_file) || !file.exists(fasta_file)) {
    stop("A reference FASTA is required when BAM SEQ fields are empty.")
  }

  keep <- !is.na(chunk$rname) & !is.na(chunk$pos) & !is.na(chunk$cigar)
  if (!any(keep)) return(character())

  rname <- as.character(chunk$rname[keep])
  pos <- as.integer(chunk$pos[keep])
  cigar <- as.character(chunk$cigar[keep])
  strand <- as.character(chunk$strand[keep])

  ref_width <- suppressWarnings(
    GenomicAlignments::cigarWidthAlongReferenceSpace(cigar)
  )
  valid <- !is.na(ref_width) & ref_width > 0L
  if (!any(valid)) return(character())

  ranges <- GenomicRanges::GRanges(
    seqnames = rname[valid],
    ranges = IRanges::IRanges(
      start = pos[valid],
      width = ref_width[valid]
    )
  )

  fasta <- Rsamtools::FaFile(fasta_file)
  open(fasta)
  on.exit(close(fasta), add = TRUE)

  seqs <- unname(as.character(Biostrings::getSeq(fasta, ranges)))
  is_minus <- strand[valid] == "-"

  if (any(is_minus)) {
    seqs[is_minus] <- as.character(
      Biostrings::reverseComplement(Biostrings::DNAStringSet(seqs[is_minus]))
    )
  }

  seqs
}

summarize_bam_read_end_nt <- function(bam_file, run_id = basename(bam_file),
                                      yield_size = 1000000L,
                                      fasta_file = NULL) {
  if (!file.exists(bam_file)) stop("BAM file does not exist: ", bam_file)

  nts <- c("A", "C", "G", "T")
  last_counts <- stats::setNames(rep(0, length(nts)), nts)
  base_counts <- stats::setNames(rep(0, length(nts)), nts)
  sequence_source <- "bam_seq"

  bam <- Rsamtools::BamFile(bam_file, yieldSize = yield_size)
  param <- Rsamtools::ScanBamParam(
    what = c("seq", "rname", "pos", "strand", "cigar"),
    flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE)
  )

  open(bam)
  on.exit(close(bam), add = TRUE)

  repeat {
    chunk <- Rsamtools::scanBam(bam, param = param)[[1]]
    if (!length(chunk$seq)) break

    seqs <- as.character(chunk$seq)
    seqs <- seqs[!is.na(seqs) & nzchar(seqs)]

    if (!length(seqs)) {
      seqs <- read_sequences_from_reference(chunk, fasta_file)
      sequence_source <- "reference_reconstructed"
    }

    chunk_dt <- read_end_nt_counts_from_sequences(seqs, run_id = run_id)
    last_counts[chunk_dt$nt] <- last_counts[chunk_dt$nt] + chunk_dt$last_nt_count
    base_counts[chunk_dt$nt] <- base_counts[chunk_dt$nt] + chunk_dt$all_read_nt_count
  }

  last_total <- sum(last_counts)
  base_total <- sum(base_counts)

  out <- data.table::data.table(
    runID = run_id,
    sequence_source = sequence_source,
    nt = nts,
    last_nt_count = as.integer(last_counts[nts]),
    last_nt_fraction = if (last_total > 0) as.numeric(last_counts[nts]) / last_total else 0,
    all_read_nt_count = as.integer(base_counts[nts]),
    all_read_nt_fraction = if (base_total > 0) as.numeric(base_counts[nts]) / base_total else 0
  )

  out[, last_vs_all_enrichment := data.table::fifelse(
    all_read_nt_fraction > 0,
    last_nt_fraction / all_read_nt_fraction,
    NA_real_
  )]

  out
}

summarize_bam_files_read_end_nt <- function(bam_files,
                                            run_ids = basename(bam_files),
                                            yield_size = 1000000L,
                                            fasta_file = NULL) {
  stopifnot(length(bam_files) == length(run_ids))

  data.table::rbindlist(Map(
    function(bam_file, run_id) {
      message("Summarizing read-end nucleotide bias: ", run_id)
      summarize_bam_read_end_nt(
        bam_file = bam_file,
        run_id = run_id,
        yield_size = yield_size,
        fasta_file = fasta_file
      )
    },
    bam_files,
    run_ids
  ))
}

plot_read_end_nt_summary <- function(summary_dt, pdf_file) {
  dir.create(dirname(pdf_file), recursive = TRUE, showWarnings = FALSE)

  pdf(pdf_file, width = 8, height = 5)
  on.exit(dev.off(), add = TRUE)

  for (run_id in unique(summary_dt$runID)) {
    x <- summary_dt[runID == run_id]
    mat <- rbind(
      last_nt = x$last_nt_fraction,
      whole_read = x$all_read_nt_fraction
    )
    colnames(mat) <- x$nt

    barplot(
      mat,
      beside = TRUE,
      ylim = c(0, max(mat, na.rm = TRUE) * 1.2),
      ylab = "Fraction",
      main = paste("Read-end nucleotide bias:", run_id),
      col = c("firebrick", "gray70")
    )
    legend(
      "topright",
      legend = c("last nucleotide", "all read nucleotides"),
      fill = c("firebrick", "gray70"),
      bty = "n"
    )
  }

  invisible(pdf_file)
}

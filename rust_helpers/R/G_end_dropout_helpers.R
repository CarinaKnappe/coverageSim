g_end_dropout_label <- function(drop_fraction) {
  paste0("drop_G_", as.integer(round(drop_fraction * 100)))
}

reference_read_sequences_from_alignment <- function(chunk, fasta_file) {
  if (is.null(fasta_file) || !file.exists(fasta_file)) {
    stop("Reference FASTA does not exist: ", fasta_file)
  }

  keep <- !is.na(chunk$rname) & !is.na(chunk$pos) & !is.na(chunk$cigar)
  seqs <- rep(NA_character_, length(chunk$pos))
  if (!any(keep)) return(seqs)

  rname <- as.character(chunk$rname[keep])
  pos <- as.integer(chunk$pos[keep])
  cigar <- as.character(chunk$cigar[keep])
  strand <- as.character(chunk$strand[keep])

  ref_width <- suppressWarnings(
    GenomicAlignments::cigarWidthAlongReferenceSpace(cigar)
  )
  valid <- !is.na(ref_width) & ref_width > 0L
  if (!any(valid)) return(seqs)

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

  reconstructed <- unname(as.character(Biostrings::getSeq(fasta, ranges)))
  is_minus <- strand[valid] == "-"
  if (any(is_minus)) {
    reconstructed[is_minus] <- as.character(
      Biostrings::reverseComplement(Biostrings::DNAStringSet(reconstructed[is_minus]))
    )
  }

  seqs[which(keep)[valid]] <- reconstructed
  seqs
}

reference_terminal_nt_from_alignment <- function(chunk, fasta_file) {
  seqs <- reference_read_sequences_from_alignment(chunk, fasta_file)
  terminal <- rep(NA_character_, length(seqs))
  ok <- !is.na(seqs) & nzchar(seqs)
  terminal[ok] <- substr(seqs[ok], nchar(seqs[ok]), nchar(seqs[ok]))
  terminal
}

make_reference_g_end_filter <- function(fasta_file, drop_fraction, seed) {
  stats_env <- new.env(parent = emptyenv())
  stats_env$total <- 0L
  stats_env$g_end <- 0L
  stats_env$removed <- 0L
  stats_env$kept <- 0L

  set.seed(seed)

  filter_fun <- function(x) {
    terminal <- reference_terminal_nt_from_alignment(x, fasta_file)
    is_g_end <- !is.na(terminal) & terminal == "G"
    remove <- is_g_end & stats::runif(length(is_g_end)) < drop_fraction
    keep <- !remove

    stats_env$total <- stats_env$total + length(keep)
    stats_env$g_end <- stats_env$g_end + sum(is_g_end)
    stats_env$removed <- stats_env$removed + sum(remove)
    stats_env$kept <- stats_env$kept + sum(keep)

    keep
  }

  list(
    rules = S4Vectors::FilterRules(list(reference_no_terminal_G = filter_fun)),
    filter_fun = filter_fun,
    stats = stats_env
  )
}

write_g_end_dropout_stats <- function(stats_env, stats_file) {
  writeLines(
    c(
      paste("total", stats_env$total, sep = "\t"),
      paste("g_end", stats_env$g_end, sep = "\t"),
      paste("removed", stats_env$removed, sep = "\t"),
      paste("kept", stats_env$kept, sep = "\t")
    ),
    stats_file
  )
}

rewrite_experiment_base <- function(exp_lines, input_base, output_base) {
  gsub(input_base, output_base, exp_lines, fixed = TRUE)
}

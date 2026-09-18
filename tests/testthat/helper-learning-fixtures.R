make_end_learning_fixture <- function(seed = 102L) {
  set.seed(seed)
  fasta <- tempfile(fileext = ".fa")
  chromosome <- paste(sample(c("A", "C", "G", "T"), 5000, TRUE), collapse = "")
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(c(chr1 = chromosome)), fasta)
  Rsamtools::indexFa(fasta)
  transcripts <- GenomicRanges::GRangesList(lapply(seq_len(4), function(i) {
    starts <- 1L + (i - 1L) * 1000L + c(0L, 360L)
    if (i %% 2L == 0L) starts <- rev(starts)
    GenomicRanges::GRanges("chr1", IRanges::IRanges(starts, width = 300L),
      strand = if (i %% 2L) "+" else "-", exon_rank = 1:2)
  }))
  names(transcripts) <- paste0("tx", 1:4)
  GenomeInfoDb::seqlengths(transcripts) <- c(chr1 = 5000L)
  models <- transcript_models(transcripts, fasta)
  cds <- GenomicRanges::GRangesList(lapply(models, function(model) {
    exons <- model$exons
    if (model$strand == "+") {
      GenomicRanges::start(exons)[1] <- GenomicRanges::start(exons)[1] + 30L
      GenomicRanges::end(exons)[2] <- GenomicRanges::end(exons)[2] - 30L
    } else {
      GenomicRanges::end(exons)[1] <- GenomicRanges::end(exons)[1] - 30L
      GenomicRanges::start(exons)[2] <- GenomicRanges::start(exons)[2] + 30L
    }
    exons
  }))
  signal <- data.table::rbindlist(lapply(models, function(model) {
    sites <- seq.int(31L, 570L, 3L)
    codons <- substring(model$sequence, sites, sites + 2L)
    data.table::data.table(transcript_id = model$transcript_id,
      signal_position = transcript_position_to_genomic(model, sites),
      score = ifelse(codons == "AAA", 500L, 100L))
  }))
  geometry <- list(source = "user", distribution = data.frame(
    fragment_length = c(28L, 29L), site_offset = c(15L, 16L), probability = c(.5, .5)))
  list(fasta = fasta, transcripts = transcripts, cds = cds, models = models,
       signal = signal, geometry = geometry)
}

make_learning_bam <- function(fixture, five = 4, three = 3) {
  geometry <- fixture$geometry
  geometry$five_prime_bias <- list(source = "user", table =
    make_synthetic_end_bias(enriched_kmer = "G", enriched_weight = five))
  geometry$three_prime_bias <- list(source = "user", table =
    make_synthetic_end_bias(enriched_kmer = "C", enriched_weight = three))
  fragments <- make_simulated_rpf_fragments(fixture$signal, fixture$models, 28:29, geometry)
  bam <- tempfile(fileext = ".bam")
  write_bam_library(simulated_rpf_alignments(fragments,
    GenomeInfoDb::seqinfo(fixture$transcripts)), bam)
  bam
}

fit_learning_fixture <- function(fixture, bam, ...) {
  learn_end_bias(bam, fixture$fasta, fixture$transcripts, fixture$cds,
                 fixture$geometry, ...)
}

relative_end_weight <- function(profile, target, reference = "A") {
  profile$table[kmer == target, weight] / profile$table[kmer == reference, weight]
}

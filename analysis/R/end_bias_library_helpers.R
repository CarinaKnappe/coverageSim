# Local benchmark adapters, not part of the coverageSim package API.
# Model library recovery conditional on a fixed sequencing depth per transcript.
end_library_weights <- function(sequence, five_factor = 1, three_factor = 1) {
  if (any(!is.finite(c(five_factor, three_factor))) ||
      any(c(five_factor, three_factor) <= 0)) stop("Factors must be positive")
  if (anyNA(sequence) || any(!grepl("^[ACGT]+$", sequence))) {
    stop("Expected biological 5-prime to 3-prime DNA sequences")
  }
  first <- substr(sequence, 1L, 1L)
  last <- substr(sequence, nchar(sequence), nchar(sequence))
  ifelse(first == "G", five_factor, 1) * ifelse(last == "C", three_factor, 1)
}

resample_end_library <- function(fragments, five_factor = 1, three_factor = 1,
                                 reads_per_tx = 5000L, seed = 42L) {
  stopifnot(reads_per_tx > 0, reads_per_tx == as.integer(reads_per_tx))
  result <- data.table::copy(data.table::as.data.table(fragments))
  stopifnot(all(is.finite(result$score)), all(result$score >= 0))
  result[, source_score := score]
  result[, library_weight := end_library_weights(sequence, five_factor, three_factor)]
  set.seed(seed)
  result[, score := {
    probability <- source_score * library_weight
    if (sum(probability) <= 0) stop("Transcript has no selectable fragments")
    as.integer(stats::rmultinom(1L, reads_per_tx, probability))
  }, by = transcript_id]
  result[]
}

transcript_genomic_positions <- function(model) {
  exons <- model$exons
  unlist(lapply(seq_along(exons), function(i) {
    if (model$strand == "+") seq.int(start(exons)[i], end(exons)[i])
    else seq.int(end(exons)[i], start(exons)[i])
  }), use.names = FALSE)
}

project_fragment_truth <- function(truth, models) {
  result <- data.table::copy(data.table::as.data.table(truth))
  result[, c("tx_start", "tx_site") := {
    model <- models[[.BY$transcript_id]]
    if (is.null(model)) stop("Unknown transcript")
    positions <- transcript_genomic_positions(model)
    list(match(five_prime_end, positions), match(signal_position, positions))
  }, by = transcript_id]
  if (anyNA(result$tx_start) || anyNA(result$tx_site)) stop("Non-exonic position")
  if (any(result$tx_start + result$site_offset != result$tx_site)) {
    stop("Projection does not preserve the A-site offset")
  }
  result[, sequence_matches := {
    reference <- models[[.BY$transcript_id]]$sequence
    substring(reference, tx_start, tx_start + fragment_length - 1L) == sequence
  }, by = transcript_id]
  if (!all(result$sequence_matches)) stop("Projected fragment sequence mismatch")
  result[, sequence_matches := NULL]
  result[]
}

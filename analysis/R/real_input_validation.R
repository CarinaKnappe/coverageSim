# Helpers for auditable real-input validation; no analyses run when sourced.

project_validation_reads <- function(reads, models) {
  points <- GenomicRanges::GRanges(reads$chromosome,
    IRanges::IRanges(reads$five_position, width = 1L), strand = reads$strand)
  exons <- GenomicRanges::GRangesList(lapply(models, `[[`, "exons"))
  hits <- suppressWarnings(GenomicRanges::findOverlaps(points, exons))
  candidates <- unique(data.table::data.table(
    read_id = S4Vectors::queryHits(hits), model_id = S4Vectors::subjectHits(hits)))
  mapped <- data.table::rbindlist(lapply(split(candidates, candidates$model_id), function(group) {
    model <- models[[group$model_id[1]]]
    r <- data.table::copy(reads[group$read_id])
    r[, read_id := group$read_id]
    r[, tx_start := vapply(five_position, function(p)
      coverageSim:::genomic_site_to_transcript(model, p), integer(1))]
    r <- r[!is.na(tx_start) & tx_start + fragment_length - 1L <= model$length]
    expected <- lapply(seq_len(nrow(r)), function(i)
      coverageSim:::transcript_fragment_alignment(model, r$tx_start[i], r$fragment_length[i]))
    compatible <- r$position == vapply(expected, `[[`, integer(1), "pos") &
      r$cigar == vapply(expected, `[[`, character(1), "cigar")
    r <- r[compatible]
    r[, transcript_id := model$transcript_id]
    r
  }))
  # Exclude ambiguous assignments rather than counting one read more than once.
  mapped[, if (.N == 1L) .SD, by = read_id]
}

calibrate_validation_offsets <- function(mapped, models, cds, min_reads = 1000) {
  starts <- vapply(models, function(model)
    min(coverageSim:::learning_cds_positions(model, cds[[model$transcript_id]])), integer(1))
  ends <- starts + ORFik::widthPerGroup(cds[names(starts)], FALSE) - 1L
  diagnostics <- data.table::rbindlist(lapply(sort(unique(mapped$fragment_length)), function(len) {
    rows <- mapped[fragment_length == len]
    expected <- if (len <= 27) 14L else if (len <= 30) 15L else 16L
    data.table::rbindlist(lapply(expected + (-1L:1L), function(offset) {
      site <- rows$tx_start + offset
      start <- starts[rows$transcript_id]
      end <- ends[rows$transcript_id]
      keep <- site >= start & site <= end
      frames <- (site[keep] - start[keep]) %% 3L
      counts <- vapply(0:2, function(frame) sum(rows$count[keep][frames == frame]), numeric(1))
      data.table::data.table(fragment_length = len, site_offset = offset,
        expected_offset = expected, total = sum(counts), frame0 = counts[1],
        frame1 = counts[2], frame2 = counts[3], fraction = counts[1] / sum(counts))
    }))
  }))
  chosen <- diagnostics[total >= min_reads & is.finite(fraction)][
    order(-fraction, abs(site_offset - expected_offset)), .SD[1], by = fragment_length]
  if (!nrow(chosen)) stop("Insufficient observations to calibrate offsets")
  chosen[, probability := total / sum(total)]
  list(distribution = chosen[, .(fragment_length, site_offset, probability)],
       diagnostics = diagnostics)
}

validation_opportunities <- function(models, cds, distribution, mapped, k = 1L) {
  opportunities <- coverageSim:::end_learning_opportunities(models, cds, distribution, k)
  # Use a whole-panel ambiguity filter before train/test subdivision.
  keys <- c("chromosome", "position", "strand", "cigar")
  data <- coverageSim:::count_end_learning_reads(opportunities, mapped)$data
  data[]
}

end_validation_scores <- function(data, fit) {
  observed <- data.table::copy(data)
  observed[, group := paste(transcript_id, fragment_length, sep = ":")]
  observed <- observed[, if (sum(count) > 0) .SD, by = group]
  five <- fit$five_prime_bias$table
  three <- fit$three_prime_bias$table
  codon <- fit$diagnostics$codon_weights
  w5 <- five$weight[match(observed$five, five$kmer)]
  w3 <- three$weight[match(observed$three, three$kmer)]
  wc <- codon$weight[match(observed$codon, codon$codon)]
  if (anyNA(c(w5, w3, wc))) stop("Missing prediction weights")
  weights <- list(neutral = rep(1, nrow(observed)), codon = wc,
    five_prime = w5, three_prime = w3, both_ends = w5*w3, all = wc*w5*w3)
  data.table::rbindlist(lapply(names(weights), function(name) {
    observed[, weight := weights[[name]]]
    observed[, log_probability := log(weight / sum(weight)), by = group]
    observed[, .(model = name, reads = sum(count),
      log_likelihood = sum(count * log_probability)), by = group]
  }))
}

partitioned_end_validation_scores <- function(opportunities, fit) {
  data.table::rbindlist(lapply(unique(opportunities$partition), function(part) {
    result <- end_validation_scores(opportunities[partition == part], fit)
    result[, partition := part]
    result
  }))
}

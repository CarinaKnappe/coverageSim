#' Learn transferable read-end preferences from aligned Ribo-seq fragments
#'
#' Fits a regularized conditional multinomial model to observed and possible
#' fragments. Transcript-by-length totals are conditioned out; codon effects
#' are fitted jointly as nuisance parameters. These are estimated sequence
#' preferences, not an assumption-free separation of protocol and biology.
#'
#' @param bam Path to a single-end genomic BAM containing complete fragments.
#' @param fasta Reference FASTA used for alignment.
#' @param transcripts Named GRangesList of transcript exons in transcript order.
#' @param cds Named GRangesList of CDS exons, with the same transcript identifiers.
#' @param fragment_geometry Geometry list with an explicit distribution containing
#'   one site_offset per fragment_length. Offsets are zero-based distances from
#'   the biological 5-prime end to the first nucleotide of the reference codon.
#' @param k End motif length, from 1 to 3. Default 1 is the most conservative fit.
#' @param by_length Fit separate end, codon and frame preferences per length
#'   (default FALSE). The pooled codon profile remains available for Step 4;
#'   length-specific residual codon and frame profiles are returned for fragment
#'   selection.
#' @param min_mapq Minimum mapping quality (default 20). NH greater than 1,
#'   secondary, supplementary, failed-QC, paired and unmapped reads are excluded.
#'   Duplicate-flagged reads are retained; missing NH is allowed.
#' @param ridge Positive penalty shrinking log preferences towards zero.
#'   The per-length fragment totals are conditioned on, not learned as unbiased
#'   length probabilities. Geometry probabilities only select supported lengths.
#' @param maxit Maximum optimizer iterations.
#' @param dmn_min_reads Minimum usable reads per transcript for estimating the
#'   Dirichlet-multinomial concentration. Default 50.
#' @param dmn_min_sites Minimum CDS codon opportunities per transcript for
#'   estimating the concentration. Default 20.
#' @param acf_max_lag Maximum codon lag used to learn local coverage
#'   autocorrelation. Default 9.
#' @return A list containing five_prime_bias, three_prime_bias, codon_bias and
#'   frame_bias, ready to place in fragment_geometry, plus sequence_bias and
#'   dmn_alpha_scale for Step 4, diagnostics and provenance. Passing sequence_bias to simNGScoverage() lets
#'   it use the learned dmn_alpha_scale automatically. Save with saveRDS().
#'   Each profile has source='learned' and strength=1. Unobserved motifs shrink
#'   towards neutral; confidence intervals are not provided. Only complete
#'   M/=/X/N alignments unambiguously compatible with one supplied transcript,
#'   with the reference site on a CDS codon boundary, are used. Soft clips,
#'   indels, ambiguous reference bases and boundary-adjusted offsets are excluded.
#' @export
learn_end_bias <- function(bam, fasta, transcripts, cds, fragment_geometry,
                           k = 1L, by_length = FALSE, min_mapq = 20,
                           ridge = 1, maxit = 500L,
                           dmn_min_reads = 50L, dmn_min_sites = 20L,
                           acf_max_lag = 9L) {
  validate_end_learning_input(bam, fasta, transcripts, cds, k, by_length,
                              min_mapq, ridge, maxit,
                              dmn_min_reads, dmn_min_sites, acf_max_lag)
  geometry <- normalize_fragment_geometry(fragment_geometry)
  if (is.null(geometry$distribution)) {
    stop("Learning requires an explicit length/offset distribution")
  }
  distribution <- validate_fragment_distribution(geometry$distribution)
  distribution <- distribution[probability > 0]
  if (anyDuplicated(distribution$fragment_length)) {
    stop("Learning currently requires one offset per fragment length")
  }
  if (any(distribution$fragment_length < k)) stop("Motifs must fit inside fragments")
  models <- transcript_models(transcripts, fasta)
  opportunities <- end_learning_opportunities(models, cds, distribution, k)
  if (!nrow(opportunities)) stop("No valid CDS fragment opportunities")
  observed <- read_end_learning_bam(bam, min_mapq)
  counted <- count_end_learning_reads(opportunities, observed$reads)
  if (!sum(counted$data$count)) stop("No usable reads match the supplied CDS geometry")
  fit <- fit_end_preferences(counted$data, k, by_length, ridge, maxit)
  frame_fit <- if (by_length) {
    learn_length_frame_bias(observed$reads, models, cds, distribution)
  } else {
    list(profile = list(source = "none"), diagnostics = NULL)
  }
  dispersion <- estimate_dmn_alpha_from_opportunities(
    counted$data, fit, dmn_min_reads, dmn_min_sites
  )
  structure <- learn_coverage_structure(dispersion$sites, acf_max_lag)
  fit$dmn_alpha_scale <- dispersion$scale
  fit$auto_correlation <- structure$kernel
  fit$coverage_qc <- structure$qc
  fit$sequence_bias <- learned_sequence_bias(
    fit$diagnostics$codon_weights, dispersion$scale
  )
  fit$codon_bias <- learned_length_codon_bias(fit$diagnostics$codon_weights)
  fit$frame_bias <- frame_fit$profile
  fit$diagnostics$dmn_alpha <- dispersion$diagnostics
  fit$diagnostics$auto_correlation <- structure$acf
  fit$diagnostics$frame_counts <- frame_fit$diagnostics
  fit$diagnostics$reads <- c(observed$diagnostics, counted$diagnostics)
  fit$provenance <- list(
    bam = normalizePath(bam), fasta = normalizePath(fasta),
    transcript_ids = names(transcripts), distribution = distribution,
    site_reference = geometry$site_reference, k = k, by_length = by_length,
    min_mapq = min_mapq, ridge = ridge, maxit = maxit,
    dmn_min_reads = dmn_min_reads, dmn_min_sites = dmn_min_sites,
    acf_max_lag = acf_max_lag,
    model = paste("conditional multinomial: transcript:length + codon +",
                  "5prime + 3prime; robust DMN moment concentration +",
                  "residual codon autocorrelation"),
    created = as.character(Sys.time()), package_version = as.character(utils::packageVersion("coverageSim"))
  )
  class(fit) <- c("covsim_end_bias_fit", "list")
  fit
}

validate_end_learning_input <- function(bam, fasta, transcripts, cds, k,
                                        by_length, min_mapq, ridge, maxit,
                                        dmn_min_reads, dmn_min_sites,
                                        acf_max_lag) {
  for (path in list(bam, fasta)) {
    if (!is.character(path) || length(path) != 1L || !file.exists(path)) {
      stop("bam and fasta must name existing files")
    }
  }
  for (ranges in list(transcripts, cds)) {
    ids <- names(ranges)
    if (!methods::is(ranges, "GRangesList") || !length(ranges) ||
        is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
      stop("transcripts and cds must be nonempty GRangesLists with unique names")
    }
  }
  if (!setequal(names(transcripts), names(cds))) stop("Transcript and CDS names must match")
  scalar <- function(x, lower, upper = Inf, integer = FALSE) {
    is.numeric(x) && length(x) == 1L && is.finite(x) &&
      x >= lower && x <= upper && (!integer || x == floor(x))
  }
  if (!scalar(k, 1, 3, TRUE)) stop("k must be an integer from 1 to 3")
  if (!is.logical(by_length) || length(by_length) != 1L || is.na(by_length)) {
    stop("by_length must be TRUE or FALSE")
  }
  if (!scalar(min_mapq, 0, 255, TRUE)) stop("min_mapq must be between 0 and 255")
  if (!scalar(ridge, .Machine$double.eps)) stop("ridge must be positive and finite")
  if (!scalar(maxit, 1, Inf, TRUE)) stop("maxit must be a positive integer")
  if (!scalar(dmn_min_reads, 1, Inf, TRUE)) {
    stop("dmn_min_reads must be a positive integer")
  }
  if (!scalar(dmn_min_sites, 2, Inf, TRUE)) {
    stop("dmn_min_sites must be an integer of at least two")
  }
  if (!scalar(acf_max_lag, 1, 100, TRUE)) {
    stop("acf_max_lag must be an integer from one to 100")
  }
}

# Require a contiguous, complete CDS in transcript coordinates, even across introns.
learning_cds_positions <- function(model, cds) {
  if (any(as.character(GenomicRanges::seqnames(cds)) !=
          as.character(GenomicRanges::seqnames(model$exons)[1])) ||
      any(as.character(GenomicRanges::strand(cds)) != model$strand)) {
    stop("CDS chromosome/strand does not match transcript ", model$transcript_id)
  }
  positions <- unlist(IRanges::IntegerList(lapply(seq_along(cds), function(i) {
    seq.int(GenomicRanges::start(cds)[i], GenomicRanges::end(cds)[i])
  })), use.names = FALSE)
  mapped <- sort(vapply(positions, function(p) genomic_site_to_transcript(model, p), integer(1)))
  if (length(mapped) != length(positions) || anyNA(mapped) ||
      length(mapped) %% 3L != 0L || any(diff(mapped) != 1L)) {
    stop("CDS must be complete and contiguous within transcript ", model$transcript_id)
  }
  mapped[seq.int(1L, length(mapped), by = 3L)]
}

end_learning_opportunities <- function(models, cds, distribution, k) {
  data.table::rbindlist(lapply(models, function(model) {
    sites <- learning_cds_positions(model, cds[[model$transcript_id]])
    data.table::rbindlist(lapply(seq_len(nrow(distribution)), function(i) {
      length_value <- distribution$fragment_length[i]
      starts <- sites - distribution$site_offset[i]
      valid <- starts >= 1L & starts + length_value - 1L <= model$length
      sites <- sites[valid]
      starts <- starts[valid]
      if (!length(sites)) return(NULL)
      sequences <- substring(model$sequence, starts, starts + length_value - 1L)
      codons <- substring(model$sequence, sites, sites + 2L)
      alignments <- lapply(starts, function(p) transcript_fragment_alignment(model, p, length_value))
      result <- data.table::data.table(
        transcript_id = model$transcript_id, fragment_length = length_value,
        site_tx = sites,
        position = vapply(alignments, `[[`, integer(1), "pos"),
        cigar = vapply(alignments, `[[`, character(1), "cigar"),
        chromosome = as.character(GenomicRanges::seqnames(model$exons)[1]),
        strand = model$strand, codon = codons,
        five = substring(sequences, 1L, k),
        three = substring(sequences, length_value - k + 1L, length_value)
      )
      result[grepl("^[ACGT]+$", sequences) & grepl("^[ACGT]{3}$", codon)]
    }))
  }))
}

dmn_alpha_moment <- function(observed, expected, nt_positions) {
  total <- sum(observed)
  sites <- length(observed)
  if (sites < 2L || total <= 1L || length(expected) != sites ||
      length(nt_positions) != 1L || nt_positions < sites ||
      any(!is.finite(c(observed, expected))) || any(observed < 0) ||
      any(expected <= 0)) {
    return(data.table::data.table(
      reads = total, sites = sites, inflation = NA_real_,
      alpha_total = NA_real_, dmn_alpha_scale = NA_real_, boundary = "invalid"
    ))
  }
  probability <- expected / sum(expected)
  pearson <- sum((observed - total * probability)^2 /
                   (total * probability))
  inflation <- pearson / (sites - 1L)
  boundary <- "interior"
  if (inflation <= 1) {
    alpha_total <- Inf
    scale <- 1e6
    boundary <- "multinomial_limit"
  } else if (inflation >= total) {
    alpha_total <- 0
    scale <- 1e-8
    boundary <- "maximum_overdispersion"
  } else {
    alpha_total <- (total - inflation) / (inflation - 1)
    scale <- alpha_total / nt_positions
  }
  data.table::data.table(
    reads = total, sites = sites, inflation = inflation,
    alpha_total = alpha_total, dmn_alpha_scale = scale, boundary = boundary
  )
}

end_profile_weights <- function(profile, motif, fragment_length) {
  table <- profile$table
  if ("fragment_length" %in% names(table)) {
    key <- paste(table$kmer, table$fragment_length, sep = ":")
    index <- match(paste(motif, fragment_length, sep = ":"), key)
  } else {
    index <- match(motif, table$kmer)
  }
  result <- table$weight[index]
  if (anyNA(result)) stop("Learned end profile is missing an observed motif")
  result
}

estimate_dmn_alpha_from_opportunities <- function(opportunities, fit,
                                                   min_reads, min_sites) {
  data <- data.table::copy(opportunities)
  codon_weight <- codon_profile_weights(
    fit$diagnostics$codon_weights, data$codon, data$fragment_length
  )
  data[, expected_weight :=
         codon_weight *
         end_profile_weights(fit$five_prime_bias, five, fragment_length) *
         end_profile_weights(fit$three_prime_bias, three, fragment_length)]
  if (any(!is.finite(data$expected_weight) | data$expected_weight <= 0)) {
    stop("Fitted coverage weights must be finite and positive")
  }
  data[, group_reads := sum(count), by = .(transcript_id, fragment_length)]
  data <- data[group_reads > 0]
  data[, conditional_probability := expected_weight / sum(expected_weight),
       by = .(transcript_id, fragment_length)]
  sites <- data[, .(
    observed = sum(count),
    expected = sum(group_reads * conditional_probability)
  ), by = .(transcript_id, site_tx)]
  diagnostics <- sites[, {
    result <- dmn_alpha_moment(observed, expected, nt_positions = 3L * .N)
    result[, usable := reads >= min_reads & sites >= min_sites &
             is.finite(dmn_alpha_scale)]
    result
  }, by = transcript_id]
  usable <- diagnostics[usable == TRUE, dmn_alpha_scale]
  if (!length(usable)) {
    warning("DMN concentration could not be estimated; using fallback 1")
    scale <- 1
  } else {
    scale <- stats::median(usable)
  }
  list(scale = scale, diagnostics = diagnostics, sites = sites)
}

learn_coverage_structure <- function(sites, max_lag) {
  autocorrelation <- estimate_residual_autocorrelation(sites, max_lag)
  list(
    kernel = autocorrelation$kernel,
    acf = autocorrelation$diagnostics,
    qc = coverage_roughness_qc(sites)
  )
}

estimate_residual_autocorrelation <- function(sites, max_lag) {
  data <- data.table::copy(sites)
  data.table::setorder(data, transcript_id, site_tx)
  data[, residual := (observed - expected) / sqrt(expected)]
  per_transcript <- data[, {
    n <- .N
    lags <- seq_len(min(max_lag, n - 1L))
    if (!length(lags)) {
      data.table::data.table(lag = integer(), correlation = numeric(),
                             pairs = integer())
    } else data.table::rbindlist(lapply(lags, function(lag) {
      left <- head(residual, -lag)
      right <- tail(residual, -lag)
      # At least three pairs are needed for a standard deviation and a correlation.
      correlation <- if (length(left) > 2L && isTRUE(stats::sd(left) > 0) &&
                         isTRUE(stats::sd(right) > 0)) {
        stats::cor(left, right)
      } else {
        NA_real_
      }
      data.table::data.table(lag = lag, correlation = correlation,
                             pairs = length(left))
    }))
  }, by = transcript_id]
  diagnostics <- per_transcript[is.finite(correlation), .(
    correlation = stats::weighted.mean(correlation, pairs),
    transcripts = .N,
    pairs = sum(pairs)
  ), by = lag]
  diagnostics <- merge(
    data.table::data.table(lag = seq_len(max_lag)), diagnostics,
    by = "lag", all.x = TRUE, sort = TRUE
  )
  diagnostics[is.na(correlation), correlation := 0]
  diagnostics[is.na(transcripts), `:=`(transcripts = 0L, pairs = 0L)]
  neighbour_weights <- pmax(diagnostics$correlation, 0)
  weights <- c(rev(neighbour_weights), 1, neighbour_weights)
  weights <- weights / sum(weights)
  names(weights) <- as.character(seq.int(-max_lag, max_lag))
  class(weights) <- c("covsim_autocorrelation", class(weights))
  list(kernel = weights, diagnostics = diagnostics)
}

coverage_roughness_qc <- function(sites) {
  longest_run <- function(x) {
    runs <- rle(x == 0)
    if (!any(runs$values)) return(0L)
    max(runs$lengths[runs$values])
  }
  data <- data.table::copy(sites)
  data.table::setorder(data, transcript_id, site_tx)
  per_transcript <- data[, .(
    reads = sum(observed),
    sites = .N,
    zero_fraction = mean(observed == 0),
    variance_to_mean = if (mean(observed) > 0) stats::var(observed) / mean(observed) else NA_real_,
    peak_fraction = if (sum(observed) > 0) max(observed) / sum(observed) else NA_real_,
    spike_fraction = mean(observed > stats::qpois(0.999, lambda = expected)),
    longest_zero_run = longest_run(observed)
  ), by = transcript_id]
  measures <- setdiff(names(per_transcript), "transcript_id")
  summary <- data.table::rbindlist(lapply(measures, function(measure) {
    values <- per_transcript[[measure]]
    values <- values[is.finite(values)]
    data.table::data.table(
      measure = measure,
      median = if (length(values)) stats::median(values) else NA_real_,
      q10 = if (length(values)) unname(stats::quantile(values, 0.1)) else NA_real_,
      q90 = if (length(values)) unname(stats::quantile(values, 0.9)) else NA_real_
    )
  }))
  list(per_transcript = per_transcript, summary = summary)
}

learned_sequence_bias <- function(codon_weights, dmn_alpha_scale) {
  weights <- data.table::copy(codon_weights)
  if ("fragment_length" %in% names(weights)) {
    weights[, normalized := weight / mean(weight), by = fragment_length]
    weights <- weights[, .(weight = exp(mean(log(normalized)))), by = codon]
  }
  alpha <- weights$weight / mean(weights$weight)
  data.table::data.table(
    variable = "learned", seqs = weights$codon, alpha = alpha,
    dmn_alpha_scale = dmn_alpha_scale
  )
}

codon_profile_weights <- function(table, codon, fragment_length) {
  if ("fragment_length" %in% names(table)) {
    key <- paste(table$codon, table$fragment_length, sep = ":")
    index <- match(paste(codon, fragment_length, sep = ":"), key)
  } else {
    index <- match(codon, table$codon)
  }
  result <- table$weight[index]
  if (anyNA(result)) stop("Learned codon profile is missing an opportunity")
  result
}

learned_length_codon_bias <- function(codon_weights) {
  if (!"fragment_length" %in% names(codon_weights)) {
    return(list(source = "none"))
  }
  table <- data.table::copy(codon_weights)
  table[, normalized := weight / mean(weight), by = fragment_length]
  pooled <- table[, .(pooled = exp(mean(log(normalized)))), by = codon]
  table <- merge(table, pooled, by = "codon", sort = FALSE)
  table[, weight := normalized / pooled]
  table[, weight := weight / mean(weight), by = fragment_length]
  list(source = "learned", strength = 1,
       table = table[, .(codon, fragment_length, weight)])
}

learn_length_frame_bias <- function(reads, models, cds, distribution,
                                    pseudocount = 0.5) {
  projected <- project_region_learning_reads(reads, models, distribution)$reads
  bounds <- data.table::rbindlist(lapply(models, function(model) {
    sites <- learning_cds_positions(model, cds[[model$transcript_id]])
    data.table::data.table(
      transcript_id = model$transcript_id,
      cds_start = sites[1L], cds_end = sites[length(sites)] + 2L
    )
  }))
  projected <- merge(projected, bounds, by = "transcript_id")
  projected <- projected[site_tx >= cds_start & site_tx <= cds_end]
  projected[, frame := as.integer((site_tx - cds_start) %% 3L)]
  lengths <- sort(unique(distribution$fragment_length))
  grid <- data.table::CJ(fragment_length = lengths, frame = 0:2, unique = TRUE)
  counts <- projected[, .(reads = sum(count)), by = .(fragment_length, frame)]
  counts <- merge(grid, counts, by = c("fragment_length", "frame"), all.x = TRUE)
  counts[is.na(reads), reads := 0]
  counts[, probability := (reads + pseudocount) / sum(reads + pseudocount),
         by = fragment_length]
  counts[, weight := probability / mean(probability), by = fragment_length]
  list(
    profile = list(source = "learned", strength = 1,
      table = counts[, .(frame, fragment_length, weight)]),
    diagnostics = counts
  )
}

# Normalize explicit match/mismatch operations into M for reference compatibility.
learning_alignment_cigar <- function(cigar) {
  distinct <- unique(cigar)
  normalized <- vapply(distinct, function(value) {
    if (!grepl("^([0-9]+[M=XN])+$", value)) return(NA_character_)
    widths <- as.integer(strsplit(value, "[M=XN]")[[1]])
    ops <- strsplit(gsub("[0-9]+", "", value), "")[[1]]
    ops[ops %in% c("=", "X")] <- "M"
    runs <- cumsum(c(TRUE, tail(ops, -1L) != head(ops, -1L)))
    paste0(as.vector(rowsum(widths, runs)), ops[!duplicated(runs)], collapse = "")
  }, character(1), USE.NAMES = FALSE)
  normalized[match(cigar, distinct)]
}

read_end_learning_bam <- function(bam, min_mapq) {
  fields <- c("rname", "pos", "cigar", "strand", "flag", "mapq")
  file <- Rsamtools::BamFile(bam, yieldSize = 500000L)
  open(file)
  on.exit(close(file))
  reads <- data.table::data.table(chromosome = character(), position = integer(),
    five_position = integer(), strand = character(), cigar = character(),
    fragment_length = integer(), count = integer())
  diagnostics <- c(bam_records = 0, flag_or_mapq_or_NH_excluded = 0, unsupported_cigar = 0)
  repeat {
    raw <- Rsamtools::scanBam(file,
      param = Rsamtools::ScanBamParam(what = fields, tag = "NH"))[[1]]
    if (!length(raw$flag)) break
    chunk <- filter_end_learning_reads(raw, min_mapq)
    diagnostics <- diagnostics + chunk$diagnostics
    reads <- data.table::rbindlist(list(reads, chunk$reads))[
      , .(count = sum(count)), by = .(
        chromosome, position, five_position, strand, cigar, fragment_length
      )]
  }
  list(reads = reads, diagnostics = diagnostics)
}

filter_end_learning_reads <- function(raw, min_mapq) {
  keep <- bitwAnd(raw$flag, 1L + 4L + 256L + 512L + 2048L) == 0L & raw$mapq >= min_mapq
  if (!is.null(raw$tag$NH)) keep <- keep & (is.na(raw$tag$NH) | raw$tag$NH <= 1L)
  cigar <- learning_alignment_cigar(raw$cigar[keep])
  reads <- data.table::data.table(chromosome = as.character(raw$rname[keep]),
    position = raw$pos[keep], strand = as.character(raw$strand[keep]), cigar = cigar)
  reads <- reads[!is.na(cigar)]
  if (nrow(reads)) {
    reads[, fragment_length := learning_cigar_width(cigar, reference = FALSE)]
    reads[, reference_width := learning_cigar_width(cigar, reference = TRUE)]
    reads[, five_position := ifelse(
      strand == "+", position, position + reference_width - 1L
    )]
    reads <- reads[, .(count = .N), by = .(
      chromosome, position, five_position, strand, cigar, fragment_length
    )]
  } else {
    reads <- data.table::data.table(
      chromosome = character(), position = integer(), five_position = integer(),
      strand = character(), cigar = character(), fragment_length = integer(),
      count = integer()
    )
  }
  list(reads = reads, diagnostics = c(
    bam_records = length(keep), flag_or_mapq_or_NH_excluded = sum(!keep),
    unsupported_cigar = sum(is.na(cigar))))
}

learning_cigar_width <- function(cigar, reference = FALSE) {
  as.integer(vapply(cigar, function(value) {
    widths <- as.integer(strsplit(value, "[MN]")[[1]])
    operations <- strsplit(gsub("[0-9]+", "", value), "")[[1]]
    sum(widths[if (reference) operations %in% c("M", "N") else operations == "M"])
  }, numeric(1)))
}

count_end_learning_reads <- function(opportunities, reads) {
  keys <- c("chromosome", "position", "strand", "cigar")
  frequencies <- reads[, .(count = sum(count)), by = keys]
  ambiguous <- opportunities[, .N, by = keys][N > 1L, ..keys]
  eligible <- opportunities[!ambiguous, on = keys]
  result <- merge(eligible, frequencies, by = keys, all.x = TRUE, sort = FALSE)
  result[is.na(count), count := 0L]
  used <- sum(result$count)
  list(data = result, diagnostics = c(used = used,
    unmatched_or_ambiguous = sum(reads$count) - used,
    ambiguous_opportunities = nrow(opportunities) - nrow(eligible)))
}

end_motif_levels <- function(k) {
  apply(expand.grid(rep(list(c("A", "C", "G", "T")), k)), 1L, paste0, collapse = "")
}

# Aggregate identical features: zero-count opportunities still contribute exposure.
end_learning_design <- function(opportunities, k, by_length) {
  data <- data.table::copy(opportunities)
  data[, group := paste(transcript_id, fragment_length, sep = ":")]
  totals <- data[, .(total = sum(count)), by = group][total > 0]
  data <- merge(data, totals, by = "group")
  if (!nrow(data)) stop("No groups with usable reads")
  data <- data[, .(count = sum(count), exposure = .N, total = total[1]),
               by = .(group, fragment_length, codon, five, three)]
  lengths <- sort(unique(data$fragment_length))
  motifs <- end_motif_levels(k)
  labels <- if (by_length) as.vector(outer(motifs, lengths, paste, sep = ":")) else motifs
  five <- if (by_length) paste(data$five, data$fragment_length, sep = ":") else data$five
  three <- if (by_length) paste(data$three, data$fragment_length, sep = ":") else data$three
  codons <- end_motif_levels(3L)
  codon_labels <- if (by_length) {
    as.vector(outer(codons, lengths, paste, sep = ":"))
  } else codons
  codon <- if (by_length) {
    paste(data$codon, data$fragment_length, sep = ":")
  } else data$codon
  index <- cbind(match(five, labels), length(labels) + match(three, labels),
                 2L * length(labels) + match(codon, codon_labels))
  list(data = data, index = index, groups = match(data$group, unique(data$group)),
       motifs = motifs, labels = labels, lengths = lengths, codons = codons,
       codon_labels = codon_labels)
}

fit_end_preferences <- function(opportunities, k, by_length, ridge, maxit) {
  design <- end_learning_design(opportunities, k, by_length)
  data <- design$data
  index <- design$index
  n_parameters <- 2L * length(design$labels) + length(design$codon_labels)
  evaluate <- function(beta, gradient = FALSE) {
    eta <- rowSums(matrix(beta[index], nrow = nrow(index))) + log(data$exposure)
    maxima <- as.numeric(tapply(eta, design$groups, max))
    shifted <- exp(eta - maxima[design$groups])
    sums <- as.numeric(rowsum(shifted, design$groups, reorder = FALSE))
    logp <- eta - maxima[design$groups] - log(sums[design$groups])
    if (!gradient) return(-sum(data$count * logp) + ridge * sum(beta^2) / 2)
    residual <- data$total * exp(logp) - data$count
    values <- rowsum(rep(residual, 3L), as.vector(index), reorder = TRUE)
    result <- ridge * beta
    result[as.integer(rownames(values))] <- result[as.integer(rownames(values))] + values[, 1]
    result
  }
  fit <- stats::optim(rep(0, n_parameters), evaluate,
    gr = function(beta) evaluate(beta, TRUE), method = "L-BFGS-B",
    control = list(maxit = maxit, factr = 1e7))
  if (fit$convergence != 0L) stop("End-bias fit did not converge: ", fit$message)
  make_profile <- function(offset) {
    table <- data.table::data.table(kmer = rep(design$motifs,
      if (by_length) length(design$lengths) else 1L),
      weight = exp(fit$par[offset + seq_along(design$labels)]))
    if (by_length) table[, fragment_length := rep(design$lengths, each = length(design$motifs))]
    list(source = "learned", table = table, strength = 1)
  }
  codon_table <- data.table::data.table(
    codon = rep(design$codons,
      if (by_length) length(design$lengths) else 1L),
    weight = exp(tail(fit$par, length(design$codon_labels)))
  )
  if (by_length) {
    codon_table[, fragment_length := rep(
      design$lengths, each = length(design$codons)
    )]
  }
  list(five_prime_bias = make_profile(0L),
       three_prime_bias = make_profile(length(design$labels)),
       diagnostics = list(convergence = fit$convergence, objective = fit$value,
         transcripts = data.table::uniqueN(opportunities[count > 0, transcript_id]),
         groups = max(design$groups), opportunities = sum(data$exposure),
         codon_weights = codon_table,
         motif_counts = opportunities[, .(reads = sum(count), opportunities = .N),
           by = .(fragment_length, five, three)]))
}

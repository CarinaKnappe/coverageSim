normalize_fragment_geometry <- function(fragment_geometry) {
  defaults <- list(
    source = "default",
    site_reference = "a_site",
    distribution = NULL,
    site_offset = NULL,
    boundary_action = "renormalize",
    five_prime_bias = list(source = "none"),
    three_prime_bias = list(source = "none")
  )
  if (is.null(fragment_geometry)) return(defaults)
  stopifnot(is.list(fragment_geometry))
  unknown <- setdiff(names(fragment_geometry), names(defaults))
  if (length(unknown)) {
    stop("Unknown fragment_geometry field(s): ", paste(unknown, collapse = ", "))
  }
  geometry <- utils::modifyList(defaults, fragment_geometry)
  if (!geometry$source %in% c("none", "default", "user", "learned")) {
    stop("fragment_geometry$source must be one of: none, default, user, learned")
  }
  if (!geometry$site_reference %in% c("p_site", "a_site")) {
    stop("fragment_geometry$site_reference must be 'p_site' or 'a_site'")
  }
  if (!is.null(geometry$site_offset) &&
      (!is.numeric(geometry$site_offset) || !length(geometry$site_offset) ||
       anyNA(geometry$site_offset) || any(geometry$site_offset < 0) ||
       any(geometry$site_offset != as.integer(geometry$site_offset)))) {
    stop("fragment_geometry$site_offset must be NULL or non-negative integers")
  }
  if (!geometry$boundary_action %in% c("error", "renormalize", "adjust_offset")) {
    stop("fragment_geometry$boundary_action must be 'error' or 'renormalize'")
  }
  if (identical(geometry$boundary_action, "adjust_offset")) {
    geometry$boundary_action <- "renormalize"
  }
  if (!is.null(geometry$site_offset)) {
    geometry$site_offset <- as.integer(geometry$site_offset)
  }
  for (bias_name in c("five_prime_bias", "three_prime_bias")) {
    bias <- geometry[[bias_name]]
    geometry[[bias_name]] <- normalize_end_bias(bias, bias_name)
  }
  geometry
}

validate_end_bias_table <- function(table, bias_name = "end bias") {
  table <- data.table::as.data.table(table)
  if (!all(c("kmer", "weight") %in% names(table))) {
    stop(bias_name, " table requires columns: kmer, weight")
  }
  if (anyNA(table$kmer) || any(!grepl("^[ACGT]+$", toupper(table$kmer))) ||
      any(!is.finite(table$weight)) || any(table$weight <= 0)) {
    stop(bias_name, " weights must be finite and strictly positive, with DNA k-mers")
  }
  table[, kmer := toupper(as.character(kmer))]
  if ("fragment_length" %in% names(table)) {
    if (anyNA(table$fragment_length) || any(table$fragment_length <= 0 |
                                             table$fragment_length != as.integer(table$fragment_length))) {
      stop(bias_name, " fragment_length must contain positive integers")
    }
    table[, fragment_length := as.integer(fragment_length)]
  } else table[, fragment_length := NA_integer_]
  if (anyDuplicated(table[, .(kmer, fragment_length)])) {
    stop(bias_name, " table contains duplicate kmer/fragment_length rows")
  }
  table[, .(kmer, weight, fragment_length)]
}

normalize_end_bias <- function(bias, bias_name) {
  if (!is.list(bias) || is.null(bias$source) || length(bias$source) != 1L ||
      !bias$source %in% c("none", "default", "user", "learned")) {
    stop("fragment_geometry$", bias_name,
         " must specify source none, default, user, or learned")
  }
  if (bias$source %in% c("none", "default")) {
    return(list(source = bias$source, table = NULL))
  }
  if (is.null(bias$table)) {
    stop("fragment_geometry$", bias_name, " requires a kmer/weight table")
  }
  list(source = bias$source, table = validate_end_bias_table(bias$table, bias_name))
}

end_bias_weight <- function(sequence, fragment_length,
                            five_bias, three_bias) {
  lookup <- function(kmer, length_value, bias) {
    if (is.null(bias$table)) return(1)
    value <- toupper(kmer)
    exact <- bias$table[kmer == value & fragment_length == length_value, weight]
    if (length(exact)) return(exact[[1L]])
    hit <- bias$table[kmer == value & is.na(fragment_length), weight]
    if (!length(hit)) 1 else hit[[1L]]
  }
  extract_kmer <- function(bias, from_end = FALSE) {
    if (is.null(bias$table)) return(NULL)
    k <- unique(nchar(bias$table$kmer))
    if (length(k) != 1L) stop("Each end-bias table must use one k-mer length")
    if (from_end) substr(sequence, nchar(sequence) - k + 1L, nchar(sequence))
    else substr(sequence, 1L, k)
  }
  lookup(extract_kmer(five_bias), fragment_length, five_bias) *
    lookup(extract_kmer(three_bias, TRUE), fragment_length, three_bias)
}

end_bias_kmer <- function(sequence, bias, from_end = FALSE) {
  if (is.null(bias$table)) return(NA_character_)
  k <- unique(nchar(bias$table$kmer))
  if (length(k) != 1L) stop("Each end-bias table must use one k-mer length")
  if (from_end) substr(sequence, nchar(sequence) - k + 1L, nchar(sequence))
  else substr(sequence, 1L, k)
}

#' Create a simple synthetic sequence-end bias profile.
#'
#' The returned table can be passed as `table` in `five_prime_bias` or
#' `three_prime_bias`. Weights are relative sampling weights; 1 is neutral.
#' @param k integer k-mer length.
#' @param enriched_kmer optional k-mer to enrich.
#' @param enriched_weight weight for `enriched_kmer`.
#' @param depleted_kmer optional k-mer to deplete.
#' @param depleted_weight weight for `depleted_kmer`.
#' @return a data.table with `kmer` and `weight` columns.
#' @export
make_synthetic_end_bias <- function(k = 1L, enriched_kmer = "G",
                                    enriched_weight = 2,
                                    depleted_kmer = NULL,
                                    depleted_weight = 0.5) {
  if (length(k) != 1L || is.na(k) || k < 1 || k != as.integer(k)) {
    stop("k must be one positive integer", call. = FALSE)
  }
  k <- as.integer(k)
  kmers <- apply(expand.grid(rep(list(c("A", "C", "G", "T")), k)), 1L,
                 paste0, collapse = "")
  result <- data.table::data.table(kmer = kmers, weight = 1)
  update <- function(kmer, new_weight) {
    if (is.null(kmer)) return(invisible(NULL))
    if (length(kmer) != 1L || nchar(kmer) != k ||
        !grepl("^[ACGT]+$", toupper(kmer)) || !is.finite(new_weight) || new_weight <= 0) {
      stop("Synthetic end-bias k-mers and weights are invalid", call. = FALSE)
    }
    value <- toupper(kmer)
    result[kmer == value, weight := new_weight]
  }
  update(enriched_kmer, enriched_weight)
  update(depleted_kmer, depleted_weight)
  result[]
}

default_fragment_distribution <- function(fragment_lengths,
                                          site_reference = "a_site") {
  length_counts <- table(as.integer(fragment_lengths))
  lengths <- as.integer(names(length_counts))
  length_weight <- stats::dnorm(lengths, mean = 29, sd = 1.5) *
    as.numeric(length_counts)
  if (!any(length_weight > 0)) length_weight <- as.numeric(length_counts)
  p_offset <- ifelse(lengths <= 27L, 11L, ifelse(lengths <= 30L, 12L, 13L))
  offset <- p_offset + ifelse(site_reference == "a_site", 3L, 0L)
  offset <- pmin(offset, lengths - 1L)
  data.table::data.table(
    fragment_length = lengths,
    site_offset = as.integer(offset),
    probability = length_weight / sum(length_weight)
  )
}

validate_fragment_distribution <- function(distribution) {
  distribution <- data.table::as.data.table(distribution)
  required <- c("fragment_length", "site_offset", "probability")
  if (!all(required %in% names(distribution))) {
    stop("Fragment distribution requires columns: ", paste(required, collapse = ", "))
  }
  distribution <- distribution[, ..required]
  if (!nrow(distribution) || anyNA(distribution) ||
      any(distribution$fragment_length <= 0) ||
      any(distribution$fragment_length != as.integer(distribution$fragment_length)) ||
      any(distribution$site_offset < 0) ||
      any(distribution$site_offset != as.integer(distribution$site_offset)) ||
      any(distribution$site_offset >= distribution$fragment_length) ||
      any(!is.finite(distribution$probability)) ||
      any(distribution$probability < 0) || sum(distribution$probability) <= 0) {
    stop("Fragment distribution contains invalid lengths, offsets, or probabilities")
  }
  distribution[, `:=`(
    fragment_length = as.integer(fragment_length),
    site_offset = as.integer(site_offset)
  )]
  distribution <- distribution[, .(probability = sum(probability)),
                               by = .(fragment_length, site_offset)]
  distribution[, probability := probability / sum(probability)]
  distribution
}

resolve_fragment_distribution <- function(fragment_lengths, geometry) {
  if (!is.null(geometry$distribution)) {
    if (geometry$source == "none") {
      stop("source = 'none' cannot be combined with a fragment distribution")
    }
    return(validate_fragment_distribution(geometry$distribution))
  }
  if (geometry$source %in% c("user", "learned")) {
    if (is.null(geometry$site_offset)) {
      stop("A user/learned source requires distribution or site_offset")
    }
  }
  if (!is.null(geometry$site_offset)) {
    distribution <- expand.grid(
      fragment_length = as.integer(fragment_lengths),
      site_offset = geometry$site_offset,
      KEEP.OUT.ATTRS = FALSE
    )
    distribution$probability <- 1
    return(validate_fragment_distribution(distribution))
  }
  if (geometry$source == "none") {
    stop("Simulated RPFs require a default, user, or learned geometry source")
  }
  default_fragment_distribution(fragment_lengths, geometry$site_reference)
}

transcript_models <- function(mrna_ranges, fasta_file) {
  sequences <- ORFik::txSeqsFromFa(mrna_ranges, fasta_file)
  models <- lapply(seq_along(mrna_ranges), function(i) {
    exons <- mrna_ranges[[i]]
    exon_rank <- S4Vectors::mcols(exons)$exon_rank
    if (!is.null(exon_rank)) exons <- exons[order(exon_rank)]
    strand_value <- as.character(unique(GenomicRanges::strand(exons)))
    if (length(strand_value) != 1L || !strand_value %in% c("+", "-")) {
      stop("Transcript has missing or inconsistent strand: ", names(mrna_ranges)[i])
    }
    list(
      transcript_id = names(mrna_ranges)[i],
      exons = exons,
      cumulative_start = cumsum(c(1L, head(GenomicRanges::width(exons), -1L))),
      length = sum(GenomicRanges::width(exons)),
      strand = strand_value,
      sequence = as.character(sequences[[i]])
    )
  })
  names(models) <- names(mrna_ranges)
  models
}

genomic_site_to_transcript <- function(model, genomic_position) {
  exon_index <- which(
    GenomicRanges::start(model$exons) <= genomic_position &
      GenomicRanges::end(model$exons) >= genomic_position
  )
  if (length(exon_index) != 1L) return(NA_integer_)
  i <- exon_index[[1]]
  within_exon <- if (model$strand == "+") {
    genomic_position - GenomicRanges::start(model$exons)[i]
  } else {
    GenomicRanges::end(model$exons)[i] - genomic_position
  }
  as.integer(model$cumulative_start[i] + within_exon)
}

transcript_position_to_genomic <- function(model, transcript_position) {
  transcript_position <- pmin(pmax(as.integer(transcript_position), 1L), model$length)
  exon_index <- findInterval(transcript_position, model$cumulative_start)
  within_exon <- transcript_position - model$cumulative_start[exon_index]
  if (model$strand == "+") {
    GenomicRanges::start(model$exons)[exon_index] + within_exon
  } else {
    GenomicRanges::end(model$exons)[exon_index] - within_exon
  }
}

append_rnase_to_simulated_rpf_table <- function(dt_range, rnase_bias,
                                                transcript_models) {
  reach <- floor(length(rnase_bias[["RFP"]]) / 2L)
  if (reach == 0L) return(dt_range)
  groups <- split(dt_range, dt_range$genes, keep.by = TRUE)
  data.table::rbindlist(lapply(groups, function(group) {
    model <- transcript_models[[group$transcript_id[1]]]
    first_tx <- genomic_site_to_transcript(model, group$start[1])
    last_tx <- genomic_site_to_transcript(model, group$start[nrow(group)])
    prefix <- group[rep(1L, reach)]
    suffix <- group[rep(nrow(group), reach)]
    prefix[, start := transcript_position_to_genomic(
      model, seq.int(first_tx - reach, first_tx - 1L)
    )]
    suffix[, start := transcript_position_to_genomic(
      model, seq.int(last_tx + 1L, last_tx + reach)
    )]
    prefix[, end := start]
    suffix[, end := start]
    data.table::rbindlist(list(prefix, group, suffix))
  }))
}

transcript_fragment_alignment <- function(model, transcript_start,
                                          fragment_length) {
  transcript_end <- transcript_start + fragment_length - 1L
  exon_tx_start <- model$cumulative_start
  exon_tx_end <- exon_tx_start + GenomicRanges::width(model$exons) - 1L
  used <- which(exon_tx_start <= transcript_end & exon_tx_end >= transcript_start)

  pieces <- lapply(used, function(i) {
    overlap_start <- max(transcript_start, exon_tx_start[i])
    overlap_end <- min(transcript_end, exon_tx_end[i])
    offset_start <- overlap_start - exon_tx_start[i]
    offset_end <- overlap_end - exon_tx_start[i]
    if (model$strand == "+") {
      genomic_start <- GenomicRanges::start(model$exons)[i] + offset_start
      genomic_end <- GenomicRanges::start(model$exons)[i] + offset_end
    } else {
      genomic_start <- GenomicRanges::end(model$exons)[i] - offset_end
      genomic_end <- GenomicRanges::end(model$exons)[i] - offset_start
    }
    c(start = genomic_start, end = genomic_end)
  })
  pieces <- do.call(rbind, pieces)
  pieces <- pieces[order(pieces[, "start"]), , drop = FALSE]
  widths <- pieces[, "end"] - pieces[, "start"] + 1L
  if (length(widths) == 1L) {
    cigar <- paste0(widths, "M")
  } else {
    introns <- pieces[-1L, "start"] -
      pieces[-nrow(pieces), "end"] - 1L
    cigar <- paste0(
      paste0(widths, "M"),
      c(paste0(introns, "N"), ""),
      collapse = ""
    )
  }
  list(
    pos = as.integer(min(pieces[, "start"])),
    end = as.integer(max(pieces[, "end"])),
    cigar = cigar
  )
}

make_simulated_rpf_fragments <- function(signal_table, transcript_models,
                                         fragment_lengths, fragment_geometry) {
  geometry <- normalize_fragment_geometry(fragment_geometry)
  candidates <- resolve_fragment_distribution(fragment_lengths, geometry)
  if (any(signal_table$score < 0 | signal_table$score != as.integer(signal_table$score))) {
    stop("Simulated RPF scores must be non-negative integers")
  }

  rows <- lapply(seq_len(nrow(signal_table)), function(i) {
    row <- signal_table[i]
    model <- transcript_models[[row$transcript_id]]
    if (is.null(model)) stop("No transcript model for: ", row$transcript_id)
    site_tx <- genomic_site_to_transcript(model, row$signal_position)
    if (is.na(site_tx)) {
      stop("Signal position is not exonic in transcript ", row$transcript_id)
    }
    valid <- candidates[
      site_tx - candidates$site_offset >= 1L &
        site_tx - candidates$site_offset + candidates$fragment_length - 1L <= model$length,
    ]
    if (!nrow(valid) && geometry$boundary_action == "renormalize") {
      valid <- candidates[fragment_length <= model$length]
      if (nrow(valid)) {
        valid[, site_offset := vapply(seq_len(.N), function(j) {
          lower <- max(0L, site_tx + fragment_length[j] - 1L - model$length)
          upper <- min(fragment_length[j] - 1L, site_tx - 1L)
          min(max(site_offset[j], lower), upper)
        }, integer(1))]
        valid <- valid[, .(probability = sum(probability)),
                       by = .(fragment_length, site_offset)]
      }
    }
    if (!nrow(valid)) {
      stop(
        "Cannot create a complete fragment for transcript ", row$transcript_id,
        " at signal position ", row$signal_position,
        "; affected count: ", row$score
      )
    }
    valid[, `:=`(
      tx_start = site_tx - site_offset,
      tx_end = site_tx - site_offset + fragment_length - 1L
    )]
    valid[, sequence := vapply(seq_len(.N), function(j) {
      substr(model$sequence, tx_start[j], tx_end[j])
    }, character(1))]
    valid[, bias_weight := vapply(seq_len(.N), function(j) {
      end_bias_weight(
        sequence[j],
        fragment_length[j], geometry$five_prime_bias, geometry$three_prime_bias
      )
    }, numeric(1))]
    valid[, probability := probability / sum(probability)]
    valid[, geometry_probability := probability]
    valid[, probability := probability * bias_weight]
    if (sum(valid$probability) <= 0) {
      stop("End-bias weights exclude every feasible fragment", call. = FALSE)
    }
    valid[, probability := probability / sum(probability)]
    allocation <- as.integer(stats::rmultinom(
      1L, size = as.integer(row$score), prob = valid$probability
    ))
    selected <- valid[allocation > 0L]
    selected[, score := allocation[allocation > 0L]]

    data.table::rbindlist(lapply(seq_len(nrow(selected)), function(j) {
      fragment_length <- selected$fragment_length[j]
      site_offset <- selected$site_offset[j]
      tx_start <- selected$tx_start[j]
      alignment <- transcript_fragment_alignment(
        model, tx_start, fragment_length
      )
      tx_end <- tx_start + fragment_length - 1L
      five_prime <- if (model$strand == "+") alignment$pos else alignment$end
      three_prime <- if (model$strand == "+") alignment$end else alignment$pos
      data.table::data.table(
        seqnames = as.character(GenomicRanges::seqnames(model$exons)[1]),
        start = alignment$pos,
        end = alignment$end,
        strand = model$strand,
        cigar = alignment$cigar,
        score = selected$score[j],
        transcript_id = row$transcript_id,
        ribosome_site = row$signal_position,
        signal_position = row$signal_position,
        fragment_start = alignment$pos,
        fragment_end = alignment$end,
        five_prime_end = five_prime,
        three_prime_end = three_prime,
        fragment_length = fragment_length,
        site_offset = site_offset,
        site_reference = geometry$site_reference,
        geometry_probability = selected$geometry_probability[j],
        sequence = selected$sequence[j],
        five_prime_kmer = end_bias_kmer(selected$sequence[j], geometry$five_prime_bias),
        three_prime_kmer = end_bias_kmer(selected$sequence[j], geometry$three_prime_bias, TRUE),
        end_bias_weight = selected$bias_weight[j],
        final_probability = selected$probability[j]
      )
    }))
  })
  result <- data.table::rbindlist(rows)
  result[, fragment_id := sprintf("fragment_%d", seq_len(.N))]
  result
}

simulated_rpf_alignments <- function(fragment_table, seqinfo) {
  GenomicAlignments::GAlignments(
    seqnames = fragment_table$seqnames,
    pos = fragment_table$start,
    cigar = fragment_table$cigar,
    strand = fragment_table$strand,
    seqinfo = seqinfo,
    score = fragment_table$score,
    fragment_id = fragment_table$fragment_id,
    transcript_id = fragment_table$transcript_id,
    ribosome_site = fragment_table$ribosome_site,
    signal_position = fragment_table$signal_position,
    fragment_start = fragment_table$fragment_start,
    fragment_end = fragment_table$fragment_end,
    five_prime_end = fragment_table$five_prime_end,
    three_prime_end = fragment_table$three_prime_end,
    fragment_length = fragment_table$fragment_length,
    site_offset = fragment_table$site_offset,
    site_reference = fragment_table$site_reference,
    geometry_probability = fragment_table$geometry_probability,
    sequence = fragment_table$sequence
  )
}

# Compatibility aliases for scripts written before the simulated-RPF naming.
make_physical_fragments <- function(...) {
  warning("make_physical_fragments() is deprecated; use make_simulated_rpf_fragments()")
  make_simulated_rpf_fragments(...)
}

physical_fragment_alignments <- function(...) {
  warning("physical_fragment_alignments() is deprecated; use simulated_rpf_alignments()")
  simulated_rpf_alignments(...)
}

append_rnase_to_physical_table <- function(...) {
  warning(paste0(
    "append_rnase_to_physical_table() is deprecated; use ",
    "append_rnase_to_simulated_rpf_table()"
  ))
  append_rnase_to_simulated_rpf_table(...)
}

write_fragment_ground_truth <- function(fragment_table, file_base,
                                        ground_truth) {
  if (identical(ground_truth, FALSE) || is.null(ground_truth)) return(NULL)
  path <- if (identical(ground_truth, TRUE)) {
    paste0(file_base, "_ground_truth.tsv")
  } else {
    if (length(ground_truth) != 1L || !is.character(ground_truth)) {
      stop("ground_truth must be FALSE, TRUE, or a directory path")
    }
    dir.create(ground_truth, recursive = TRUE, showWarnings = FALSE)
    file.path(ground_truth, paste0(basename(file_base), "_ground_truth.tsv"))
  }
  columns <- c(
    "fragment_id", "transcript_id", "ribosome_site", "signal_position",
    "fragment_start", "fragment_end", "five_prime_end", "three_prime_end",
    "fragment_length", "site_offset", "site_reference", "geometry_probability",
    "strand", "cigar", "sequence", "score"
  )
  data.table::fwrite(fragment_table[, ..columns], path, sep = "\t")
  path
}

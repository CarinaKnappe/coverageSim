# Fill in the defaults and reject settings that cannot mean anything.
#
# Callers supply only the fields they care about, so this completes the rest and
# checks the result once, up front, rather than letting a bad value surface deep
# inside a loop. Unknown field names are rejected rather than ignored: a
# misspelled field would otherwise be silently dropped and the default used, so
# the simulation would quietly not do what was asked.
normalize_fragment_geometry <- function(fragment_geometry) {
  defaults <- list(
    source = "default",
    site_reference = "a_site",
    distribution = NULL,
    site_offset = NULL,
    boundary_action = "renormalize",
    five_prime_bias = list(source = "none"),
    three_prime_bias = list(source = "none"),
    codon_bias = list(source = "none"),
    frame_bias = list(source = "none")
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
  if (identical(geometry$boundary_action, "adjust_offset")) {
    stop("fragment_geometry$boundary_action 'adjust_offset' no longer exists: ",
         "offsets are not shifted at transcript ends any more. Use 'renormalize' ",
         "(positions without room for a fragment get no reads) or 'error'")
  }
  if (!geometry$boundary_action %in% c("error", "renormalize")) {
    stop("fragment_geometry$boundary_action must be 'error' or 'renormalize'")
  }
  if (!is.null(geometry$site_offset)) {
    geometry$site_offset <- as.integer(geometry$site_offset)
  }
  for (bias_name in c("five_prime_bias", "three_prime_bias")) {
    bias <- geometry[[bias_name]]
    geometry[[bias_name]] <- normalize_end_bias(bias, bias_name)
  }
  geometry$codon_bias <- normalize_fragment_feature_bias(
    geometry$codon_bias, "codon_bias", "codon"
  )
  geometry$frame_bias <- normalize_fragment_feature_bias(
    geometry$frame_bias, "frame_bias", "frame"
  )
  geometry
}

normalize_fragment_feature_bias <- function(bias, bias_name, feature) {
  if (!is.list(bias) || is.null(bias$source) || length(bias$source) != 1L ||
      !bias$source %in% c("none", "default", "user", "learned")) {
    stop("fragment_geometry$", bias_name,
         " must specify source none, default, user, or learned")
  }
  strength <- if (is.null(bias$strength)) 1 else bias$strength
  if (!is.numeric(strength) || length(strength) != 1L ||
      !is.finite(strength) || strength < 0) {
    stop(bias_name, " strength must be one finite non-negative number")
  }
  if (bias$source %in% c("none", "default")) {
    return(list(source = bias$source, table = NULL, strength = strength))
  }
  table <- data.table::copy(data.table::as.data.table(bias$table))
  required <- c(feature, "fragment_length", "weight")
  if (!all(required %in% names(table)) || !nrow(table)) {
    stop(bias_name, " table requires columns: ", paste(required, collapse = ", "))
  }
  if (feature == "codon") {
    table[, codon := toupper(as.character(codon))]
    if (anyNA(table$codon) || any(!grepl("^[ACGT]{3}$", table$codon))) {
      stop(bias_name, " codons must be three-base DNA strings")
    }
  } else {
    if (anyNA(table$frame) || any(!table$frame %in% 0:2)) {
      stop(bias_name, " frames must be 0, 1, or 2")
    }
    table[, frame := as.integer(frame)]
  }
  if (anyNA(table$fragment_length) ||
      any(table$fragment_length <= 0 | !is.finite(table$fragment_length) |
          table$fragment_length != as.integer(table$fragment_length)) ||
      anyNA(table$weight) || any(!is.finite(table$weight) | table$weight <= 0)) {
    stop(bias_name, " lengths and weights must be finite and positive")
  }
  table[, fragment_length := as.integer(fragment_length)]
  if (anyDuplicated(table[, c(feature, "fragment_length"), with = FALSE])) {
    stop(bias_name, " table contains duplicate feature/length rows")
  }
  list(source = bias$source, table = table[, c(feature, "fragment_length", "weight"), with = FALSE],
       strength = strength)
}

fragment_feature_weight <- function(value, fragment_length, bias, feature) {
  if (is.null(bias$table) || identical(bias$strength, 0)) return(1)
  index <- bias$table[[feature]] == value &
    bias$table$fragment_length == fragment_length
  hit <- bias$table$weight[index]
  if (!length(hit)) 1 else hit[[1L]] ^ bias$strength
}

validate_end_bias_table <- function(table, bias_name = "end bias") {
  table <- data.table::copy(data.table::as.data.table(table))
  if (!all(c("kmer", "weight") %in% names(table))) {
    stop(bias_name, " table requires columns: kmer, weight")
  }
  if (anyNA(table$kmer) || any(!grepl("^[ACGT]+$", toupper(table$kmer))) ||
      any(!is.finite(table$weight)) || any(table$weight <= 0)) {
    stop(bias_name, " weights must be finite and strictly positive, with DNA k-mers")
  }
  table[, kmer := toupper(as.character(kmer))]
  if ("fragment_length" %in% names(table)) {
    if (any(!is.na(table$fragment_length) &
            (!is.finite(table$fragment_length) | table$fragment_length <= 0 |
             table$fragment_length != as.integer(table$fragment_length)))) {
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
  strength <- if (is.null(bias$strength)) 1 else bias$strength
  if (!is.numeric(strength) || length(strength) != 1L ||
      !is.finite(strength) || strength < 0) {
    stop(bias_name, " strength must be one finite non-negative number")
  }
  if (bias$source %in% c("none", "default")) {
    return(list(source = bias$source, table = NULL, strength = strength))
  }
  if (is.null(bias$table)) {
    stop("fragment_geometry$", bias_name, " requires a kmer/weight table")
  }
  list(source = bias$source, table = validate_end_bias_table(bias$table, bias_name),
       strength = strength)
}

# How much a fragment's own end sequences change its chance of being sequenced.
#
# Ligation and digestion are not indifferent to sequence, so two fragments from
# the same position can survive library preparation at different rates purely
# because of the bases at their ends. The two ends are looked up independently
# and their weights multiplied, since each end is handled by its own enzymatic
# step.
#
# A table may give weights per fragment length or once for all lengths; the
# length-specific entry wins where both exist. A k-mer the table does not mention
# is treated as neutral rather than as impossible. `strength` scales the whole
# effect as an exponent, so 0 disables the bias and 2 doubles its reach on the
# log scale.
#
# Takes any number of fragments at once and returns one weight each. The lookup
# is a plain vector match rather than a data.table filter per fragment, because
# it runs for every candidate fragment of every read: filtering the table once
# per fragment cost about 4.5 ms, which made a 300-gene sample with an end bias
# take more than twice as long as one without.
end_bias_weight <- function(sequence, fragment_length,
                            five_bias, three_bias) {
  lookup <- function(kmer, bias) {
    if (is.null(bias$table) || identical(bias$strength, 0)) {
      return(rep(1, length(sequence)))
    }
    value <- toupper(kmer)
    table_kmer <- bias$table$kmer
    table_length <- bias$table$fragment_length
    table_weight <- bias$table$weight
    # A weight given for this fragment length wins over one given for all
    # lengths; a k-mer the table does not mention stays neutral.
    by_length <- !is.na(table_length)
    weight <- table_weight[by_length][match(
      paste(value, fragment_length),
      paste(table_kmer[by_length], table_length[by_length])
    )]
    general <- is.na(weight)
    weight[general] <- table_weight[!by_length][
      match(value[general], table_kmer[!by_length])
    ]
    weight[is.na(weight)] <- 1
    weight ^ bias$strength
  }
  extract_kmer <- function(bias, from_end = FALSE) {
    if (is.null(bias$table)) return(NULL)
    k <- unique(nchar(bias$table$kmer))
    if (length(k) != 1L) stop("Each end-bias table must use one k-mer length")
    if (from_end) substr(sequence, nchar(sequence) - k + 1L, nchar(sequence))
    else substr(sequence, 1L, k)
  }
  lookup(extract_kmer(five_bias), five_bias) *
    lookup(extract_kmer(three_bias, TRUE), three_bias)
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
#' Library preparation is not indifferent to sequence: ligation and digestion
#' favour some bases at a fragment's ends over others, so the observed ends are
#' not a faithful record of where the enzyme cut. This builds a deliberately
#' simple version of that effect, for testing whether an analysis notices it.
#'
#' The result is passed as `table` in `five_prime_bias` or `three_prime_bias`.
#' Weights are relative: 1 is neutral, above 1 makes a fragment ending in that
#' k-mer more likely to be kept, below 1 less likely.
#' @param k integer, how many bases at the end the bias depends on. 1 makes it
#'   depend on the final base alone, 2 on the final two, and so on.
#' @param enriched_kmer optional k-mer to favour; every other k-mer stays
#'   neutral.
#' @param enriched_weight how strongly to favour it. 4 means fragments ending
#'   in that k-mer are four times as likely to survive as any other.
#' @param depleted_kmer optional k-mer to disfavour.
#' @param depleted_weight how strongly to disfavour it; below 1.
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

# The biologically plausible default P-site (ribosome peptidyl-site) offset
# from a fragment's 5' end, by fragment length, plus 3 nt more for the
# A-site (one codon further into the ribosome). This is the simulator's own
# default geometry (default_fragment_distribution() below); elsewhere it is
# only an initial guess used to center a search window when learning real
# offsets (learn_fragment_geometry(), and the rust_helpers analysis scripts
# via coverageSim:::default_ribosome_site_offset(), kept in sync with this
# rather than an independently duplicated copy of the same formula).
default_ribosome_site_offset <- function(fragment_length, site_reference) {
  p_offset <- ifelse(fragment_length <= 27L, 11L,
                     ifelse(fragment_length <= 30L, 12L, 13L))
  as.integer(p_offset + ifelse(site_reference == "a_site", 3L, 0L))
}

default_fragment_distribution <- function(fragment_lengths,
                                          site_reference = "a_site") {
  length_counts <- table(as.integer(fragment_lengths))
  lengths <- as.integer(names(length_counts))
  length_weight <- stats::dnorm(lengths, mean = 29, sd = 1.5) *
    as.numeric(length_counts)
  if (!any(length_weight > 0)) length_weight <- as.numeric(length_counts)
  offset <- default_ribosome_site_offset(lengths, site_reference)
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

# Shared core for building the per-transcript position/strand model used by
# fragment-geometry projection and placement, learned end-bias, and
# region-proportion/fragment-geometry learning: exon coordinates, each
# exon's cumulative transcript-position start, transcript length, and
# strand. fasta_file, when given, also fetches and attaches the transcript
# sequence (needed to look up sequence/codon motifs during simulation and
# end-bias learning; the geometry-learning callers don't need it and don't
# pass it). check_chromosome additionally rejects a transcript whose exons
# span more than one chromosome, which only the geometry-learning callers
# have historically required.
# Precompute, once per transcript, everything the fragment code keeps needing.
#
# Placing fragments means repeatedly converting between genome and transcript
# coordinates, which needs the exons in transcript order, where each one starts
# along the transcript, and the total length. Working that out per fragment would
# repeat it millions of times, so it is done once here and carried along.
#
# Exons are sorted by rank rather than by position, because on the minus strand
# the first exon of the transcript has the highest genomic coordinate.
build_transcript_models <- function(ranges, fasta_file = NULL,
                                    check_chromosome = FALSE) {
  sequences <- if (!is.null(fasta_file)) ORFik::txSeqsFromFa(ranges, fasta_file)
  models <- lapply(seq_along(ranges), function(i) {
    exons <- ranges[[i]]
    exon_rank <- S4Vectors::mcols(exons)$exon_rank
    if (!is.null(exon_rank)) exons <- exons[order(exon_rank)]
    strand_value <- as.character(unique(GenomicRanges::strand(exons)))
    strand_ok <- length(strand_value) == 1L && strand_value %in% c("+", "-")
    if (check_chromosome) {
      chromosome <- unique(as.character(GenomicRanges::seqnames(exons)))
      if (!strand_ok || length(chromosome) != 1L) {
        stop("Transcript has inconsistent chromosome or strand: ", names(ranges)[i])
      }
    } else if (!strand_ok) {
      stop("Transcript has missing or inconsistent strand: ", names(ranges)[i])
    }
    model <- list(
      transcript_id = names(ranges)[i],
      exons = exons,
      cumulative_start = cumsum(c(1L, head(GenomicRanges::width(exons), -1L))),
      length = sum(GenomicRanges::width(exons)),
      strand = strand_value
    )
    if (!is.null(fasta_file)) model$sequence <- as.character(sequences[[i]])
    model
  })
  names(models) <- names(ranges)
  models
}

transcript_models <- function(mrna_ranges, fasta_file) {
  # fasta_file is required here (unlike build_transcript_models()'s optional
  # default): every caller needs the attached sequence, so fail fast on a
  # missing one instead of silently building sequence-less models.
  if (is.null(fasta_file)) stop("fasta_file is required by transcript_models()")
  build_transcript_models(mrna_ranges, fasta_file = fasta_file)
}

# How many nucleotides on each side of a footprint the RNase kernel reaches:
# the RFP kernel's half-width, rounded down. Shared by both fragment_mode's
# RNase-extension helpers (this file's append_rnase_to_simulated_rpf_table()
# and Sim_reads_from_counts_helpers.R's append_rnase_to_dt()) so the two
# fragment modes always extend by the same amount for the same kernel.
rnase_kernel_reach <- function(rnase_bias) {
  floor(length(rnase_bias[["RFP"]]) / 2L)
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
  reach <- rnase_kernel_reach(rnase_bias)
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
    if ("region_position" %in% names(group)) {
      prefix[, region_position := seq.int(
        group$region_position[1L] - reach, group$region_position[1L] - 1L
      )]
      suffix[, region_position := seq.int(
        group$region_position[nrow(group)] + 1L,
        group$region_position[nrow(group)] + reach
      )]
    }
    data.table::rbindlist(list(prefix, group, suffix))
  }))
}

# Translate a fragment's position along the transcript into genome coordinates.
#
# A fragment is defined by where it sits in the spliced transcript, but an
# alignment file describes the genome, where the same fragment may be split
# across exons with introns in between. This walks the exons the fragment covers
# and writes the result as a CIGAR string: stretches that match the genome as M,
# the skipped introns as N.
#
# Minus-strand transcripts run against the genome, so their coordinates are
# measured from the exon's far end; the pieces are then sorted into ascending
# genomic order, because that is the order an alignment file expects regardless
# of which strand the gene is on.
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

# Genome coordinates and CIGAR of many fragments at once: the same result as
# transcript_fragment_alignment() for each fragment, without a loop over them.
# `models` are transcript models, `model_index` says which one each fragment
# belongs to, and the fragments are given in transcript coordinates.
fragment_alignments <- function(models, model_index, transcript_start,
                                fragment_length) {
  n <- length(transcript_start)
  if (!n) return(list(pos = integer(0), end = integer(0), cigar = character(0)))
  # All exons of all models in one table, in transcript order. `base` moves
  # each model's transcript coordinates past those of the models before it,
  # so one sorted vector of exon starts serves every fragment.
  exon_count <- vapply(models, function(model) length(model$exons), integer(1))
  model_length <- vapply(models, function(model) as.numeric(model$length), numeric(1))
  base <- c(0, cumsum(model_length))[seq_along(models)]
  exon_model <- rep(seq_along(models), exon_count)
  exon_tx_start <- unlist(lapply(models, function(model) model$cumulative_start))
  exon_width <- unlist(lapply(models, function(model) GenomicRanges::width(model$exons)))
  exon_genome_start <- unlist(lapply(models, function(model) GenomicRanges::start(model$exons)))
  exon_genome_end <- unlist(lapply(models, function(model) GenomicRanges::end(model$exons)))
  exon_minus <- rep(vapply(models, function(model) model$strand != "+", logical(1)),
                    exon_count)
  exon_key <- base[exon_model] + exon_tx_start

  # The first and last exon each fragment touches; it covers the ones between.
  transcript_end <- transcript_start + fragment_length - 1L
  first <- findInterval(base[model_index] + transcript_start, exon_key)
  last <- findInterval(base[model_index] + transcript_end, exon_key)
  piece_count <- last - first + 1L
  fragment <- rep(seq_len(n), piece_count)
  exon <- first[fragment] + sequence(piece_count) - 1L

  # The part of each exon the fragment covers, in genome coordinates.
  offset_start <- pmax(transcript_start[fragment], exon_tx_start[exon]) - exon_tx_start[exon]
  offset_end <- pmin(transcript_end[fragment], exon_tx_start[exon] + exon_width[exon] - 1L) -
    exon_tx_start[exon]
  minus <- exon_minus[exon]
  piece_start <- ifelse(minus, exon_genome_end[exon] - offset_end,
                        exon_genome_start[exon] + offset_start)
  piece_end <- ifelse(minus, exon_genome_end[exon] - offset_start,
                      exon_genome_start[exon] + offset_end)

  # Pieces left to right on the genome; a gap between two is an intron (N).
  order_pieces <- order(fragment, piece_start)
  fragment <- fragment[order_pieces]
  piece_start <- piece_start[order_pieces]
  piece_end <- piece_end[order_pieces]
  first_piece <- !duplicated(fragment)
  last_piece <- !duplicated(fragment, fromLast = TRUE)
  intron <- as.integer(c(piece_start[-1L], NA) - piece_end - 1L)
  piece_text <- paste0(as.integer(piece_end - piece_start + 1L), "M",
                       ifelse(last_piece, "", paste0(intron, "N")))
  cigar <- piece_text[first_piece]
  end <- piece_end[last_piece]
  spliced <- which(piece_count > 1L)
  if (length(spliced)) {
    in_spliced <- fragment %in% spliced
    cigar[spliced] <- vapply(split(piece_text[in_spliced], fragment[in_spliced]),
                             paste, character(1), collapse = "")
    # The rightmost end, also if one exon lies inside another.
    end[spliced] <- vapply(split(piece_end[in_spliced], fragment[in_spliced]),
                           max, numeric(1))
  }
  list(
    pos = as.integer(piece_start[first_piece]),
    end = as.integer(end),
    cigar = cigar
  )
}

has_active_fragment_bias <- function(geometry) {
  bias_names <- c("five_prime_bias", "three_prime_bias", "codon_bias", "frame_bias")
  any(vapply(geometry[bias_names], function(bias) {
    !is.null(bias$table) && bias$strength > 0
  }, logical(1)))
}

# Transcript positions of several genome positions on one transcript; NA where
# a position is not inside exactly one exon. The same arithmetic as
# genomic_site_to_transcript(), for all positions at once.
genomic_sites_to_transcript <- function(model, genomic_positions) {
  exon_start <- GenomicRanges::start(model$exons)
  exon_end <- GenomicRanges::end(model$exons)
  inside <- outer(genomic_positions, exon_start, ">=") &
    outer(genomic_positions, exon_end, "<=")
  exon_index <- max.col(inside, ties.method = "first")
  within_exon <- if (model$strand == "+") {
    genomic_positions - exon_start[exon_index]
  } else {
    exon_end[exon_index] - genomic_positions
  }
  site <- as.integer(model$cumulative_start[exon_index] + within_exon)
  site[rowSums(inside) != 1L] <- NA_integer_
  site
}

# Codon or frame weights for many fragments: each distinct value and length is
# looked up once with fragment_feature_weight().
fragment_feature_weights <- function(values, fragment_lengths, bias, feature) {
  if (is.null(bias$table) || identical(bias$strength, 0)) {
    return(rep(1, length(values)))
  }
  key <- paste(values, fragment_lengths)
  first <- !duplicated(key)
  weights <- vapply(which(first), function(j) {
    fragment_feature_weight(values[j], fragment_lengths[j], bias, feature)
  }, numeric(1))
  weights[match(key, key[first])]
}

# Every fragment each signal position could give rise to, with its weights.
#
# Each position is combined with every length-and-offset pair of `candidates`,
# and the pairs whose fragment fits inside the transcript are kept. Positions
# near a transcript end where nothing fits are handled by `boundary_action`
# (see make_simulated_rpf_fragments()). Rows come in signal-row order, then in
# the order of `candidates`.
fragment_candidates <- function(signal_table, transcript_models, candidates,
                                geometry) {
  n_signal <- nrow(signal_table)
  transcript_id <- signal_table$transcript_id
  used_ids <- unique(transcript_id)
  used_models <- lapply(seq_along(used_ids), function(k) {
    transcript_models[[used_ids[k]]]
  })
  model_index <- match(transcript_id, used_ids)
  has_model <- !vapply(used_models, is.null, logical(1))[model_index]

  site_tx <- rep(NA_integer_, n_signal)
  rows_per_model <- split(seq_len(n_signal), model_index)
  for (k in as.integer(names(rows_per_model))) {
    if (is.null(used_models[[k]])) next
    rows <- rows_per_model[[as.character(k)]]
    site_tx[rows] <- genomic_sites_to_transcript(
      used_models[[k]], signal_table$signal_position[rows]
    )
  }
  transcript_length <- vapply(used_models, function(model) {
    if (is.null(model)) NA_real_ else as.numeric(model$length)
  }, numeric(1))[model_index]

  # Every position against every length-and-offset pair; keep what fits.
  n_pairs <- nrow(candidates)
  signal_row <- rep(seq_len(n_signal), each = n_pairs)
  candidate_index <- rep(seq_len(n_pairs), times = n_signal)
  tx_start <- site_tx[signal_row] - candidates$site_offset[candidate_index]
  fits <- !is.na(tx_start) & tx_start >= 1L &
    tx_start + candidates$fragment_length[candidate_index] - 1L <=
      transcript_length[signal_row]
  pairs <- data.table::data.table(
    signal_row = signal_row[fits],
    fragment_length = candidates$fragment_length[candidate_index[fits]],
    site_offset = candidates$site_offset[candidate_index[fits]],
    probability = candidates$probability[candidate_index[fits]]
  )

  # With "renormalize", a position where no fragment fits simply has no
  # candidates; make_simulated_rpf_fragments() moves or drops its reads.
  nothing_fits <- has_model & !is.na(site_tx) &
    tabulate(pairs$signal_row, nbins = n_signal) == 0L &
    geometry$boundary_action == "error"

  # Stop at the first position that cannot be placed, with the same message as
  # when positions were placed one at a time.
  problem <- which(!has_model | is.na(site_tx) | nothing_fits)
  if (length(problem)) {
    i <- problem[1L]
    if (!has_model[i]) stop("No transcript model for: ", transcript_id[i])
    if (is.na(site_tx[i])) {
      stop("Signal position is not exonic in transcript ", transcript_id[i])
    }
    stop(
      "Cannot create a complete fragment for transcript ", transcript_id[i],
      " at signal position ", signal_table$signal_position[i],
      "; affected count: ", signal_table$score[i]
    )
  }

  transcript_sequence <- vapply(used_models, function(model) model$sequence,
                                character(1))[model_index[pairs$signal_row]]
  pairs[, `:=`(
    tx_start = site_tx[signal_row] - site_offset,
    tx_end = site_tx[signal_row] - site_offset + fragment_length - 1L
  )]
  pairs[, sequence := substr(transcript_sequence, tx_start, tx_end)]
  pairs[, codon := substr(transcript_sequence, site_tx[signal_row],
                          site_tx[signal_row] + 2L)]
  pairs[, frame := if ("region_position" %in% names(signal_table)) {
    as.integer((signal_table$region_position[signal_row] - 1L) %% 3L)
  } else 0L]
  pairs[, end_bias_weight := end_bias_weight(
    sequence, fragment_length,
    geometry$five_prime_bias, geometry$three_prime_bias
  )]
  pairs[, codon_bias_weight := fragment_feature_weights(
    codon, fragment_length, geometry$codon_bias, "codon"
  )]
  pairs[, frame_bias_weight := fragment_feature_weights(
    frame, fragment_length, geometry$frame_bias, "frame"
  )]
  pairs[, bias_weight := end_bias_weight * codon_bias_weight * frame_bias_weight]
  pairs[, probability := probability / sum(probability), by = signal_row]
  pairs[, geometry_probability := probability]
  pairs
}

# Preserve region budgets when present; direct callers supply counts per transcript.
# Decide which of the possible fragments actually get sequenced.
#
# Each signal position could give rise to fragments of several lengths and
# several offsets. Every such candidate is weighted by three things multiplied
# together: how likely that length-and-offset combination is at all, how much the
# sequence at its ends helps or hinders it surviving library preparation, and how
# much signal the position carries. The region's read budget is then drawn over
# the candidates in those proportions, so the total comes out exactly right
# rather than approximately.
#
# `joint` decides what competes with what. Drawing per transcript lets fragments
# compete across the whole transcript, which is what happens in a real library;
# drawing per position keeps each position's count fixed, which is useful when
# the point is to reproduce a given profile exactly.
#
# A group none of whose positions has room for a complete fragment cannot
# receive its reads; they are left out and listed in the "dropped_reads"
# attribute of the result.
allocate_fragment_candidates <- function(candidates, signal_table, joint) {
  groups <- if (joint) fragment_sampling_groups(signal_table) else {
    seq_len(nrow(signal_table))
  }
  # The read budget of each group, from all of its positions.
  budgets <- data.table::data.table(sampling_group = groups,
                                    transcript_id = signal_table$transcript_id)
  if ("read_count" %in% names(signal_table)) {
    budgets[, read_count := signal_table$read_count]
    budgets <- budgets[, .(transcript_id = transcript_id[1],
                           total = if (data.table::uniqueN(read_count) == 1L) {
                             read_count[1]
                           } else NA_real_), by = sampling_group]
  } else {
    budgets[, score := signal_table$score]
    budgets <- budgets[, .(transcript_id = transcript_id[1],
                           total = round(sum(score))), by = sampling_group]
  }
  if (anyNA(budgets$total) || any(budgets$total < 0) ||
      any(budgets$total != floor(budgets$total))) {
    stop("Invalid fragment sampling budget")
  }
  dropped <- budgets[total > 0 & !sampling_group %in% groups[candidates$signal_row],
                     .(sampling_group = as.character(sampling_group),
                       transcript_id = as.character(transcript_id),
                       reads = as.numeric(total))]

  candidates[, sampling_group := groups[signal_row]]
  candidates[, probability := geometry_probability * bias_weight *
               signal_table$score[signal_row]]
  if (any(!is.finite(candidates$probability))) {
    stop("End-bias strength produces non-finite weights; reduce strength")
  }
  candidates[, budget := budgets$total[match(sampling_group, budgets$sampling_group)]]
  # No signal positions at all: there is nothing to draw, and the table stays
  # empty instead of reaching the draw with no budget.
  if (!nrow(candidates)) {
    candidates[, budget := NULL]
    candidates[, score := integer(0)]
  } else {
    candidates[, c("probability", "score") := {
      if (sum(probability) <= 0) stop("End-bias weights exclude every feasible fragment")
      probability <- probability / max(probability)
      p <- probability / sum(probability)
      list(p, as.integer(stats::rmultinom(1L, budget[1], p)))
    }, by = sampling_group]
    candidates[, budget := NULL]
  }
  data.table::setattr(candidates, "dropped_reads", dropped)
  candidates
}

# The groups whose read budgets are kept: the region of a transcript when the
# simulation labelled it (sampling_group), otherwise the whole transcript.
fragment_sampling_groups <- function(signal_table) {
  groups <- if ("sampling_group" %in% names(signal_table)) {
    signal_table$sampling_group
  } else signal_table$transcript_id
  if (anyNA(groups)) stop("sampling_group and transcript_id must not be missing")
  groups
}

# Reads at a position where none of the configured fragment lengths and
# offsets fits inside the transcript.
#
# This happens close to either end of an mRNA, where a ribosome has no room for
# the footprints the distribution allows. Such reads go to the other positions
# of the same region of the same transcript, in proportion to the reads those
# positions already have, so the region keeps its read count. If no other
# position of the region has both reads and room (positions without reads were
# removed before this step), the region loses these reads; they are returned in
# `dropped`. Used when each position keeps its own count (its `read_count` if
# given, else its `score`); when reads are drawn per group (end bias),
# allocate_fragment_candidates() does this by itself.
move_unplaceable_reads <- function(signal_table, placeable) {
  column <- if ("read_count" %in% names(signal_table)) "read_count" else "score"
  count <- signal_table[[column]]
  stuck <- which(!placeable & count > 0)
  if (!length(stuck)) return(list(signal_table = signal_table, dropped = NULL))
  signal_table <- data.table::copy(signal_table)
  groups <- fragment_sampling_groups(signal_table)
  involved <- which(groups %in% groups[stuck])
  rows_per_group <- split(involved, factor(groups[involved],
                                           levels = unique(groups[stuck])))
  dropped <- list()
  for (rows in rows_per_group) {
    from <- rows[!placeable[rows] & count[rows] > 0]
    to <- rows[placeable[rows] & count[rows] > 0]
    reads <- sum(count[from])
    count[from] <- count[from] * 0L
    if (length(to)) {
      count[to] <- count[to] + as.vector(stats::rmultinom(1L, reads, count[to]))
    } else {
      dropped[[length(dropped) + 1L]] <- data.table::data.table(
        sampling_group = as.character(groups[from[1]]),
        transcript_id = as.character(signal_table$transcript_id[from[1]]),
        reads = as.numeric(reads)
      )
    }
  }
  data.table::set(signal_table, j = column, value = count)
  list(signal_table = signal_table, dropped = data.table::rbindlist(dropped))
}

# Turn signal positions into whole fragments with real coordinates.
#
# Up to here a read is just a position: where the ribosome was. A sequenced
# fragment is a stretch of RNA around that position, and how far it extends in
# each direction depends on its length and on where inside it the ribosome site
# falls. This picks a length and an offset for every read, checks the fragment
# fits inside its transcript, and hands back fragments ready to be written out.
#
# Transcript ends are the awkward case: near them, some length-and-offset
# combinations would run past the edge, and letting a fragment hang over the
# edge would invent sequence. `boundary_action` chooses: "error" stops, and
# "renormalize" keeps only the combinations that do fit and redistributes the
# probability among them, so positions near an end draw from a narrower set of
# geometries. Where none of the configured combinations fits (a ribosome too
# close to the end for the footprints the distribution allows), the position
# gets no reads: they go to other positions of the same region, or are left out
# if no other position of the region has both signal and room.
# Reads left out are listed in the "dropped_reads" attribute of the result.
make_simulated_rpf_fragments <- function(signal_table, transcript_models,
                                         fragment_lengths, fragment_geometry) {
  geometry <- normalize_fragment_geometry(fragment_geometry)
  candidates <- resolve_fragment_distribution(fragment_lengths, geometry)
  if (any(!is.finite(signal_table$score) | signal_table$score < 0) ||
      (!"read_count" %in% names(signal_table) &&
       any(signal_table$score != as.integer(signal_table$score)))) {
    stop("Simulated RPF scores must be non-negative integers")
  }

  candidate_table <- fragment_candidates(
    signal_table, transcript_models, candidates, geometry
  )
  joint <- has_active_fragment_bias(geometry)
  dropped <- NULL
  if (!joint) {
    placeable <- tabulate(candidate_table$signal_row,
                          nbins = nrow(signal_table)) > 0L
    moved <- move_unplaceable_reads(signal_table, placeable)
    signal_table <- moved$signal_table
    dropped <- moved$dropped
  }
  allocated <- allocate_fragment_candidates(candidate_table, signal_table, joint)
  if (joint) dropped <- attr(allocated, "dropped_reads")
  # Only the candidates that were drawn become fragments, often a small part of
  # all candidates. They stay in the order of the candidate table (signal row,
  # then candidate order), which sets the fragment ids.
  drawn <- allocated[score > 0L]
  transcript_id <- signal_table$transcript_id[drawn$signal_row]
  signal_position <- signal_table$signal_position[drawn$signal_row]
  # Each transcript's model is looked up once, the same way as above.
  used_ids <- unique(transcript_id)
  used_models <- lapply(seq_along(used_ids), function(k) {
    transcript_models[[used_ids[k]]]
  })
  model_index <- match(transcript_id, used_ids)
  strand <- vapply(used_models, function(model) model$strand,
                   character(1))[model_index]
  chromosome <- vapply(used_models, function(model) {
    as.character(GenomicRanges::seqnames(model$exons)[1])
  }, character(1))[model_index]
  alignments <- fragment_alignments(used_models, model_index, drawn$tx_start,
                                    drawn$fragment_length)
  start <- alignments$pos
  end <- alignments$end
  # `sequence` is transcript-oriented (5'->3' of the mRNA), matching the
  # biological fragment/read direction; it is what codon/end-motif biases and
  # the ground-truth table use. SAM/BAM SEQ must instead be reported relative
  # to the reference (genome plus-strand), so for minus-strand fragments we
  # store the reverse complement separately rather than reusing `sequence`
  # when writing alignments.
  reference_sequence <- drawn$sequence
  minus <- strand != "+"
  reference_sequence[minus] <- as.character(Biostrings::reverseComplement(
    Biostrings::DNAStringSet(drawn$sequence[minus])
  ))
  result <- data.table::data.table(
    seqnames = chromosome,
    start = start,
    end = end,
    strand = strand,
    cigar = alignments$cigar,
    score = drawn$score,
    sampling_group = drawn$sampling_group,
    transcript_id = transcript_id,
    ribosome_site = signal_position,
    signal_position = signal_position,
    fragment_start = start,
    fragment_end = end,
    five_prime_end = ifelse(minus, end, start),
    three_prime_end = ifelse(minus, start, end),
    fragment_length = drawn$fragment_length,
    site_offset = drawn$site_offset,
    site_reference = geometry$site_reference,
    geometry_probability = drawn$geometry_probability,
    sequence = drawn$sequence,
    reference_sequence = reference_sequence,
    codon = drawn$codon,
    frame = drawn$frame,
    five_prime_kmer = end_bias_kmer(drawn$sequence, geometry$five_prime_bias),
    three_prime_kmer = end_bias_kmer(drawn$sequence, geometry$three_prime_bias, TRUE),
    end_bias_weight = drawn$end_bias_weight,
    codon_bias_weight = drawn$codon_bias_weight,
    frame_bias_weight = drawn$frame_bias_weight,
    fragment_bias_weight = drawn$bias_weight,
    final_probability = drawn$probability
  )
  result[, fragment_id := sprintf("fragment_%d", seq_len(.N))]
  if (is.null(dropped) || !nrow(dropped)) {
    dropped <- data.table::data.table(sampling_group = character(0),
                                      transcript_id = character(0),
                                      reads = numeric(0))
  }
  data.table::setattr(result, "dropped_reads", dropped)
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
    sequence = fragment_table$sequence,
    reference_sequence = fragment_table$reference_sequence,
    codon = fragment_table$codon,
    frame = fragment_table$frame,
    end_bias_weight = fragment_table$end_bias_weight,
    codon_bias_weight = fragment_table$codon_bias_weight,
    frame_bias_weight = fragment_table$frame_bias_weight,
    fragment_bias_weight = fragment_table$fragment_bias_weight
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
    "strand", "cigar", "sequence", "score", "sampling_group",
    "codon", "frame", "five_prime_kmer", "three_prime_kmer",
    "end_bias_weight", "codon_bias_weight", "frame_bias_weight",
    "fragment_bias_weight", "final_probability"
  )
  data.table::fwrite(fragment_table[, ..columns], path, sep = "\t")
  path
}

# Next to the library file: the reads of each transcript region that were left
# out because no position of the region with reads has room for a complete
# fragment (see make_simulated_rpf_fragments()). Written only when something
# was left out.
write_dropped_reads <- function(dropped, file_base) {
  path <- paste0(file_base, "_dropped_reads.tsv")
  if (!nrow(dropped)) {
    # Not left over from an earlier run into the same folder.
    if (file.exists(path)) file.remove(path)
    return(NULL)
  }
  group <- as.character(dropped$sampling_group)
  data.table::fwrite(data.table::data.table(
    transcript_id = dropped$transcript_id,
    region = ifelse(grepl(":", group, fixed = TRUE), sub(":.*", "", group), NA_character_),
    reads = dropped$reads
  ), path, sep = "\t")
  path
}

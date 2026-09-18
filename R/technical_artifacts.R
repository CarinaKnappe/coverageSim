normalize_technical_artifacts <- function(technical_artifacts) {
  defaults <- list(
    duplication_rate = 0,
    duplicate_copies = 1L,
    multimapping_rate = 0,
    secondary_alignments = 1L
  )
  if (is.null(technical_artifacts)) return(defaults)
  if (!is.list(technical_artifacts)) stop("technical_artifacts must be a list")
  unknown <- setdiff(names(technical_artifacts), names(defaults))
  if (length(unknown)) {
    stop("Unknown technical_artifacts field(s): ", paste(unknown, collapse = ", "))
  }
  result <- utils::modifyList(defaults, technical_artifacts)
  for (name in c("duplication_rate", "multimapping_rate")) {
    value <- result[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 0 || value > 1) {
      stop(name, " must be one finite probability from zero to one")
    }
  }
  for (name in c("duplicate_copies", "secondary_alignments")) {
    value <- result[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 1 || value != as.integer(value)) {
      stop(name, " must be one positive integer")
    }
    result[[name]] <- as.integer(value)
  }
  result
}

has_active_technical_artifacts <- function(technical_artifacts) {
  technical_artifacts$duplication_rate > 0 ||
    technical_artifacts$multimapping_rate > 0
}

expand_fragment_records <- function(fragment_table) {
  index <- rep.int(seq_len(nrow(fragment_table)), fragment_table$score)
  records <- data.table::copy(fragment_table[index])
  records[, `:=`(
    qname = sprintf("molecule_%d", seq_len(.N)),
    flag = ifelse(strand == "+", 0L, 16L),
    mapq = 30L, nh = 1L, is_duplicate = FALSE,
    is_secondary = FALSE, source_molecule = seq_len(.N)
  )]
  records
}

simulate_alignment_artifacts <- function(fragment_table, technical_artifacts) {
  settings <- normalize_technical_artifacts(technical_artifacts)
  primary <- expand_fragment_records(fragment_table)
  original_count <- nrow(primary)

  duplicate_count <- stats::rbinom(
    1L, original_count, settings$duplication_rate
  )
  duplicates <- primary[0]
  if (duplicate_count > 0L) {
    selected <- sample.int(original_count, duplicate_count, replace = FALSE)
    duplicate_index <- rep(selected, each = settings$duplicate_copies)
    duplicates <- data.table::copy(primary[duplicate_index])
    duplicates[, `:=`(
      qname = paste0(qname, "_dup", sequence(.N), "_", seq_len(.N)),
      flag = bitwOr(flag, 1024L), is_duplicate = TRUE
    )]
  }

  alignment_key <- paste(primary$seqnames, primary$start, primary$cigar,
                         primary$strand, sep = ":")
  alternatives <- split(seq_len(original_count), primary$fragment_length)
  eligible <- vapply(seq_len(original_count), function(i) {
    any(alignment_key[alternatives[[as.character(primary$fragment_length[i])]]] !=
          alignment_key[i])
  }, logical(1))
  requested_multi <- stats::rbinom(
    1L, original_count, settings$multimapping_rate
  )
  selected_multi <- if (requested_multi > 0L && any(eligible)) {
    sample(which(eligible), min(requested_multi, sum(eligible)), replace = FALSE)
  } else integer()
  secondary <- primary[0]
  if (length(selected_multi)) {
    primary[selected_multi, nh := 1L + settings$secondary_alignments]
    secondary <- data.table::rbindlist(lapply(selected_multi, function(i) {
      choices <- alternatives[[as.character(primary$fragment_length[i])]]
      choices <- choices[alignment_key[choices] != alignment_key[i]]
      chosen <- sample(choices, settings$secondary_alignments, replace =
                         length(choices) < settings$secondary_alignments)
      result <- data.table::copy(primary[chosen])
      result[, `:=`(
        qname = primary$qname[i],
        sequence = primary$sequence[i],
        flag = bitwOr(ifelse(strand == "+", 0L, 16L), 256L),
        nh = 1L + settings$secondary_alignments,
        is_secondary = TRUE,
        source_molecule = primary$source_molecule[i]
      )]
      result
    }))
  }
  records <- data.table::rbindlist(
    list(primary, duplicates, secondary), use.names = TRUE, fill = TRUE
  )
  attr(records, "diagnostics") <- c(
    biological_molecules = original_count,
    pcr_duplicate_records = nrow(duplicates),
    requested_multimapping_molecules = requested_multi,
    multimapping_molecules = length(selected_multi),
    secondary_alignment_records = nrow(secondary),
    output_records = nrow(records)
  )
  records
}

write_artifact_sam <- function(records, path, seqinfo) {
  header <- paste(
    "@SQ", paste0("SN:", GenomeInfoDb::seqnames(seqinfo)),
    paste0("LN:", GenomeInfoDb::seqlengths(seqinfo)), sep = "\t"
  )
  lines <- paste(
    records$qname, records$flag, records$seqnames, records$start,
    records$mapq, records$cigar, "*", 0L, 0L, records$reference_sequence, "*",
    paste0("NH:i:", records$nh), sep = "\t"
  )
  writeLines(c(header, lines), path)
  path
}

write_artifact_library <- function(records, file_base, format, seqinfo) {
  sam_path <- paste0(file_base, ".sam")
  write_artifact_sam(records, sam_path, seqinfo)
  if (format == "sam") return(c(default = sam_path, sam = sam_path))
  bam_path <- Rsamtools::asBam(sam_path, destination = file_base, overwrite = TRUE)
  unlink(sam_path)
  c(default = bam_path, bam = bam_path)
}

write_artifact_ground_truth <- function(records, file_base, ground_truth) {
  if (identical(ground_truth, FALSE) || is.null(ground_truth)) return(NULL)
  path <- if (identical(ground_truth, TRUE)) {
    paste0(file_base, "_artifact_truth.tsv")
  } else {
    dir.create(ground_truth, recursive = TRUE, showWarnings = FALSE)
    file.path(ground_truth, paste0(basename(file_base), "_artifact_truth.tsv"))
  }
  columns <- c(
    "qname", "source_molecule", "fragment_id", "seqnames", "start",
    "strand", "cigar", "fragment_length", "flag", "nh",
    "is_duplicate", "is_secondary"
  )
  data.table::fwrite(records[, ..columns], path, sep = "\t")
  path
}

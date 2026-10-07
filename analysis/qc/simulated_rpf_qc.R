# Biological QC report for CoverageSim simulated RPF output.
#
# Required environment variable:
#   COVSIM_QC_INPUT  Dataset directory containing genome/, reads/, and
#                    simulated_region_counts.tsv.
# Optional:
#   COVSIM_QC_OUT    Output directory (default: <input>/qc).

map_genomic_to_transcript <- function(position, exons, strand_value) {
  exon_rank <- S4Vectors::mcols(exons)$exon_rank
  if (!is.null(exon_rank)) exons <- exons[order(exon_rank)]
  cumulative_start <- cumsum(c(1L, head(GenomicRanges::width(exons), -1L)))
  result <- rep.int(NA_integer_, length(position))
  for (i in seq_along(exons)) {
    hit <- position >= GenomicRanges::start(exons)[i] &
      position <= GenomicRanges::end(exons)[i]
    within <- if (strand_value == "+") {
      position[hit] - GenomicRanges::start(exons)[i]
    } else {
      GenomicRanges::end(exons)[i] - position[hit]
    }
    result[hit] <- cumulative_start[i] + within
  }
  result
}

transcript_layout <- function(txdb_file) {
  txdb <- ORFik::loadTxdb(txdb_file)
  mrna <- ORFik::loadRegion(txdb, "mrna")
  leader <- ORFik::loadRegion(txdb, "leaders")
  cds <- ORFik::loadRegion(txdb, "cds")
  trailer <- ORFik::loadRegion(txdb, "trailers")
  models <- lapply(seq_along(mrna), function(i) {
    exons <- mrna[[i]]
    list(
      exons = exons,
      strand = as.character(unique(GenomicRanges::strand(exons))),
      transcript_length = sum(GenomicRanges::width(exons))
    )
  })
  names(models) <- names(mrna)
  list(
    models = models,
    leader_length = stats::setNames(ORFik::widthPerGroup(leader, FALSE), names(leader)),
    cds_length = stats::setNames(ORFik::widthPerGroup(cds, FALSE), names(cds)),
    trailer_length = stats::setNames(ORFik::widthPerGroup(trailer, FALSE), names(trailer))
  )
}

read_fragment_truth <- function(reads_dir) {
  paths <- sort(list.files(
    reads_dir, pattern = "_ground_truth[.]tsv$", full.names = TRUE
  ))
  if (!length(paths)) stop("No ground-truth TSV files found in: ", reads_dir)
  data.table::rbindlist(lapply(paths, function(path) {
    result <- data.table::fread(path)
    result[, sample := sub("_ground_truth[.]tsv$", "", basename(path))]
    result
  }), use.names = TRUE)
}

annotate_truth_positions <- function(truth, layout) {
  truth <- data.table::copy(truth)
  truth[, c("tx_position", "transcript_length") := {
    model <- layout$models[[transcript_id[1L]]]
    if (is.null(model)) stop("Missing transcript model: ", transcript_id[1L])
    list(
      map_genomic_to_transcript(ribosome_site, model$exons, model$strand),
      model$transcript_length
    )
  }, by = transcript_id]
  if (anyNA(truth$tx_position)) stop("Some ribosome sites could not be mapped")

  truth[, `:=`(
    leader_length = unname(layout$leader_length[transcript_id]),
    cds_length = unname(layout$cds_length[transcript_id]),
    trailer_length = unname(layout$trailer_length[transcript_id])
  )]
  truth[, region := data.table::fcase(
    tx_position <= leader_length, "leader",
    tx_position <= leader_length + cds_length, "cds",
    default = "trailer"
  )]
  truth[, region_position := data.table::fcase(
    region == "leader", tx_position,
    region == "cds", tx_position - leader_length,
    default = tx_position - leader_length - cds_length
  )]
  truth[, region_length := data.table::fcase(
    region == "leader", leader_length,
    region == "cds", cds_length,
    default = trailer_length
  )]
  truth
}

summarize_simulated_rpf_qc <- function(truth, dataset_dir) {
  reads_dir <- file.path(dataset_dir, "reads")
  bam_qc <- truth[, .(truth_reads = sum(score)), by = sample]
  bam_qc[, bam_reads := vapply(sample, function(id) {
    Rsamtools::countBam(file.path(reads_dir, paste0(id, ".bam")))$records
  }, numeric(1))]
  bam_qc[, count_difference := bam_reads - truth_reads]

  list(
    bam_qc = bam_qc,
    length_offset = truth[, .(reads = sum(score)),
                          by = .(sample, fragment_length, site_offset)],
    strand = truth[, .(reads = sum(score)), by = .(sample, strand)],
    splicing = truth[, .(reads = sum(score)),
                     by = .(sample, spliced = grepl("N", cigar, fixed = TRUE))],
    frame = truth[region == "cds", .(reads = sum(score)),
                  by = .(sample, frame = (region_position - 1L) %% 3L)],
    start_profile = truth[
      region == "cds" & region_position >= 1L & region_position <= 121L,
      .(reads = sum(score)),
      by = .(sample, relative_position = region_position - 1L)
    ],
    stop_profile = truth[
      region == "cds" & region_position - cds_length >= -120L &
        region_position - cds_length <= 0L,
      .(reads = sum(score)),
      by = .(sample, relative_position = region_position - cds_length)
    ],
    metagene = truth[, .(reads = sum(score)), by = .(
      sample,
      region,
      bin = pmin(100L, floor((region_position - 1L) / region_length * 100L) + 1L)
    )],
    end_nt = summarize_fragment_end_nt(truth),
    gene_counts = truth[, .(reads = sum(score)), by = .(sample, transcript_id)]
  )
}

write_qc_tables <- function(summaries, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (name in names(summaries)) {
    data.table::fwrite(
      summaries[[name]], file.path(output_dir, paste0(name, ".tsv")), sep = "\t"
    )
  }
}

qc_theme <- function() {
  ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

condition_from_sample <- function(sample) sub("^RFP_([^_]+)_.*$", "\\1", sample)

summarize_fragment_end_nt <- function(truth) {
  five_prime <- truth[, .(reads = sum(score)),
                      by = .(sample, nucleotide = substr(sequence, 1L, 1L))]
  five_prime[, end := "5-prime"]
  three_prime <- truth[, .(reads = sum(score)), by = .(
    sample,
    nucleotide = substr(sequence, nchar(sequence), nchar(sequence))
  )]
  three_prime[, end := "3-prime"]
  data.table::rbindlist(list(five_prime, three_prime), use.names = TRUE)
}

plot_qc_report <- function(summaries, region_counts, output_file) {
  grDevices::pdf(output_file, width = 10, height = 7, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  counts <- summaries$bam_qc
  print(ggplot2::ggplot(counts, ggplot2::aes(sample, bam_reads, fill = sample)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::labs(title = "Reads per library", x = NULL, y = "BAM reads") + qc_theme())

  regions <- data.table::melt(region_counts, id.vars = "region",
                              variable.name = "sample", value.name = "reads")
  regions[, fraction := reads / sum(reads), by = sample]
  print(ggplot2::ggplot(regions, ggplot2::aes(sample, fraction, fill = region)) +
    ggplot2::geom_col() +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::labs(title = "Requested signal by transcript region", x = NULL, y = NULL) +
    qc_theme())

  print(ggplot2::ggplot(summaries$length_offset,
                        ggplot2::aes(factor(fragment_length), reads, fill = sample)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::labs(title = "Simulated RPF length distribution",
                  x = "Fragment length (nt)", y = "Reads") + qc_theme())

  print(ggplot2::ggplot(summaries$length_offset,
                        ggplot2::aes(factor(fragment_length), site_offset,
                                     size = reads, colour = sample)) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::scale_size_area() +
    ggplot2::labs(title = "A-site offset from biological 5-prime end",
                  x = "Fragment length (nt)", y = "A-site offset (nt)") + qc_theme())

  frame <- summaries$frame
  frame[, fraction := reads / sum(reads), by = sample]
  print(ggplot2::ggplot(frame, ggplot2::aes(factor(frame), fraction, fill = sample)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::labs(title = "CDS A-site reading-frame periodicity",
                  x = "CDS frame", y = "Fraction of CDS reads") + qc_theme())

  print(ggplot2::ggplot(summaries$start_profile,
                        ggplot2::aes(relative_position, reads, colour = sample)) +
    ggplot2::geom_line() +
    ggplot2::labs(title = "A-site coverage downstream of CDS start",
                  x = "Position relative to CDS start (nt)", y = "Reads") + qc_theme())

  print(ggplot2::ggplot(summaries$stop_profile,
                        ggplot2::aes(relative_position, reads, colour = sample)) +
    ggplot2::geom_line() +
    ggplot2::labs(title = "A-site coverage upstream of CDS end",
                  x = "Position relative to final CDS nucleotide", y = "Reads") + qc_theme())

  print(ggplot2::ggplot(summaries$metagene,
                        ggplot2::aes(bin, reads, colour = sample)) +
    ggplot2::geom_line() + ggplot2::facet_wrap(~region, scales = "free_y") +
    ggplot2::labs(title = "Scaled transcript-region metagene profiles",
                  x = "Relative region bin (1-100)", y = "Reads") + qc_theme())

  end_nt <- summaries$end_nt
  end_nt[, fraction := reads / sum(reads), by = .(sample, end)]
  print(ggplot2::ggplot(end_nt,
                        ggplot2::aes(nucleotide, fraction, fill = sample)) +
    ggplot2::geom_col(position = "dodge") + ggplot2::facet_wrap(~end) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::labs(title = "Fragment-end nucleotide composition",
                  x = "Nucleotide", y = "Fraction") + qc_theme())

  splicing <- summaries$splicing
  splicing[, fraction := reads / sum(reads), by = sample]
  print(ggplot2::ggplot(splicing,
                        ggplot2::aes(sample, fraction, fill = spliced)) +
    ggplot2::geom_col() + ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::labs(title = "Reads crossing an exon junction (CIGAR N)",
                  x = NULL, y = "Fraction") + qc_theme())

  genes <- data.table::copy(summaries$gene_counts)
  genes[, condition := condition_from_sample(sample)]
  genes <- genes[, .(reads = mean(reads)), by = .(condition, transcript_id)]
  genes <- data.table::dcast(genes, transcript_id ~ condition, value.var = "reads", fill = 0)
  if (all(c("WT", "KO") %in% names(genes))) {
    print(ggplot2::ggplot(genes, ggplot2::aes(log10(WT + 1), log10(KO + 1))) +
      ggplot2::geom_point(alpha = 0.25, size = 0.8) +
      ggplot2::geom_abline(slope = 1, intercept = 0, colour = "red", linetype = 2) +
      ggplot2::labs(title = "Mean transcript-level signal: WT versus KO",
                    x = "log10(mean WT reads + 1)", y = "log10(mean KO reads + 1)") +
      qc_theme())
  }
}

run_simulated_rpf_qc <- function(dataset_dir, output_dir = file.path(dataset_dir, "qc")) {
  required <- c("data.table", "ggplot2", "scales", "Rsamtools", "ORFik")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
  dataset_dir <- normalizePath(dataset_dir, mustWork = TRUE)
  txdb_file <- file.path(dataset_dir, "genome", "human_flavoured_sim.gtf.db")
  region_file <- file.path(dataset_dir, "simulated_region_counts.tsv")
  if (!file.exists(txdb_file) || !file.exists(region_file)) {
    stop("Dataset lacks TxDb or simulated_region_counts.tsv")
  }
  message("Reading fragment ground truth...")
  truth <- read_fragment_truth(file.path(dataset_dir, "reads"))
  stopifnot(all(truth$site_reference == "a_site"))
  message("Mapping A-sites to transcript coordinates...")
  truth <- annotate_truth_positions(truth, transcript_layout(txdb_file))
  message("Summarizing QC metrics...")
  summaries <- summarize_simulated_rpf_qc(truth, dataset_dir)
  region_counts <- data.table::fread(region_file)
  write_qc_tables(summaries, output_dir)
  plot_qc_report(
    summaries, region_counts, file.path(output_dir, "simulated_rpf_qc.pdf")
  )
  invisible(list(summaries = summaries, output_dir = output_dir))
}

if (identical(Sys.getenv("COVSIM_QC_AUTORUN", unset = "true"), "true")) {
  input_dir <- Sys.getenv("COVSIM_QC_INPUT", unset = "")
  if (!nzchar(input_dir)) stop("Set COVSIM_QC_INPUT to the simulation dataset directory")
  output_dir <- Sys.getenv("COVSIM_QC_OUT", unset = file.path(input_dir, "qc"))
  run_simulated_rpf_qc(input_dir, output_dir)
}

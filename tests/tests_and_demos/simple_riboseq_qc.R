# Simple, browser-like coverage plots for simulated CoverageSim RPFs.
#
# Environment variables:
#   COVSIM_QC_INPUT   simulation dataset directory (required)
#   COVSIM_QC_OUT     output directory (default: <input>/qc)
#   COVSIM_QC_SAMPLE  sample to display (default: RFP_WT_1)

load_simulated_rpf_qc_helpers <- function() {
  old <- Sys.getenv("COVSIM_QC_AUTORUN", unset = NA_character_)
  on.exit({
    if (is.na(old)) Sys.unsetenv("COVSIM_QC_AUTORUN") else
      Sys.setenv(COVSIM_QC_AUTORUN = old)
  })
  Sys.setenv(COVSIM_QC_AUTORUN = "false")
  candidates <- c(
    file.path("tests", "tests_and_demos", "simulated_rpf_qc.R"),
    file.path("..", "tests_and_demos", "simulated_rpf_qc.R")
  )
  helper <- candidates[file.exists(candidates)][1L]
  if (is.na(helper)) stop("Cannot locate simulated_rpf_qc.R")
  sys.source(helper, envir = environment(load_simulated_rpf_qc_helpers))
}

fragment_coverage_vectors <- function(truth, transcript_length) {
  starts <- truth$tx_position - truth$site_offset
  ends <- starts + truth$fragment_length - 1L
  if (any(starts < 1L | ends > transcript_length)) {
    stop("Fragment geometry exceeds transcript bounds")
  }
  delta <- numeric(transcript_length + 1L)
  additions <- data.table::data.table(position = starts, value = truth$score)[,
    .(value = sum(value)), by = position
  ]
  removals <- data.table::data.table(position = ends + 1L, value = truth$score)[,
    .(value = sum(value)), by = position
  ]
  delta[additions$position] <- delta[additions$position] + additions$value
  delta[removals$position] <- delta[removals$position] - removals$value
  fragment_coverage <- cumsum(delta)[seq_len(transcript_length)]

  a_site <- numeric(transcript_length)
  sites <- truth[, .(value = sum(score)), by = tx_position]
  a_site[sites$tx_position] <- sites$value
  list(fragment = fragment_coverage, a_site = a_site)
}

select_example_transcripts <- function(truth, n_per_strand = 1L) {
  totals <- truth[, .(
    reads = sum(score),
    spliced_reads = sum(score[grepl("N", cigar, fixed = TRUE)])
  ), by = .(transcript_id, strand)]
  totals[spliced_reads > 0L][order(-reads), head(.SD, n_per_strand), by = strand]
}

uorf_transcript_ranges <- function(uorfs, transcript_id, model) {
  selected <- uorfs[names(uorfs) == transcript_id]
  if (!length(selected)) return(data.table::data.table(start = integer(), end = integer()))
  ranges <- unlist(selected, use.names = FALSE)
  first <- map_genomic_to_transcript(
    GenomicRanges::start(ranges), model$exons, model$strand
  )
  last <- map_genomic_to_transcript(
    GenomicRanges::end(ranges), model$exons, model$strand
  )
  data.table::data.table(start = pmin(first, last), end = pmax(first, last))
}

draw_annotation <- function(transcript_length, leader_length, cds_length,
                            uorf_ranges, transcript_id, strand_value) {
  graphics::plot(NA, xlim = c(1, transcript_length), ylim = c(0, 1),
                 axes = FALSE, xlab = "", ylab = "",
                 main = paste0(transcript_id, "  strand ", strand_value))
  cds_start <- leader_length + 1L
  cds_end <- leader_length + cds_length
  graphics::rect(1, 0.25, leader_length, 0.55, col = "grey75", border = NA)
  graphics::rect(cds_start, 0.25, cds_end, 0.55, col = "tomato", border = NA)
  graphics::rect(cds_end + 1L, 0.25, transcript_length, 0.55,
                 col = "goldenrod", border = NA)
  if (nrow(uorf_ranges)) {
    graphics::rect(uorf_ranges$start, 0.65, uorf_ranges$end, 0.9,
                   col = "darkseagreen4", border = NA)
  }
  graphics::axis(1)
  graphics::legend("topright", legend = c("leader", "CDS", "trailer", "uORF"),
                   fill = c("grey75", "tomato", "goldenrod", "darkseagreen4"),
                   horiz = TRUE, bty = "n", cex = 0.8)
}

draw_coverage_bars <- function(values, title, colour, xlim = NULL,
                               cds_start = NULL, xlab = "") {
  if (is.null(xlim)) xlim <- c(1L, length(values))
  positions <- seq.int(max(1L, floor(xlim[1])), min(length(values), ceiling(xlim[2])))
  ymax <- max(values[positions], 1)
  graphics::plot(positions, values[positions], type = "h", lend = 1,
                 col = colour, lwd = 1, xlim = xlim, ylim = c(0, ymax * 1.08),
                 xlab = xlab, ylab = "Read count", main = title)
  if (!is.null(cds_start)) graphics::abline(v = cds_start, col = "red", lwd = 2)
}

plot_transcript_page <- function(truth, layout, uorfs, tx_id, sample) {
  model <- layout$models[[tx_id]]
  selected <- truth[transcript_id == tx_id]
  selected[, tx_position := map_genomic_to_transcript(
    ribosome_site, model$exons, model$strand
  )]
  coverage <- fragment_coverage_vectors(selected, model$transcript_length)
  leader_length <- unname(layout$leader_length[tx_id])
  cds_length <- unname(layout$cds_length[tx_id])
  cds_start <- leader_length + 1L
  uorf_ranges <- uorf_transcript_ranges(uorfs, tx_id, model)

  graphics::par(mfrow = c(4, 1), mar = c(3.2, 4.2, 2.4, 1), oma = c(0, 0, 2, 0))
  draw_annotation(model$transcript_length, leader_length, cds_length,
                  uorf_ranges, tx_id, model$strand)
  draw_coverage_bars(
    coverage$fragment,
    "Simulated RPF coverage: how many RPFs cover each nucleotide",
    "steelblue4", cds_start = cds_start
  )
  draw_coverage_bars(
    coverage$a_site,
    "A-site coverage: one ribosome position per RPF",
    "firebrick3", cds_start = cds_start
  )
  draw_coverage_bars(
    coverage$a_site,
    "A-site coverage around the CDS start (red line)",
    "firebrick3", xlim = c(max(1L, cds_start - 60L), cds_start + 150L),
    cds_start = cds_start, xlab = "Transcript position (nt)"
  )
  graphics::mtext(paste(sample, "- direct transcript-level Ribo-seq view"),
                  outer = TRUE, cex = 1.1, font = 2)
}

run_simple_riboseq_qc <- function(dataset_dir, output_dir, sample = "RFP_WT_1") {
  load_simulated_rpf_qc_helpers()
  dataset_dir <- normalizePath(dataset_dir, mustWork = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  truth_file <- file.path(dataset_dir, "reads", paste0(sample, "_ground_truth.tsv"))
  truth <- data.table::fread(truth_file)
  layout <- transcript_layout(file.path(
    dataset_dir, "genome", "human_flavoured_sim.gtf.db"
  ))
  uorfs <- readRDS(file.path(dataset_dir, "genome", "true_uORFs.rds"))
  examples <- select_example_transcripts(truth)
  data.table::fwrite(examples, file.path(output_dir, "simple_qc_transcripts.tsv"), sep = "\t")

  pdf_file <- file.path(output_dir, "simple_riboseq_coverage.pdf")
  grDevices::pdf(pdf_file, width = 12, height = 10, onefile = TRUE)
  for (transcript_id in examples$transcript_id) {
    plot_transcript_page(truth, layout, uorfs, transcript_id, sample)
  }
  grDevices::dev.off()

  for (transcript_id in examples$transcript_id) {
    png_file <- file.path(output_dir, paste0("coverage_", transcript_id, ".png"))
    grDevices::png(png_file, width = 1500, height = 1400, res = 130)
    plot_transcript_page(truth, layout, uorfs, transcript_id, sample)
    grDevices::dev.off()
  }
  invisible(list(pdf = pdf_file, examples = examples))
}

if (identical(Sys.getenv("COVSIM_SIMPLE_QC_AUTORUN", unset = "true"), "true")) {
  dataset_dir <- Sys.getenv("COVSIM_QC_INPUT", unset = "")
  if (!nzchar(dataset_dir)) stop("Set COVSIM_QC_INPUT")
  output_dir <- Sys.getenv("COVSIM_QC_OUT", unset = file.path(dataset_dir, "qc"))
  sample <- Sys.getenv("COVSIM_QC_SAMPLE", unset = "RFP_WT_1")
  run_simple_riboseq_qc(dataset_dir, output_dir, sample)
}

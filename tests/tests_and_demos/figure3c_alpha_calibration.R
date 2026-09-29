repo_dir <- normalizePath(
  Sys.getenv("COVSIM_REPO", unset = getwd()), mustWork = TRUE
)
devtools::load_all(repo_dir)
data.table::setDTthreads(4)

source_run <- normalizePath(Sys.getenv(
  "COVSIM_HUMAN_LEARNING_RUN",
  unset = file.path(repo_dir, "rust_helpers", "runs", "2026-09-11_human-learning")
), mustWork = TRUE)
out_dir <- Sys.getenv(
  "COVSIM_FIGURE3C_OUTPUT",
  unset = file.path(tempdir(), "coverageSim_figure3c_alpha_calibration")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Record which code produced this run. Without it, a stored result cannot be
# traced back to the source that made it -- and because load_all() above reads
# the working directory, unsaved edits count too, which is what this records.
write_code_version(out_dir)

alpha_scales <- as.numeric(strsplit(Sys.getenv(
  "COVSIM_FIGURE3C_ALPHA_SCALES",
  unset = "1,0.01,0.003,0.001,0.0005,0.00025,0.0001,0.00005,0.00002"
), ",", fixed = TRUE)[[1]])
if (!length(alpha_scales) || any(!is.finite(alpha_scales) | alpha_scales <= 0)) {
  stop("COVSIM_FIGURE3C_ALPHA_SCALES must contain positive finite numbers")
}
simulation_replicates <- as.integer(Sys.getenv(
  "COVSIM_FIGURE3C_REPLICATES", unset = "5"
))
if (length(simulation_replicates) != 1L || is.na(simulation_replicates) ||
    simulation_replicates < 1L) {
  stop("COVSIM_FIGURE3C_REPLICATES must be one positive integer")
}

prepared <- readRDS(file.path(source_run, "prepared.rds"))
mapped <- data.table::as.data.table(readRDS(file.path(source_run, "mapped.rds")))
geometry <- data.table::fread(file.path(source_run, "learned_geometry.tsv"))
fit <- readRDS(file.path(source_run, "codon_only_fit.rds"))

cds_bounds <- data.table::rbindlist(lapply(names(prepared$models), function(id) {
  positions <- coverageSim:::learning_cds_positions(
    prepared$models[[id]], prepared$cds[[id]]
  )
  data.table::data.table(
    transcript_id = id,
    cds_start = positions[1L],
    cds_end = positions[length(positions)] + 2L,
    cds_width = length(positions) * 3L
  )
}))

mapped <- mapped[transcript_id %in% prepared$test_ids]
mapped[, site_offset := geometry$site_offset[match(fragment_length, geometry$fragment_length)]]
mapped <- mapped[!is.na(site_offset)]
mapped[, site_tx := tx_start + site_offset]
mapped <- merge(mapped, cds_bounds, by = "transcript_id")
mapped <- mapped[site_tx >= cds_start & site_tx <= cds_end]
real_coverage <- mapped[, .(score = sum(count)), by = .(transcript_id, site_tx)]
real_counts <- real_coverage[, .(count = sum(score)), by = transcript_id]
ids <- prepared$test_ids[prepared$test_ids %in% real_counts$transcript_id]
if (length(ids) < 50L) stop("Too few held-out transcripts have supported real reads")

counts <- real_counts$count[match(ids, real_counts$transcript_id)]

codon_weights <- data.table::copy(fit$diagnostics$codon_weights)
seq_bias <- codon_weights[, .(
  variable = "heldout_codon_fit", seqs = codon, alpha = weight * 1000
)]
seq_bias <- data.table::rbindlist(list(
  seq_bias,
  data.table::data.table(
    variable = "heldout_codon_fit",
    seqs = c("###", "%%%", "&&&"),
    alpha = mean(seq_bias$alpha)
  )
))

relative_peaks <- function(coverage, label, scale = NA_real_, replicate = NA_integer_) {
  coverage <- merge(coverage, cds_bounds, by = "transcript_id")
  coverage <- coverage[transcript_id %in% ids & site_tx >= cds_start & site_tx <= cds_end]
  peaks <- coverage[order(site_tx), .SD[which.max(score)], by = transcript_id]
  peaks[, relative_position := 100 * (site_tx - cds_start) / pmax(1, cds_width - 1L)]
  peaks[, `:=`(profile = label, dmn_alpha_scale = scale,
               simulation_replicate = replicate)]
  peaks[, .(profile, dmn_alpha_scale, simulation_replicate,
            transcript_id, site_tx, score,
            cds_start, cds_end, cds_width, relative_position)]
}

# Build the same final Step-4 alpha rows used by simNGScoverage.  Figure 3C
# only needs A-site coverage, so fragment construction and BAM writing would
# add runtime without changing this statistic.
codons <- coverageSim:::translate_orf_seq(
  prepared$cds[ids], prepared$fasta, is.sorted = TRUE, as = "codon",
  start.as.hash = TRUE, startp1.as.per = TRUE, stopm1.as.amp = TRUE,
  return.as.list = TRUE
)
alpha_lookup <- stats::setNames(seq_bias$alpha, seq_bias$seqs)
codon_alpha <- lapply(codons, function(x) unname(alpha_lookup[x]))
if (anyNA(unlist(codon_alpha, use.names = FALSE))) {
  stop("Learned codon profile does not cover every held-out CDS codon")
}
set.seed(20260918)
base_alpha <- coverageSim:::sim_sequence_bias(
  quote(rep(1, x)),
  lengths = cds_bounds[match(ids, transcript_id), cds_width],
  alpha_matrix = codon_alpha,
  seq_acf = shapes(9),
  rnase_acf = c(0.5, 2, 1, 10, 2, 1, 0.5)
)

sample_step4 <- function(scale, seed) {
  alpha <- coverageSim:::scale_dmn_alpha(base_alpha, scale)
  row_lengths <- lengths(alpha)
  alpha_matrix <- coverageSim:::pack_alpha_rows(
    alpha,
    coverageSim:::list_to_mat(
      cds_bounds[match(ids, transcript_id), cds_width], rnase_length = 6L
    )
  )
  set.seed(seed)
  sampled <- extraDistr::rdirmnom(
    n = length(ids), size = counts, alpha = alpha_matrix
  )
  sampled <- split(
    coverageSim:::flatten_sample_rows(sampled, row_lengths),
    rep(ids, row_lengths)
  )
  coverage <- data.table::rbindlist(lapply(ids, function(id) {
    values <- sampled[[id]]
    width <- cds_bounds[transcript_id == id, cds_width]
    # RNase convolution adds three positions at each side. Figure 3C is CDS-only.
    values <- values[seq.int(4L, width + 3L)]
    data.table::data.table(
      transcript_id = id,
      site_tx = cds_bounds[transcript_id == id, cds_start] + seq_len(width) - 1L,
      score = values
    )
  }))
  list(coverage = coverage, total_reads = sum(unlist(sampled, use.names = FALSE)))
}

peak_tables <- list(relative_peaks(real_coverage, "Real human reads"))
simulation_totals <- list()
for (scale in alpha_scales) {
  label <- paste0("Simulated (alpha scale ", format(scale, scientific = TRUE), ")")
  for (replicate in seq_len(simulation_replicates)) {
    simulation <- sample_step4(scale, seed = 20260918L + replicate - 1L)
    simulated <- simulation$coverage
    observed_total <- simulation$total_reads
    expected_total <- sum(counts)
    if (observed_total != expected_total) {
      stop("Read budget changed for alpha scale ", scale,
           ", replicate ", replicate)
    }
    key <- paste(label, replicate, sep = "_")
    peak_tables[[key]] <- relative_peaks(
      simulated, label, scale, replicate = replicate
    )
    simulation_totals[[key]] <- data.table::data.table(
      profile = label, dmn_alpha_scale = scale,
      simulation_replicate = replicate,
      expected_reads = expected_total, observed_reads = observed_total
    )
  }
}

peaks <- data.table::rbindlist(peak_tables, use.names = TRUE)
bin_distribution <- function(values, bins = 50L) {
  index <- pmin(bins, pmax(1L, floor(values / 100 * bins) + 1L))
  tabulate(index, nbins = bins) / length(values)
}
real_distribution <- bin_distribution(
  peaks[profile == "Real human reads", relative_position]
)
summary <- data.table::rbindlist(lapply(setdiff(unique(peaks$profile), "Real human reads"), function(label) {
  simulated <- peaks[profile == label]
  distribution <- bin_distribution(simulated$relative_position)
  probabilities <- seq(0, 1, length.out = 101L)
  simulated_quantiles <- stats::quantile(
    simulated$relative_position, probabilities, names = FALSE
  )
  real_quantiles <- stats::quantile(
    peaks[profile == "Real human reads", relative_position],
    probabilities, names = FALSE
  )
  data.table::data.table(
    profile = label,
    dmn_alpha_scale = unique(simulated$dmn_alpha_scale),
    transcripts = data.table::uniqueN(simulated$transcript_id),
    simulation_replicates = data.table::uniqueN(simulated$simulation_replicate),
    total_variation_to_real = 0.5 * sum(abs(distribution - real_distribution)),
    mean_absolute_quantile_distance = mean(abs(simulated_quantiles - real_quantiles)),
    ks_distance = suppressWarnings(stats::ks.test(
      simulated$relative_position,
      peaks[profile == "Real human reads", relative_position], exact = FALSE
    )$statistic),
    first_10_percent = mean(simulated$relative_position <= 10),
    median_peak_position = stats::median(simulated$relative_position)
  )
}))
real_summary <- data.table::data.table(
  profile = "Real human reads", dmn_alpha_scale = NA_real_,
  transcripts = sum(peaks$profile == "Real human reads"),
  simulation_replicates = NA_integer_,
  total_variation_to_real = 0,
  mean_absolute_quantile_distance = 0,
  ks_distance = 0,
  first_10_percent = mean(
    peaks[profile == "Real human reads", relative_position] <= 10
  ),
  median_peak_position = stats::median(
    peaks[profile == "Real human reads", relative_position]
  )
)
summary <- data.table::rbindlist(list(real_summary, summary))
best <- summary[!is.na(dmn_alpha_scale)][which.min(mean_absolute_quantile_distance)]
summary[, selected := profile == best$profile]

data.table::fwrite(peaks, file.path(out_dir, "figure3c_peak_positions.tsv"), sep = "\t")
data.table::fwrite(summary, file.path(out_dir, "figure3c_alpha_summary.tsv"), sep = "\t")
data.table::fwrite(
  data.table::rbindlist(simulation_totals),
  file.path(out_dir, "simulation_read_budgets.tsv"), sep = "\t"
)

plot_profiles <- function(plot_data, title, binwidth) {
  plot_data[, profile := factor(profile, levels = unique(profile))]
  ggplot2::ggplot(plot_data, ggplot2::aes(x = relative_position, group = profile)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density * binwidth)),
      binwidth = binwidth, boundary = 0, fill = "grey45", color = "grey45"
    ) +
    ggplot2::facet_grid(profile ~ ., scales = "free_y") +
    ggplot2::scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
    ggplot2::labs(
      title = title,
      subtitle = paste(
        length(ids), "held-out, well-covered human CDSs;",
        simulation_replicates, "simulation draws"
      ),
      x = "Position in CDS [%]", y = "Fraction of maximum peaks"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(strip.background = ggplot2::element_blank(),
                   strip.text.y = ggplot2::element_text(angle = 0))
}

primary <- peaks[profile %in% c("Real human reads", best$profile)]
primary[, profile := factor(profile, levels = c("Real human reads", best$profile))]
primary_plot <- plot_profiles(
  primary, "Figure 3C-style maximum-peak positions", binwidth = 1
)
diagnostic_order <- c("Real human reads", summary[!is.na(dmn_alpha_scale)][
  order(-dmn_alpha_scale), profile])
diagnostic <- peaks[profile %in% diagnostic_order]
diagnostic[, profile := factor(profile, levels = diagnostic_order)]
diagnostic_plot <- plot_profiles(
  diagnostic, "DMN alpha-scale sensitivity", binwidth = 2
)

ggplot2::ggsave(file.path(out_dir, "figure3c_real_vs_selected.png"), primary_plot,
                width = 8, height = 5.5, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "figure3c_real_vs_selected.pdf"), primary_plot,
                width = 8, height = 5.5)
ggplot2::ggsave(file.path(out_dir, "figure3c_all_alpha_scales.png"), diagnostic_plot,
                width = 8, height = 11, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "figure3c_all_alpha_scales.pdf"), diagnostic_plot,
                width = 8, height = 11)

data.table::fwrite(data.table::data.table(
  field = c("created", "source_run", "real_bam", "transcripts",
            "simulation_replicates", "selected_scale", "selection_metric", "seed"),
  value = c(as.character(Sys.time()), source_run, prepared$bam, length(ids),
            simulation_replicates, best$dmn_alpha_scale,
            "minimum mean absolute empirical-quantile distance to real peaks",
            "20260918 + replicate - 1")
), file.path(out_dir, "manifest.tsv"), sep = "\t")

print(summary)
message("Selected dmn_alpha_scale: ", best$dmn_alpha_scale)
message("Outputs: ", normalizePath(out_dir))

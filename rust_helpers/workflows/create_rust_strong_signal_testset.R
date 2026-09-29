## Create small, deterministic RUST input datasets with a deliberately strong
## frame-0 signal. The generated data are local analysis artefacts and are
## written below the local coverageSim_data/runs/ directory.

repo_dir <- normalizePath(
  Sys.getenv("COVSIM_REPO", unset = getwd()), mustWork = TRUE
)
devtools::load_all(repo_dir)
library(coverageSim)
library(ORFik)
library(SummarizedExperiment)

set.seed(42)

output_root <- normalizePath(
  Sys.getenv(
    "COVSIM_RUST_TESTSET",
    unset = file.path(
      Sys.getenv("COVSIM_DATA_ROOT", unset = file.path(dirname(repo_dir), "coverageSim_data")),
      "runs", "2026-09-09_rust-frame-signal"
    )
  ),
  mustWork = FALSE
)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

# Record which code produced this run. Without it, a stored result cannot be
# traced back to the source that made it -- and because load_all() above reads
# the working directory, unsaved edits count too, which is what this records.
write_code_version(output_root, repo_dir = repo_dir)

sim_genome <- simGenome(
  n = 1L,
  out_dir = file.path(output_root, "genome"),
  genome_name = "rust_strong_signal",
  leader_length = 60L,
  cds_length = 300L,
  trailer_length = 60L,
  cds_exons = 2L,
  cds_intron_length = 30L,
  max_uorfs = 0L,
  debug_on = FALSE
)
cds <- ORFik::loadRegion(sim_genome["txdb"], "cds")

# Use a fixed count instead of a negative-binomial draw. This makes every
# scenario directly comparable and leaves enough reads for a clear frame peak.
count_table <- SummarizedExperiment::SummarizedExperiment(
  assays = list(
    gene = matrix(
      100000L, nrow = length(cds), ncol = 1L,
      dimnames = list(names(cds), "RFP_strong")
    )
  ),
  rowRanges = cds,
  colData = S4Vectors::DataFrame(
    libtype = factor("RFP"),
    condition = factor("strong_signal"),
    replicate = "1",
    row.names = "RFP_strong"
  )
)
count_table <- simCountTablesRegions(
  count_table,
  regionsToSample = "cds",
  region_proportion = list(cds = list(RFP = 1)),
  sampling = c(RFP = "MN")
)

# Keep signal away from both CDS boundaries so all requested footprints are
# complete. The expression is evaluated with x = CDS width by coverageSim.
strong_frame0 <- quote({
  result <- rep(0, x)
  usable <- 31:(x - 30)
  result[usable[((usable - 1) %% 3) == 0]] <- 1
  result
})

scenarios <- list(
  single_length_28 = data.frame(
    fragment_length = 28L, site_offset = 15L, probability = 1
  ),
  mixed_lengths_correct_offsets = data.frame(
    fragment_length = c(27L, 28L, 29L),
    site_offset = c(14L, 15L, 15L),
    probability = c(1 / 3, 1 / 3, 1 / 3)
  ),
  mixed_lengths_wrong_common_offset = data.frame(
    fragment_length = c(27L, 28L, 29L),
    site_offset = c(13L, 15L, 17L),
    probability = c(1 / 3, 1 / 3, 1 / 3)
  )
)

cds_start <- min(GenomicRanges::start(cds[[1L]]))
summary_rows <- list()

for (scenario_name in names(scenarios)) {
  scenario_dir <- file.path(output_root, scenario_name)
  if (dir.exists(scenario_dir)) unlink(scenario_dir, recursive = TRUE)
  reads_dir <- file.path(scenario_dir, "reads")
  exp_dir <- file.path(scenario_dir, "experiment")
  dir.create(reads_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(exp_dir, recursive = TRUE, showWarnings = FALSE)

  geometry <- list(
    source = "user",
    site_reference = "a_site",
    distribution = scenarios[[scenario_name]],
    boundary_action = "error",
    five_prime_bias = list(source = "none"),
    three_prime_bias = list(source = "none")
  )

  sim_exp <- simNGScoverage(
    simGenome = sim_genome,
    count_table = count_table,
    out_dir = reads_dir,
    exp_name = paste0("rust_", scenario_name),
    exp_save_dir = exp_dir,
    ideal_coverage = list(cds = list(RFP = strong_frame0)),
    rnase_bias = list(RFP = 1),
    auto_correlation = list(cds = list(RFP = 1)),
    read_lengths_per = list(RFP = scenarios[[scenario_name]]$fragment_length),
    sampling = list(cds = list(RFP = "MN")),
    seq_bias = load_seq_bias(bias = "R2"),
    fragment_geometry = geometry,
    ground_truth = TRUE,
    libFormats = list(RFP = "bam"),
    validate = TRUE
  )

  bam <- ORFik::filepath(sim_exp, "default")[[1L]]
  truth <- data.table::fread(
    file.path(reads_dir, "RFP_strong_ground_truth.tsv")
  )
  truth[, frame := (signal_position - cds_start) %% 3L]
  # This is what a codon-wise Rust pass using one common offset would see.
  # The length-specific column uses the simulated geometry offset and should
  # remain concentrated in frame 0.
  truth[, rust_frame_common_15 := (five_prime_end + 15L - cds_start) %% 3L]
  truth[, rust_frame_length_specific :=
          (five_prime_end + site_offset - cds_start) %% 3L]
  frame_counts <- truth[, .(reads = sum(score)), by = frame][order(frame)]
  common_counts <- truth[, .(reads = sum(score)), by = rust_frame_common_15][
    order(rust_frame_common_15)
  ]
  geometry_counts <- truth[, .(reads = sum(score)), by = rust_frame_length_specific][
    order(rust_frame_length_specific)
  ]
  lengths <- truth[, .(reads = sum(score)), by = fragment_length][order(fragment_length)]
  data.table::fwrite(frame_counts, file.path(scenario_dir, "frame_counts.tsv"), sep = "\t")
  data.table::fwrite(
    common_counts, file.path(scenario_dir, "rust_common_offset_frame_counts.tsv"),
    sep = "\t"
  )
  data.table::fwrite(
    geometry_counts, file.path(scenario_dir, "rust_length_specific_frame_counts.tsv"),
    sep = "\t"
  )
  data.table::fwrite(lengths, file.path(scenario_dir, "length_counts.tsv"), sep = "\t")
  frame_count <- function(table, column, value) {
    result <- table[get(column) == value, reads]
    if (length(result)) result[[1L]] else 0L
  }
  data.table::fwrite(
    data.table::data.table(
      scenario = scenario_name,
      bam = bam,
      total_reads = sum(truth$score),
      true_frame0_reads = frame_count(frame_counts, "frame", 0L),
      common_offset_15_frame0 = frame_count(common_counts, "rust_frame_common_15", 0L),
      common_offset_15_frame1 = frame_count(common_counts, "rust_frame_common_15", 1L),
      common_offset_15_frame2 = frame_count(common_counts, "rust_frame_common_15", 2L),
      length_specific_frame0 = frame_count(
        geometry_counts, "rust_frame_length_specific", 0L
      )
    ),
    file.path(scenario_dir, "summary.tsv"), sep = "\t"
  )
  summary_rows[[scenario_name]] <- data.table::data.table(
    scenario = scenario_name,
    bam = bam,
    total_reads = sum(truth$score),
    true_frame0_reads = frame_count(frame_counts, "frame", 0L),
    common_offset_15_frame0 = frame_count(common_counts, "rust_frame_common_15", 0L),
    common_offset_15_frame1 = frame_count(common_counts, "rust_frame_common_15", 1L),
    common_offset_15_frame2 = frame_count(common_counts, "rust_frame_common_15", 2L),
    length_specific_frame0 = frame_count(
      geometry_counts, "rust_frame_length_specific", 0L
    )
  )
}

data.table::fwrite(
  data.table::rbindlist(summary_rows),
  file.path(output_root, "summary.tsv"), sep = "\t"
)
message("Wrote RUST strong-signal testset to: ", output_root)

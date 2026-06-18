rm(list = ls(all.names = TRUE))
gc(reset = TRUE)

repo_dir <- normalizePath(
  Sys.getenv("COVSIM_REPO", unset = getwd()),
  mustWork = TRUE
)

devtools::load_all(repo_dir)

# Helper code is kept outside the coverageSim package so this technical-bias
# dataset generator does not become part of the simulator core.
source(file.path(repo_dir, "rust_helpers", "R", "G_end_dropout_helpers.R"))

input_base_raw <- Sys.getenv(
  "COVSIM_G_DROP_INPUT",
  unset = file.path("tests", "tests_and_demos", "10000a_different_gene_counts")
)
input_base <- normalizePath(input_base_raw, mustWork = TRUE)

output_parent <- normalizePath(
  Sys.getenv("COVSIM_G_DROP_OUTPUT_PARENT", unset = dirname(input_base)),
  mustWork = TRUE
)

exp_name <- "human_flavoured_riboseq"
drop_fractions <- c(0.30, 0.60, 1.00)
seed <- 42L

samtools <- Sys.which("samtools")
if (samtools == "") {
  stop("samtools is required for BAM filtering and indexing")
}

input_exp <- file.path(input_base, "experiment", paste0(exp_name, ".csv"))
if (!file.exists(input_exp)) {
  stop("Experiment CSV not found: ", input_exp)
}

input_bams <- list.files(
  file.path(input_base, "reads"),
  pattern = "\\.bam$",
  full.names = TRUE
)
if (!length(input_bams)) {
  stop("No BAM files found in: ", file.path(input_base, "reads"))
}

link_genome_files <- function(input_base, output_base) {
  # The dropout datasets reuse the same genome/annotation. Symlinks avoid
  # copying large FASTA/GTF/TxDb files for each dropout fraction.
  input_genome <- file.path(input_base, "genome")
  output_genome <- file.path(output_base, "genome")
  dir.create(output_genome, recursive = TRUE, showWarnings = FALSE)

  genome_files <- list.files(input_genome, full.names = TRUE, all.files = FALSE)
  for (src in genome_files) {
    dest <- file.path(output_genome, basename(src))
    if (file.exists(dest)) {
      next
    }
    ok <- file.symlink(normalizePath(src, mustWork = TRUE), dest)
    if (!ok) {
      stop("Could not create genome symlink: ", dest)
    }
  }
}

read_stats <- function(path) {
  x <- readLines(path)
  values <- as.integer(sub(".*\\t", "", x))
  names(values) <- sub("\\t.*", "", x)
  as.list(values)
}

filter_bam_g_end <- function(input_bam, output_bam, drop_fraction, seed, stats_file) {
  # Stream BAM -> SAM -> AWK filter -> BAM so large read files are not loaded
  # into R memory. Only reads whose sequence ends in G are randomly removed.
  awk_file <- tempfile(fileext = ".awk")
  writeLines(g_end_dropout_awk_program(), awk_file)
  on.exit(unlink(awk_file), add = TRUE)

  cmd <- paste(
    "set -euo pipefail;",
    shQuote(samtools), "view -h", shQuote(input_bam), "|",
    "awk",
    "-v", shQuote(paste0("drop=", drop_fraction)),
    "-v", shQuote(paste0("seed=", seed)),
    "-v", shQuote(paste0("stats=", stats_file)),
    "-f", shQuote(awk_file), "|",
    shQuote(samtools), "view -b -o", shQuote(output_bam), "-;",
    shQuote(samtools), "index", shQuote(output_bam)
  )

  status <- system2("bash", c("-c", cmd))
  if (!identical(status, 0L)) {
    stop("BAM filtering failed for: ", input_bam)
  }
}

summary_list <- list()
exp_lines <- readLines(input_exp)

# Create one complete dataset per dropout level so each can be analysed by the
# same RUST/plotting scripts without special-case paths.
for (drop_fraction in drop_fractions) {
  label <- g_end_dropout_label(drop_fraction)
  output_base <- file.path(output_parent, paste0(basename(input_base), "_", label))
  output_reads <- file.path(output_base, "reads")
  output_exp_dir <- file.path(output_base, "experiment")

  dir.create(output_reads, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_exp_dir, recursive = TRUE, showWarnings = FALSE)
  link_genome_files(input_base, output_base)

  message("Creating ", output_base)

  for (input_bam in input_bams) {
    run_name <- sub("\\.bam$", "", basename(input_bam))
    output_bam <- file.path(output_reads, basename(input_bam))
    stats_file <- file.path(output_reads, paste0(run_name, "_", label, "_stats.tsv"))

    if (file.exists(output_bam) || file.exists(paste0(output_bam, ".bai"))) {
      stop("Refusing to overwrite existing BAM/index: ", output_bam)
    }

    filter_bam_g_end(
      input_bam = input_bam,
      output_bam = output_bam,
      drop_fraction = drop_fraction,
      seed = seed + as.integer(round(drop_fraction * 1000)) + match(input_bam, input_bams),
      stats_file = stats_file
    )

    stats <- read_stats(stats_file)
    summary_list[[length(summary_list) + 1L]] <- data.frame(
      dataset = basename(output_base),
      drop_fraction = drop_fraction,
      run = run_name,
      total_reads = stats$total,
      g_end_reads = stats$g_end,
      removed_reads = stats$removed,
      kept_reads = stats$kept,
      stringsAsFactors = FALSE
    )
  }

  output_exp <- file.path(output_exp_dir, paste0(exp_name, ".csv"))
  output_exp_lines <- rewrite_experiment_base(exp_lines, input_base_raw, output_base)
  output_exp_lines <- rewrite_experiment_base(output_exp_lines, input_base, output_base)
  writeLines(output_exp_lines, output_exp)
}

summary_dt <- data.table::rbindlist(summary_list)
summary_file <- file.path(output_parent, paste0(basename(input_base), "_G_end_dropout_summary.csv"))
data.table::fwrite(summary_dt, summary_file)

message("Wrote summary: ", summary_file)
message("Done.")

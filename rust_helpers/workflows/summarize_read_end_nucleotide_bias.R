rm(list = ls(all.names = TRUE))
gc(reset = TRUE)

repo_dir <- normalizePath(
  Sys.getenv("COVSIM_REPO", unset = getwd()),
  mustWork = TRUE
)

devtools::load_all(repo_dir)

library(data.table)
library(Biostrings)
library(Rsamtools)

source(file.path(repo_dir, "rust_helpers", "R", "Read_end_bias_helpers.R"))

input_base <- normalizePath(
  path.expand(Sys.getenv("RUST_READ_END_INPUT", unset = "~/rust/10000a_different_gene_counts_drop_G_100")),
  mustWork = TRUE
)

output_dir <- file.path(input_base, "read_end_nt_summary")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

bam_files <- list.files(
  file.path(input_base, "reads"),
  pattern = "\\.bam$",
  full.names = TRUE
)

if (!length(bam_files)) {
  stop("No BAM files found in: ", file.path(input_base, "reads"))
}

run_ids <- sub("\\.bam$", "", basename(bam_files))
yield_size <- as.integer(Sys.getenv("RUST_READ_END_YIELD_SIZE", unset = "1000000"))

fasta_file <- Sys.getenv("RUST_READ_END_FASTA", unset = "")
if (!nzchar(fasta_file)) {
  fasta_candidates <- list.files(
    file.path(input_base, "genome"),
    pattern = "\\.(fa|fasta)$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(fasta_candidates) == 1L) {
    fasta_file <- fasta_candidates[[1]]
  }
}

if (nzchar(fasta_file)) {
  fasta_file <- normalizePath(path.expand(fasta_file), mustWork = TRUE)
  message("Using reference FASTA fallback: ", fasta_file)
} else {
  fasta_file <- NULL
  message("No reference FASTA fallback set. This is OK only if BAM SEQ fields are present.")
}

summary_dt <- summarize_bam_files_read_end_nt(
  bam_files = bam_files,
  run_ids = run_ids,
  yield_size = yield_size,
  fasta_file = fasta_file
)

csv_file <- file.path(output_dir, "read_end_nt_vs_whole_read_nt_frequency.csv")
pdf_file <- file.path(output_dir, "read_end_nt_vs_whole_read_nt_frequency.pdf")

data.table::fwrite(summary_dt, csv_file)
plot_read_end_nt_summary(summary_dt, pdf_file)

message("Saved read-end nucleotide summary CSV: ", csv_file)
message("Saved read-end nucleotide summary PDF: ", pdf_file)
message("Done.")

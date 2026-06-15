rm(list = ls(all.names = TRUE))
gc(reset = TRUE)

devtools::load_all(".")

library(coverageSim)
library(ORFik)
library(SummarizedExperiment)

# Genome-only workflow:
# - uses the real human FASTA and TxDb/GTF as the sequence/annotation substrate
# - does not read or learn from any real Ribo-seq BAM
# - simulates synthetic RFP counts and coverage on human CDSs
#
# Optional environment variables:
#   COVSIM_REPO                    coverageSim repository path
#   COVSIM_REAL_BASE               real_human_riboseq folder
#   COVSIM_HUMAN_GENOME_TOP_N        limit to top N CDSs after simulated counts
#   COVSIM_HUMAN_GENOME_TARGET_READS target total reads across all samples
#                                  default: 40000000; set to 0 to disable
#   COVSIM_HUMAN_GENOME_OUT_TAG      output folder tag

data.table::setDTthreads(4)

repo_dir <- normalizePath(
  Sys.getenv("COVSIM_REPO", unset = getwd()),
  mustWork = TRUE
)

base_dir <- normalizePath(
  Sys.getenv(
    "COVSIM_REAL_BASE",
    unset = file.path(
      repo_dir,
      "tests",
      "tests_and_demos",
      "real_human_riboseq"
    )
  ),
  mustWork = TRUE
)

genome_dir <- file.path(base_dir, "genome")

genome_file <- file.path(genome_dir, "GRCh38.primary_assembly.genome.fa")
txdb_file <- file.path(
  genome_dir,
  "gencode.v49.primary_assembly.basic.annotation.gtf.db"
)
gtf_file <- file.path(
  genome_dir,
  "gencode.v49.primary_assembly.basic.annotation.gtf"
)

if (!file.exists(gtf_file)) {
  gtf_gz_file <- paste0(gtf_file, ".gz")
  if (file.exists(gtf_gz_file)) {
    gtf_file <- gtf_gz_file
  }
}

if (!file.exists(genome_file)) stop("Genome FASTA does not exist: ", genome_file)
if (!file.exists(txdb_file)) stop("TxDb does not exist: ", txdb_file)

if (!file.exists(gtf_file)) {
  warning(
    "GTF file does not exist: ", gtf_file,
    "\nContinuing because this workflow uses the TxDb directly."
  )
}

sim_genome <- c(
  genome = genome_file,
  gtf = gtf_file,
  txdb = txdb_file
)

out_tag <- Sys.getenv(
  "COVSIM_HUMAN_GENOME_OUT_TAG",
  unset = "human_genome_only"
)

out_base <- file.path(base_dir, "coverageSim_from_human_genome_only", out_tag)
out_reads_dir <- file.path(out_base, "reads")
out_exp_dir <- file.path(out_base, "experiment")

dir.create(out_reads_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_exp_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading CDS regions from TxDb...")
cds <- ORFik::loadRegion(txdb_file, "cds")

cds_widths <- ORFik::widthPerGroup(cds, FALSE)
complete_cds <- (cds_widths %% 3 == 0) & (cds_widths > 50)

message(
  "Keeping codon-complete CDS transcripts: ",
  sum(complete_cds),
  " / ",
  length(complete_cds)
)

cds <- cds[complete_cds]
rm(cds_widths, complete_cds)
gc(reset = TRUE)

set.seed(42)

message("Simulating synthetic RFP count table on human CDSs...")
gene_count_table <- simCountTables(
  cds,
  libtypes = "RFP",
  conditions = c("WT", "KO"),
  replicates = 2,
  interceptMean = 4,
  interceptSD = 1.5,
  betaSD = 0.2,
  betaLibSD = c(RFP = 1),
  print_statistics = FALSE,
  plot_PCA = FALSE
)

top_n <- as.integer(Sys.getenv("COVSIM_HUMAN_GENOME_TOP_N", unset = "0"))

if (!is.na(top_n) && top_n > 0L && top_n < nrow(gene_count_table)) {
  message("Keeping top ", top_n, " CDSs by simulated RFP counts.")
  count_sums <- rowSums(SummarizedExperiment::assay(gene_count_table))
  keep <- order(count_sums, decreasing = TRUE)[seq_len(top_n)]
  gene_count_table <- gene_count_table[keep, ]
}

message("Distributing simulated counts to CDS regions...")
region_count_table <- simCountTablesRegions(
  count_table = gene_count_table,
  regionsToSample = "cds",
  region_proportion = list(cds = list(RFP = 1))
)

target_reads <- as.numeric(Sys.getenv(
  "COVSIM_HUMAN_GENOME_TARGET_READS",
  unset = "40000000"
))

current_reads <- sum(SummarizedExperiment::assay(region_count_table, "cds"))

if (!is.na(target_reads) && target_reads > 0 && current_reads > target_reads) {
  scale <- target_reads / current_reads
  message(
    "Downscaling simulated counts from ",
    current_reads,
    " to approximately ",
    round(target_reads),
    " total reads across all samples."
  )
  region_count_table <- scale_count_table(
    region_count_table,
    scale = scale,
    min_count = 0L
  )
}

message("Simulated reads requested per sample:")
print(colSums(SummarizedExperiment::assay(region_count_table, "cds")))
message(
  "Total simulated reads requested: ",
  sum(SummarizedExperiment::assay(region_count_table, "cds"))
)

rm(cds, gene_count_table)
gc(reset = TRUE)

exp_name <- paste0("human_genome_only_covsim_", out_tag)

message("Simulating coverage on human genome/annotation...")
sim_exp <- simNGScoverage(
  simGenome = sim_genome,
  count_table = region_count_table,
  out_dir = out_reads_dir,
  exp_name = exp_name,
  exp_save_dir = out_exp_dir,
  read_lengths_per = list(RFP = 28:30),
  ideal_coverage = list(
    cds = list(RFP = quote(rep(c(1, 0, 0), length.out = x)))
  ),
  auto_correlation = list(
    cds = list(RFP = shapes(9))
  ),
  sampling = list(
    cds = list(RFP = "DMN")
  ),
  seq_bias = load_seq_bias(type = "codon", shift = "a-site", bias = "all"),
  libFormats = list(RFP = "bam"),
  validate = TRUE,
  debug_coverage = FALSE
)

message("Simulated coverageSim BAM files:")
print(ORFik::filepath(sim_exp, "default"))
message("Experiment CSV directory: ", out_exp_dir)
message("Done.")

rm(list = ls(all.names = TRUE))
gc(reset = TRUE)

repo_dir <- normalizePath(
  Sys.getenv("COVSIM_REPO", unset = getwd()),
  mustWork = TRUE
)

devtools::load_all(repo_dir)

library(coverageSim)
library(ORFik)
# Server/full-data settings. Keep all usable CDSs and all learned counts, but
# avoid holding unnecessary large objects once each learning step is complete.
data.table::setDTthreads(4)

# Keep real-input learning outside the package namespace. These helpers are
# analysis adapters; coverageSim itself should not depend on them.
source(file.path(repo_dir, "rust_helpers", "R", "Learn_from_input.R"))

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

real_exp_dir <- file.path(base_dir, "experiment")
out_base <- file.path(base_dir, "coverageSim_from_real_input")
out_reads_dir <- file.path(out_base, "server_full", "reads")
out_exp_dir <- file.path(out_base, "server_full", "experiment")

dir.create(out_reads_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_exp_dir, recursive = TRUE, showWarnings = FALSE)

exp_name <- "human_real_riboseq"
original_exp <- file.path(real_exp_dir, paste0(exp_name, ".csv"))
stopifnot(file.exists(original_exp))

fixed_exp_dir <- tempfile("fixed-real-exp-")
dir.create(fixed_exp_dir)

exp_lines <- readLines(original_exp)
old_bases <- c(
  "~/forks/coverageSim/tests/tests_and_demos/real_human_riboseq",
  "/home/carink/forks/coverageSim/tests/tests_and_demos/real_human_riboseq",
  "~/coverageSim/tests/tests_and_demos/real_human_riboseq",
  "/home/carink/coverageSim/tests/tests_and_demos/real_human_riboseq"
)

# The experiment CSV may have been written on a laptop or server with a
# different base path. Rewrite only the known base prefixes and keep the CSV
# structure itself unchanged.
for (old_base in old_bases) {
  exp_lines <- gsub(old_base, base_dir, exp_lines, fixed = TRUE)
}

writeLines(exp_lines, file.path(fixed_exp_dir, paste0(exp_name, ".csv")))

df <- ORFik::read.experiment(
  exp_name,
  in.dir = fixed_exp_dir,
  validate = FALSE
)

sim_genome <- c(
  genome = file.path(base_dir, "genome", "GRCh38.primary_assembly.genome.fa"),
  gtf = file.path(
    base_dir,
    "genome",
    "gencode.v49.primary_assembly.basic.annotation.gtf"
  ),
  txdb = file.path(
    base_dir,
    "genome",
    "gencode.v49.primary_assembly.basic.annotation.gtf.db"
  )
)

stopifnot(all(file.exists(sim_genome)))

# Load the real BAM once. The object is reused for read-length learning,
# count-table learning, and codon-bias learning, then removed before simulation.
reads <- ORFik::fimport(ORFik::filepath(df, "default")[1])

cds <- ORFik::loadRegion(df, "cds")

tx_filt <- ORFik::filterTranscripts(
  df,
  minFiveUTR = 0,
  minCDS = 231,
  minThreeUTR = 0,
  by = "tx",
  longestPerGene = TRUE
)

cds <- cds[names(cds) %in% tx_filt]
rm(tx_filt)
gc(reset = TRUE)

# Codon-bias estimation requires complete CDS triplets. Incomplete CDS ranges
# are skipped instead of trying to infer the missing bases.
complete_cds <- ORFik::widthPerGroup(cds, FALSE) %% 3 == 0
message(
  "Keeping codon-complete CDS transcripts: ",
  sum(complete_cds),
  " / ",
  length(complete_cds)
)
cds <- cds[complete_cds]
rm(complete_cds)
gc(reset = TRUE)

# Learn expression first; geometry is calibrated only on retained CDSs below.
message("Learning CDS count table from real BAM...")
count_table <- learn_cds_count_table(
  cds,
  reads,
  sample_name = "RFP_mock_1",
  condition = "mock",
  replicate = "1",
  min_reads = 1L
)

message("CDSs with real signal used for simulation: ", nrow(count_table))

top_n <- as.integer(Sys.getenv("COVSIM_LEARN_TOP_N", unset = "0"))
if (!is.na(top_n) && top_n > 0L) {
  count_table <- subset_top_cds_counts(count_table, top_n)
  message("Keeping top learned CDSs: ", nrow(count_table))
}

target_reads <- as.numeric(Sys.getenv("COVSIM_LEARN_TARGET_READS", unset = "0"))
current_reads <- sum(SummarizedExperiment::assay(count_table, "cds"))
if (!is.na(target_reads) && target_reads > 0 && current_reads > target_reads) {
  count_table <- scale_count_table(
    count_table, scale = target_reads / current_reads, min_count = 1L
  )
}
message(
  "Simulated reads requested from learned CDS counts: ",
  sum(SummarizedExperiment::assay(count_table, "cds"))
)

cds_learned <- SummarizedExperiment::rowRanges(count_table)

geometry_max_reads <- as.integer(Sys.getenv(
  "COVSIM_LEARN_GEOMETRY_MAX_READS", unset = "200000"
))
geometry_min_per_length <- as.integer(Sys.getenv(
  "COVSIM_LEARN_GEOMETRY_MIN_PER_LENGTH", unset = "100"
))

message("Learning length-specific A-site fragment geometry...")
fragment_geometry_distribution <- learn_fragment_geometry(
  reads = reads,
  cds = cds,
  site_reference = "a_site",
  min_length = 25L,
  max_length = 34L,
  min_reads_per_length = geometry_min_per_length,
  max_observations = geometry_max_reads
)
geometry_file <- file.path(out_base, "learned_fragment_geometry.tsv")
diagnostics_file <- file.path(out_base, "learned_fragment_geometry_diagnostics.tsv")
data.table::fwrite(fragment_geometry_distribution, geometry_file, sep = "\t")
data.table::fwrite(
  attr(fragment_geometry_distribution, "offset_diagnostics"),
  diagnostics_file, sep = "\t"
)
message("Learned A-site fragment geometry:")
print(fragment_geometry_distribution)

# Geometry benefits from all qualified CDSs; codon-bias learning below is
# intentionally restricted to the retained simulation transcripts.
rm(cds)
gc(reset = TRUE)

seq_bias <- learn_codon_seq_bias(
  cds_learned,
  reads,
  fa_file = sim_genome["genome"],
  geometry_distribution = fragment_geometry_distribution,
  min_tx_reads = 20L,
  alpha_scale = 100
)

message("Top learned codon alpha values:")
print(seq_bias[order(-alpha)][1:10])

# The large read object is no longer needed after counts, fragment geometry,
# and codon bias are learned. Free it before simulation builds sequence tables.
rm(reads, cds_learned)
gc(reset = TRUE)

set.seed(42)

# From here on the learned counts, joint A-site geometry, and codon bias are
# passed to coverageSim through its regular public inputs.
sim_exp <- simNGScoverage(
  simGenome = sim_genome,
  count_table = count_table,
  out_dir = out_reads_dir,
  exp_name = "human_real_learned_covsim_server_full",
  exp_save_dir = out_exp_dir,
  read_lengths_per = list(
    RFP = fragment_geometry_distribution$fragment_length
  ),
  ideal_coverage = list(
    cds = list(RFP = quote(rep(c(1, 0, 0), length.out = x)))
  ),
  auto_correlation = list(
    cds = list(RFP = shapes(9))
  ),
  seq_bias = seq_bias,
  fragment_geometry = list(
    source = "learned",
    site_reference = "a_site",
    distribution = fragment_geometry_distribution[, .(
      fragment_length, site_offset, probability
    )],
    boundary_action = "renormalize",
    five_prime_bias = list(source = "none"),
    three_prime_bias = list(source = "none")
  ),
  ground_truth = TRUE,
  sampling = list(
    cds = list(RFP = "DMN")
  ),
  libFormats = list(RFP = "bam"),
  validate = TRUE,
  debug_coverage = FALSE
)

message("Simulated coverageSim BAM:")
print(ORFik::filepath(sim_exp, "default"))
message("Done.")

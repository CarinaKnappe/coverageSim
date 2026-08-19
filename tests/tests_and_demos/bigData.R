# Create (big) simulated data sets

# STEPS:
# Step 1: what the genes look like
# Step 2: how much each gene is expressed
# Step 3: where counts go within each transcript
# Step 4: what the read coverage looks like

# Create simulated data set: 3000 genes, A-site shifted, high-depth different gene counts

devtools::load_all(".")

library(ORFik)
library(SummarizedExperiment)

n_genes <- as.integer(Sys.getenv("COVSIM_BIGDATA_N_GENES", unset = "3000"))
if (is.na(n_genes) || n_genes < 1L) {
  stop("COVSIM_BIGDATA_N_GENES must be a positive integer")
}

run_tag <- Sys.getenv(
  "COVSIM_BIGDATA_RUN_TAG",
  unset = paste0(n_genes, "_genes_physical_asite")
)
out_root <- Sys.getenv(
  "COVSIM_BIGDATA_OUT_ROOT",
  unset = "tests/tests_and_demos"
)
out_base <- file.path(out_root, run_tag)

dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_base, "genome"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_base, "reads"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_base, "experiment"), recursive = TRUE, showWarnings = FALSE)

seed <- as.integer(Sys.getenv("COVSIM_BIGDATA_SEED", unset = "42"))
set.seed(seed)

n_reps  <- 2
conds   <- c("WT", "KO")

leader_length  <- rep(120, n_genes)
cds_length     <- rep(600, n_genes)
trailer_length <- rep(120, n_genes)
cds_intron_len <- rep(90, n_genes)

sim_genome <- simGenome(
  n = n_genes,
  out_dir = file.path(out_base, "genome"),
  genome_name = "human_flavoured_sim",
  leader_length = leader_length,
  cds_length = cds_length,
  trailer_length = trailer_length,
  cds_exons = 2,
  cds_intron_length = cds_intron_len,
  # uORFs are deliberately enabled so their geometry and 10% signal fraction
  # are exercised in both the smoke test and the full run.
  max_uorfs = 2,
  uorfs_can_overlap = TRUE,
  uorfs_can_overlap_cds = 0,
  uorf_max_length = 90,
  export_txdb = TRUE,
  debug_on = FALSE
)

cds <- ORFik::loadRegion(sim_genome["txdb"], "cds")

gene_counts <- simCountTables(
  n = cds,
  libtypes = "RFP",
  conditions = conds,
  replicates = n_reps,
  # Higher depth for RUST: 4 was too sparse after transcript QC.
  # Keeping the variance/effect settings unchanged isolates read-depth effects.
  interceptMean = 6,
  interceptSD = 1.5,
  betaSD = 0.2,
  plot_PCA = FALSE,
  print_statistics = TRUE
)

sample_totals <- colSums(SummarizedExperiment::assay(gene_counts))
print(sample_totals)

region_counts <- simCountTablesRegions(
  count_table = gene_counts,
  regionsToSample = c("leader", "cds", "trailer", "uorf"),
  region_proportion = list(
    leader  = list(RFP = 0.03),
    cds     = list(RFP = 0.84),
    trailer = list(RFP = 0.03),
    uorf    = list(RFP = 0.10)
  )
)

region_totals <- data.frame(
  region = c("leader", "cds", "trailer", "uorf"),
  do.call(rbind, lapply(
    c("leader", "cds", "trailer", "uorf"),
    function(region) {
      colSums(SummarizedExperiment::assay(region_counts, region))
    }
  )),
  check.names = FALSE
)
colnames(region_totals)[-1L] <- colnames(region_counts)
data.table::fwrite(
  region_totals,
  file.path(out_base, "simulated_region_counts.tsv"),
  sep = "\t"
)

exp <- simNGScoverage(
  simGenome = sim_genome,
  count_table = region_counts,
  out_dir = file.path(out_base, "reads"),
  exp_name = "human_flavoured_riboseq",
  exp_save_dir = file.path(out_base, "experiment"),
  read_lengths_per = list(RFP = 28:30),
  seq_bias = load_seq_bias(
    type = "codon",
    shift = "a-site",
    bias = "all"
  ),
  auto_correlation = list(
    cds = list(RFP = shapes(9)),
    uorf = list(RFP = shapes(9))
  ),
  fragment_geometry = list(
    source = "default",
    site_reference = "a_site",
    boundary_action = "renormalize",
    five_prime_bias = list(source = "none"),
    three_prime_bias = list(source = "none")
  ),
  ground_truth = TRUE,
  libFormats = list(RFP = "bam"),
  validate = TRUE,
  debug_coverage = FALSE
)

bam_files <- ORFik::filepath(exp, "default")
print(bam_files)

message("Ground-truth files:")
print(list.files(file.path(out_base, "reads"), pattern = "ground.truth", full.names = TRUE))

count_mat <- SummarizedExperiment::assay(gene_counts)
sample_info <- as.data.frame(SummarizedExperiment::colData(gene_counts))

grDevices::pdf(file.path(out_base, "simulation_diagnostics.pdf"))

wt_cols <- which(sample_info$condition == "WT")
ko_cols <- which(sample_info$condition == "KO")

wt_mean <- rowMeans(count_mat[, wt_cols, drop = FALSE])
ko_mean <- rowMeans(count_mat[, ko_cols, drop = FALSE])

plot(
  wt_mean, ko_mean,
  pch = 16, cex = 0.5, col = rgb(0, 0, 1, 0.25),
  xlab = "Mean WT gene count",
  ylab = "Mean KO gene count",
  main = "Gene-level Ribo-seq counts: WT vs KO"
)
abline(0, 1, col = "red", lty = 2)

region_plot_totals <- rbind(
  leader  = colSums(SummarizedExperiment::assay(region_counts, "leader")),
  cds     = colSums(SummarizedExperiment::assay(region_counts, "cds")),
  trailer = colSums(SummarizedExperiment::assay(region_counts, "trailer")),
  uorf    = colSums(SummarizedExperiment::assay(region_counts, "uorf"))
)

region_props <- sweep(
  region_plot_totals,
  2,
  colSums(region_plot_totals),
  "/"
)

barplot(
  region_props,
  beside = FALSE,
  col = c("steelblue", "tomato", "goldenrod", "darkseagreen4"),
  names.arg = paste(sample_info$condition, sample_info$replicate, sep = "_"),
  las = 2,
  ylab = "Fraction of counts",
  main = "Region composition across samples"
)

legend(
  "topright",
  legend = rownames(region_props),
  fill = c("steelblue", "tomato", "goldenrod", "darkseagreen4"),
  bty = "n"
)

grDevices::dev.off()

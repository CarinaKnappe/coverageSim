# Load the local package version, not an installed one
devtools::load_all(".")
library(coverageSim)
library(ORFik)
library(SummarizedExperiment)
library(ggplot2)

set.seed(1)

# STEPS:
# Step 1: what the genes look like
# Step 2: how much each gene is expressed
# Step 3: where counts go within each transcript
# Step 4: what the read coverage looks like

# ---------------------------
# STEP 1: Create a synthetic genome / what the genes look like
# User settings are mainly here:
# - n: number of genes
# - leader_length / cds_length / trailer_length: transcript structure
# - cds_exons / cds_intron_length: exon-intron structure (pls mind: very short leaders can make uORF sampling hard or impossible))
# - max_uorfs: max number of uORFs per transcript
# - uorfs_can_overlap / uorfs_can_overlap_cds: overlap rules
#
# Alternatively: you can provide your own FASTA/GTF files and TxDb object instead of simulating a genome.
# ---------------------------

n_genes = 1200
sim_genome <- simGenome(
  n = n_genes,                              # how many genes
  out_dir = tempdir(),                 # where FASTA/GTF/TxDb are written
  genome_name = "my_demo_genome",
  leader_length = rep(80, n_genes),         # 5' UTR length (must be divisible by 3 if uORFs are simulated)
  cds_length = rep(300, n_genes),           # CDS length (must be divisible by 3)
  trailer_length = rep(120, n_genes),       # 3' UTR length
  cds_exons = 2,                       # number of CDS exons
  cds_intron_length = rep(60, n_genes),     # intron length between CDS exons
  max_uorfs = 2,                       # max uORFs per transcript
  uorfs_can_overlap = TRUE,            # uORFs allowed to overlap each other?
  uorfs_can_overlap_cds = 1,            # 0 = no CDS overlap, 1 = can overlap, 2 = must overlap
  debug_on = TRUE
  )

# Plot after step 1:
# Plot the 12 simulated genes as horizontal region tracks
txdb <- loadTxdb(sim_genome["txdb"])
leaders  <- loadRegion(txdb, "leaders")
cds      <- loadRegion(txdb, "cds")
trailers <- loadRegion(txdb, "trailers")
par(mar = c(5, 8, 4, 3))
plot(
  NULL,
  xlim = c(1, max(c(widthPerGroup(leaders, FALSE),
                    widthPerGroup(cds, FALSE),
                    widthPerGroup(trailers, FALSE))) + 250),
  ylim = c(0.5, n_genes+0.5),
  xlab = "Transcript position (nt)",
  ylab = "Gene",
  yaxt = "n",
  main = "Step 1: Simulated transcript models"
)

axis(2, at = 1:n_genes, labels = names(leaders), las = 2)

for (i in seq_len(n_genes)) {
  y <- i
  l_len <- width(leaders[i])
  c_len <- widthPerGroup(cds[i], FALSE)
  t_len <- width(trailers[i])

  rect(1, y - 0.25, l_len, y + 0.25, col = "steelblue", border = NA)
  rect(l_len + 1, y - 0.25, l_len + c_len, y + 0.25, col = "tomato", border = NA)
  rect(l_len + c_len + 1, y - 0.25, l_len + c_len + t_len, y + 0.25,
       col = "goldenrod", border = NA)
}

legend(
  "topright",
  legend = c("leader", "CDS", "trailer"),
  fill = c("steelblue", "tomato", "goldenrod"),
  bty = "n"
)

# ---------------------------
# STEP 2: Simulate gene-level counts / how much each gene is expressed in each sample?
# simCountTables asks: for each gene, how many counts would I observe in each sample?
# User settings are mainly here:
# - libtypes, conditions, replicates: experimental design
# - interceptMean: overall expression level. Impacts the number of gene counts: higher = more counts overall
#   controls how many counts a gene tends to get on average across samples.
# - interceptSD: between-gene variability, higher = more variation between genes
# - betaSD: between-condition effect size, higher = stronger
# ---------------------------
# Load CDS annotation from the simulated TxDb
cds <- ORFik::loadRegion(sim_genome["txdb"], "cds")
# (head(cds))

gene_counts <- simCountTables(
  n = cds,
  libtypes = "RFP",                    # only ribosome profiling here (other options: RNA, TCP, RCP)
  conditions = c("WT", "KO"),          # two conditions / groups
  replicates = 2,                      # two replicates / repeated samplings per each group
  interceptMean = 8,                   # overall expression level of genes
  interceptSD = 3.5,                   # expression level variation between genes
  betaSD = 1.8,                        # condition differences
  plot_PCA = FALSE,
  print_statistics = FALSE
)

# Plot after step 2: library sizes
par(mar = c(5, 8, 4, 3))
barplot(
  colSums(assay(gene_counts)),
  main = "Step 2: Gene-level library sizes",
  ylab = "Total counts",
  las = 2,
  col = "grey40"
)

# plot mean gene count per group:
par(mar = c(5, 8, 4, 2))
count_mat <- assay(gene_counts)
sample_info <- as.data.frame(colData(gene_counts))

wt_cols <- which(sample_info$condition == "WT")
ko_cols <- which(sample_info$condition == "KO")

wt_mean <- rowMeans(count_mat[, wt_cols, drop = FALSE])
ko_mean <- rowMeans(count_mat[, ko_cols, drop = FALSE])

plot(
  wt_mean, ko_mean,
  pch = 16,
  col = "steelblue",
  xlab = "Mean WT gene count", # points above the line: higher in KO
  ylab = "Mean KO gene count", # points below the line: higher in WT
  main = "Step 2: Gene counts per group"
)
abline(0, 1, col = "red", lty = 2) #

# -> now we have: how many counts per gene in each sample,
# but we don't know where those counts go within the transcript (leader, cds, trailer, uorf)
# (...and what the read coverage looks like. We will define that in the next steps.)

# ---------------------------
# STEP 3: Split gene counts across regions / where counts go within each transcript
# User settings are mainly here:
# - regionsToSample: which transcript regions get signal
# - region_proportion: expected fraction of counts in each region
# ---------------------------
region_counts <- simCountTablesRegions(
  count_table = gene_counts,
  regionsToSample = c("leader", "cds", "trailer", "uorf"),
  region_proportion = list(
    leader  = list(RFP = 0.10), # 10% of counts of ribosome footprints in the leader
    cds     = list(RFP = 0.75), # 75% of counts of ribosome footprints in the CDS
    trailer = list(RFP = 0.05), # 05% of counts of ribosome footprints in the trailer
    uorf    = list(RFP = 0.10)  # 10% of counts of ribosome footprints in uORFs
  )
)

# Plot after step 3: region composition for the first sample
sample1_regions <- c(
  leader  = sum(assay(region_counts, "leader")[, 1]),
  cds     = sum(assay(region_counts, "cds")[, 1]),
  trailer = sum(assay(region_counts, "trailer")[, 1]),
  uorf    = sum(assay(region_counts, "uorf")[, 1])
)

barplot(
  sample1_regions,
  main = "Step 3: Region-level counts in sample 1",
  ylab = "Counts",
  col = c("steelblue", "tomato", "goldenrod", "darkseagreen4")
)

# plot for all samples per group:
# sample_info <- as.data.frame(colData(region_counts))
# sample_names <- colnames(region_counts)
#
# region_totals <- rbind(
#   leader  = colSums(assay(region_counts, "leader")),
#   cds     = colSums(assay(region_counts, "cds")),
#   trailer = colSums(assay(region_counts, "trailer")),
#   uorf    = colSums(assay(region_counts, "uorf"))
# )
#
# region_props <- sweep(region_totals, 2, colSums(region_totals), "/")
#
# barplot(
#   region_props,
#   beside = FALSE,
#   col = c("steelblue", "tomato", "goldenrod", "darkseagreen4"),
#   names.arg = paste(sample_info$condition, sample_info$replicate, sep = "_"),
#   las = 2,
#   ylab = "Fraction of counts",
#   main = "Step 3: Region composition across all samples"
# )
#
# legend(
#   "topright",
#   legend = rownames(region_props),
#   fill = c("steelblue", "tomato", "goldenrod", "darkseagreen4"),
#   bty = "n"
# )
# look similar because the groups are still expected to have broadly similar fractions in leader, cds, trailer, uorf (however can differ in overall gene abundance)

# ---------------------------
# STEP 4: Simulate read coverage / what the read coverage looks like
# User settings are mainly here:
# - read_lengths_per: read lengths
# - libFormats: output format
# - ideal_coverage, rnase_bias, auto_correlation, seq_bias: coverage shape
# ---------------------------
exp <- simNGScoverage(
  simGenome = sim_genome,
  count_table = region_counts[, 1],    # simulate one sample to keep it simple!
  # count_table = region_counts,       # simulate all
  exp_name = "demo_sim",
  exp_save_dir = tempdir(),
  read_lengths_per = list(RFP = 28:30),
  libFormats = list(RFP = "ofst"),     # export in OFST format (other options: "bam", "fastq")
  validate = TRUE                      # consistent with ORFik?
)

# Import simulated reads and plot coverage on the first CDS
reads <- ORFik::fimport(ORFik::filepath(exp, "default")[1])
tx1 <- cds[1]
# tx1 <- leaders[1]

cov_dt <- ORFik::coveragePerTiling(
  tx1,
  reads,
  is.sorted = TRUE,
  as.data.table = TRUE
)
#
 plot(
   cov_dt$position, cov_dt$count,
   type = "h",
   main = "Step 4: Simulated coverage on first CDS",
   xlab = "Position in CDS",
   ylab = "Coverage"
 )
#
# # PLOT WITH FRAMES
 cov_dt[, frame := position %% 3]
 ggplot(cov_dt, aes(x = position, y = count, fill = as.factor(frame))) +
   geom_col(position = 'identity')
#
# # PLOT WITH FRAMES 0,1,2:
cov_dt[, frame := factor(position %% 3, levels = c(0, 1, 2),
                          labels = c("Frame 0", "Frame 1", "Frame 2"))]

ggplot(cov_dt, aes(x = position, y = count, fill = frame)) +
  geom_col(position = "identity") +
  labs(
    title = "Step 4: Simulated coverage on first CDS by frame",
    x = "Position in CDS",
    y = "Coverage",
    fill = "Reading frame"
  )

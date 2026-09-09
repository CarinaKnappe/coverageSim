devtools::load_all("~/forks/coverageSim")
library(coverageSim)
library(ORFik)
library(SummarizedExperiment)

genome_dir <- "~/forks/coverageSim_tests/tests/individual_tests/genome"
exp_dir <- "~/forks/coverageSim_tests/tests/individual_tests/experiments"

dir.create(genome_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(exp_dir, showWarnings = FALSE, recursive = TRUE)

exp_name <- paste0("stress_test_", format(Sys.time(), "%Y%m%d_%H%M%S"))

set.seed(1)

cat("Starting user stress test\n")

t1 <- system.time({
  sim_genome <- simGenome(
    n = 100,
    max_uorfs = 2,
    cds_length = sample(c(300, 330, 360), 100, replace = TRUE),
    out_dir = genome_dir,
    genome_name = "stress_test_genome",
    debug_on = FALSE
  )
})

cat("simGenome done\n")
print(t1)

t2 <- system.time({
  gene_count_table <- simCountTables(
    loadRegion(sim_genome["txdb"], "cds"),
    libtypes = "RFP",
    print_statistics = FALSE,
    plot_PCA = FALSE,
    interceptMean = 10
  )
})

cat("simCountTables done\n")
print(t2)

gene_reads_per_sample <- colSums(assay(gene_count_table))
cat("\nExpected ribo-seq reads from count table\n")
print(gene_reads_per_sample)
cat("Total expected ribo-seq reads:", sum(gene_reads_per_sample), "\n")

t3 <- system.time({
  region_count_table <- simCountTablesRegions(
    gene_count_table,
    regionsToSample = c("cds", "uorf")
  )
})

cat("simCountTablesRegions done\n")
print(t3)

region_reads_per_sample <- colSums(
  assay(region_count_table, "cds") + assay(region_count_table, "uorf")
)
cat("\nExpected sampled reads across tested regions (cds + uorf)\n")
print(region_reads_per_sample)
cat("Total expected sampled reads:", sum(region_reads_per_sample), "\n")

t4 <- system.time({
  experiment <- simNGScoverage(
    sim_genome,
    region_count_table,
    exp_name = exp_name,
    exp_save_dir = exp_dir,
    validate = TRUE
  )
})

cat("simNGScoverage done\n")
print(t4)

default_files <- ORFik::filepath(experiment, "default")
written_reads_per_sample <- setNames(
  sapply(default_files, function(path) length(ORFik::fimport(path))),
  basename(default_files)
)

cat("\nImported reads from written coverage files\n")
print(written_reads_per_sample)
cat("Total imported reads:", sum(written_reads_per_sample), "\n")

t5 <- system.time({
  ORFik::convert_to_bigWig(experiment)
})

cat("convert_to_bigWig done\n")
print(t5)

cat("\nOutput files\n")
print(ORFik::filepath(experiment, "default"))
print(ORFik::filepath(experiment, "bigwig"))

cat("\nTotal timing summary\n")
print(rbind(
  simGenome = t1,
  simCountTables = t2,
  simCountTablesRegions = t3,
  simNGScoverage = t4,
  convert_to_bigWig = t5
))

cat("\nTest summary\n")
cat("Libraries tested:", length(default_files), "\n")
cat("Total ribo-seq reads simulated:", sum(gene_reads_per_sample), "\n")
cat("Total ribo-seq reads written:", sum(written_reads_per_sample), "\n")

experiment

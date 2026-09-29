# Local, reproducible benchmark; source from the coverageSim repository.
devtools::load_all(".")
library(ORFik)
library(data.table)
set.seed(20260909)
root <- file.path(getwd(), "rust_helpers", "runs", "2026-09-09_300tx_codon-bias_run2")
if (dir.exists(root)) stop("Output already exists: ", root)
dir.create(root, recursive = TRUE)
genome <- simGenome(n = 300L, out_dir = file.path(root, "genome"),
  genome_name = "rust_large", cds_length = rep(c(900L, 1200L, 1500L), 100L),
  leader_length = rep(60L, 300L), trailer_length = rep(60L, 300L),
  cds_exons = 2L, cds_intron_length = rep(30L, 300L),
  max_uorfs = 0L, debug_on = FALSE)
cds <- loadRegion(genome["txdb"], "cds")
counts <- SummarizedExperiment::SummarizedExperiment(
  assays = list(gene = matrix(5000L, length(cds), 1L,
    dimnames = list(names(cds), "RFP_benchmark"))), rowRanges = cds,
  colData = S4Vectors::DataFrame(libtype = factor("RFP"),
    condition = factor("benchmark"), replicate = "1", row.names = "RFP_benchmark"))
counts <- simCountTablesRegions(counts, regionsToSample = "cds",
  region_proportion = list(cds = list(RFP = 1)), sampling = c(RFP = "MN"))
targets <- c("AAA", "GAA", "CCA", "CGT", "TTC", "GGT")
profile <- data.table(seqs = names(Biostrings::GENETIC_CODE), alpha = 1000)
scenarios <- c("control_28", "strong_28", "strong_mixed")
Sys.setenv(RUST_AUTORUN = "false", RUST_PLOT_AUTORUN = "false")
rust <- new.env(); sys.source("/home/carink/rust_scripts/RUST/run_rust.R", rust)
plots <- new.env(); sys.source("/home/carink/rust_scripts/RUST/plot_rust.R", plots)
for (scenario in scenarios) {
  folder <- file.path(root, scenario)
  dir.create(file.path(folder, "experiment"), recursive = TRUE)
  dir.create(file.path(folder, "reads"), recursive = TRUE)
  bias <- copy(profile)
  if (scenario != "control_28") bias[seqs %in% targets, alpha := 20000]
  distribution <- if (scenario == "strong_mixed") {
    data.frame(fragment_length = 27:29, site_offset = c(13L, 15L, 17L), probability = rep(1/3, 3))
  } else data.frame(fragment_length = 28L, site_offset = 15L, probability = 1)
  experiment <- simNGScoverage(simGenome = genome, count_table = counts,
    out_dir = file.path(folder, "reads"), exp_name = scenario,
    exp_save_dir = file.path(folder, "experiment"),
    ideal_coverage = list(cds = list(RFP = quote(rep(1, x)))),
    rnase_bias = list(RFP = 1), auto_correlation = list(cds = list(RFP = NULL)),
    read_lengths_per = list(RFP = distribution$fragment_length),
    sampling = list(cds = list(RFP = "DMN")), seq_bias = bias,
    fragment_geometry = list(source = "user", site_reference = "a_site",
      distribution = distribution, boundary_action = "renormalize",
      five_prime_bias = list(source = "none"), three_prime_bias = list(source = "none")),
    ground_truth = TRUE, libFormats = list(RFP = "bam"), validate = TRUE)
  fwrite(bias, file.path(folder, "input_codon_weights.tsv"), sep = "\t")
  offsets <- if (scenario == "strong_mixed") c(correct = "27:13,28:15,29:17", common15 = "27:15,28:15,29:15") else c(correct = "28:15")
  for (analysis in names(offsets)) {
    Sys.setenv(RUST_SIM_BASE = folder, RUST_EXPERIMENT = scenario,
      RUST_OUTPUT_DIR = file.path(folder, paste0("rust_", analysis)),
      RUST_SITE_OFFSETS = offsets[[analysis]], RUST_MIN_READS = "100",
      RUST_MIN_COVERAGE = "0.1", RUST_OVERWRITE = "false")
    rust$main()
    plots$main()
  }
  gc()
}
writeLines(c("300 transcripts; CDS body lengths 900/1200/1500 nt; 5,000 reads per transcript.",
  "Seed 20260909. DMN sampling, codon alpha 1000 baseline / 20000 targets.",
  paste("Target codons:", paste(targets, collapse = ", ")),
  "No autocorrelation or end bias. Full fragments; length-specific A-site offsets.",
  "Strong mixed BAM analysed twice, with correct offsets and common offset 15.",
  "RUST filters: minCDS 231; trim first 120/last 60 nt; >=100 reads; >=0.1 covered fraction.",
  "Synthetic positive control, not calibrated biological bias. One replicate per scenario."),
  file.path(root, "README.txt"))
writeLines(capture.output(sessionInfo()), file.path(root, "sessionInfo.txt"))

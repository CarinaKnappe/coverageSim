devtools::load_all(".")
library(ORFik)
library(data.table)
source("analysis/R/end_bias_library_helpers.R")
data_root <- Sys.getenv("COVSIM_DATA_ROOT", unset = file.path(dirname(getwd()), "coverageSim_data"))
root <- file.path(data_root, "runs", "2026-09-10_original-rust_read-ends")
if (dir.exists(root)) stop("Output already exists: ", root)
dir.create(root, recursive = TRUE)
source_root <- file.path(data_root, "runs", "2026-09-09_300tx_codon-bias_run2")
fa <- file.path(source_root, "genome/rust_large.fasta")
db <- file.path(source_root, "genome/rust_large.gtf.db")
mrna <- loadRegion(db, "mrna")
cds <- loadRegion(db, "cds")
models <- transcript_models(mrna, fa)
truth <- fread(file.path(source_root, "control_28/reads/RFP_benchmark_ground_truth.tsv"))
stopifnot(sum(truth$score) == 1500000L, all(truth$fragment_length == 28L), all(truth$site_offset == 15L))
truth <- project_fragment_truth(truth, models)
annotation <- rbindlist(lapply(names(models), function(id) {
  model <- models[[id]]
  positions <- transcript_genomic_positions(model)
  cds_positions <- unlist(lapply(seq_along(cds[[id]]), function(i) {
    seq.int(start(cds[[id]])[i], end(cds[[id]])[i])
  }), use.names = FALSE)
  tx_cds <- match(cds_positions, positions)
  stopifnot(!anyNA(tx_cds))
  data.table(transcript_id = id, length = model$length,
             cds_start = min(tx_cds) - 1L, cds_end = max(tx_cds), strand = model$strand)
}))
fwrite(annotation, file.path(root, "transcript_annotation.tsv"), sep = "\t")
# Original RUST uses Python slices: zero-based CDS start and exclusive end.
fasta <- file(file.path(root, "transcripts.fa"), "w")
for (i in seq_len(nrow(annotation))) {
  row <- annotation[i]
  writeLines(c(paste0(">", row$transcript_id, "\t", row$cds_start, "\t", row$cds_end),
               models[[row$transcript_id]]$sequence), fasta)
}
close(fasta)
scenarios <- data.table(scenario = c("control", "five_prime", "three_prime", "both_ends"),
  five_factor = c(1,4,1,4), three_factor = c(1,1,4,4), seed = 20260910L + 0:3)
fwrite(scenarios, file.path(root, "scenario_settings.tsv"), sep = "\t")
summary <- list()
for (i in seq_len(nrow(scenarios))) {
  setting <- scenarios[i]
  dir <- file.path(root, setting$scenario)
  dir.create(dir)
  reads <- resample_end_library(truth, setting$five_factor, setting$three_factor,
                                reads_per_tx = 5000L, seed = setting$seed)
  reads[, first_nt := substr(sequence, 1L, 1L)]
  reads[, last_nt := substr(sequence, nchar(sequence), nchar(sequence))]
  fwrite(reads, file.path(dir, "fragment_truth.tsv"), sep = "\t")
  summary[[i]] <- reads[, .(scenario = setting$scenario, total_reads = sum(score),
    transcripts = uniqueN(transcript_id), five_G_fraction = sum(score[first_nt == "G"])/sum(score),
    three_C_fraction = sum(score[last_nt == "C"])/sum(score),
    both_GC_fraction = sum(score[first_nt == "G" & last_nt == "C"])/sum(score))]
}
fwrite(rbindlist(summary), file.path(root, "end_composition.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(root, "sessionInfo.txt"))
message("End-bias library tables ready: ", root)

# Fresh CovSim count draws and codon-weight model for the local fragment workflow.
devtools::load_all(".")
library(ORFik)
library(data.table)
root <- file.path(getwd(), "analysis/benchmarks/rust/runs/2026-09-10_natural-read-ends")
if (dir.exists(root)) stop("Output exists: ", root)
dir.create(root, recursive = TRUE)
previous <- file.path(getwd(), "analysis/benchmarks/rust/runs/2026-09-10_original-rust_read-ends")
for (name in c("transcripts.fa", "transcript_annotation.tsv")) {
  stopifnot(file.copy(file.path(previous, name), file.path(root, name)))
}
genome_dir <- file.path(getwd(), "analysis/benchmarks/rust/runs/2026-09-09_300tx_codon-bias_run2/genome")
cds <- loadRegion(file.path(genome_dir, "rust_large.gtf.db"), "cds")
set.seed(202609101L)
# Two NB draws with identical underlying expression: betaSD=0.
counts <- simCountTables(n = cds, libtypes = "RFP", conditions = c("draw1", "draw2"),
  replicates = 1L, betaSD = 0, betaLibSD = c(RFP = 0),
  interceptMean = 10, interceptSD = 2.1, dispMeanRel = function(x) 4/x + .12,
  plot_PCA = FALSE, print_statistics = FALSE)
raw <- SummarizedExperiment::assay(counts)
allocation <- function(x, n = 1500000L) {
  expected <- (x + 1) / sum(x + 1) * n
  out <- floor(expected)
  remaining <- n - sum(out)
  if (remaining > 0) {
    idx <- order(expected - out, decreasing = TRUE)[seq_len(remaining)]
    out[idx] <- out[idx] + 1L
  }
  as.integer(out)
}
count_dt <- data.table(transcript_id = names(cds), raw_rep1 = raw[,1], raw_rep2 = raw[,2],
                       rep1 = allocation(raw[,1]), rep2 = allocation(raw[,2]))
fwrite(count_dt, file.path(root, "transcript_counts.tsv"), sep = "\t")
profile <- load_seq_bias(type = "codon", shift = "a-site", bias = "R10")
fwrite(profile, file.path(root, "R10_source_profile.tsv"), sep = "\t")
regular <- profile[grepl("^[ACGT]{3}$", seqs)]
regular[, relative_weight := pmax(.6, pmin(1.8, (alpha / median(alpha))^.35))]
fwrite(regular, file.path(root, "applied_codon_weights.tsv"), sep = "\t")
sequences <- as.character(txSeqsFromFa(cds, file.path(genome_dir,"rust_large.fasta")))
annotation <- fread(file.path(root, "transcript_annotation.tsv"))
site_tables <- lapply(seq_along(cds), function(i) {
  seq <- sequences[i]
  codons <- substring(seq, seq.int(1L,nchar(seq),3L), seq.int(3L,nchar(seq),3L))
  weight <- regular$relative_weight[match(codons, regular$seqs)]
  stopifnot(all(codons[is.na(weight)] %in% c("TAA","TAG","TGA")))
  weight[is.na(weight)] <- 1
  nt_weights <- coverageSim:::sim_sequence_bias(quote(rep(1,x)), nchar(seq),
    list(weight), seq_acf = NULL, rnase_acf = NULL)[[1]]
  codon_weight <- nt_weights[seq.int(1L,length(nt_weights),3L)]
  start0 <- annotation$cds_start[match(names(cds)[i], annotation$transcript_id)]
  data.table(transcript_id = names(cds)[i], site0 = start0 + seq.int(0L,nchar(seq)-1L,3L),
    codon = codons, weight = codon_weight / mean(codon_weight))
})
fwrite(rbindlist(site_tables), file.path(root, "covsim_site_weights.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(root, "sessionInfo.txt"))
message("Fresh CovSim count and coverage inputs ready: ", root)

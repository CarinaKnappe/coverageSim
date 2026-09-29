# Prepare unique real-human Ribo-seq reads for CHOROS.
source_file <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = TRUE), error = function(e) NA_character_)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (!is.na(source_file)) dirname(source_file) else if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
source(file.path(script_dir, "choros_utils.R"))

library(data.table)
library(GenomicFeatures)
library(ORFik)
library(Biostrings)

base_dir <- Sys.getenv("CHOROS_REAL_BASE", unset = "")
if (!nzchar(base_dir)) stop("Set CHOROS_REAL_BASE to real_human_riboseq.")
base_dir <- normalizePath(base_dir, mustWork = TRUE)
output_dir <- Sys.getenv("CHOROS_OUTPUT_DIR", unset = file.path(base_dir, "choros_real_human"))
prepared_dir <- file.path(output_dir, "prepared")
dir.create(prepared_dir, recursive = TRUE, showWarnings = FALSE)

run_id <- "SRR32491292"
ofst_file <- file.path(base_dir, "reads", paste0(run_id, ".sorted.unique_nh1.ofst"))
bam_file <- file.path(base_dir, "reads", paste0(run_id, ".sorted.unique_nh1.bam"))
fasta_file <- file.path(base_dir, "genome", "GRCh38.primary_assembly.genome.fa")
txdb_file <- file.path(base_dir, "genome", "gencode.v49.primary_assembly.basic.annotation.gtf.db")
required <- c(OFST = ofst_file, BAM = bam_file, FASTA = fasta_file, TxDb = txdb_file)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing input files:\n", paste(missing, collapse = "\n"))

experiment_dir <- tempfile("choros-real-human-")
dir.create(experiment_dir)
ORFik::create.experiment(
  dir = dirname(bam_file), exper = "choros_real_human", txdb = txdb_file,
  fa = fasta_file, saveDir = experiment_dir, organism = "Homo sapiens",
  author = "CHOROS real-human control", libtype = "RFP",
  condition = "mock", rep = "1", files = bam_file
)
df <- ORFik::read.experiment("choros_real_human", in.dir = experiment_dir, validate = FALSE)
uniqueMappers(df) <- TRUE

message("Loading and aligning transcript annotations...")
tx <- loadRegion(df, "tx")
cds <- loadRegion(df, "cds")
leaders <- loadRegion(df, "leaders")
keep <- filterTranscripts(df, minFiveUTR = 1L, minCDS = 231L,
                          minThreeUTR = 0L, by = "tx", longestPerGene = TRUE)
keep <- keep[keep %in% names(tx) & keep %in% names(cds) & keep %in% names(leaders)]
tx <- set_unknown_circular_to_false(tx[keep])
cds <- cds[keep]
leaders <- leaders[keep]
if (!length(tx)) stop("No aligned transcripts remain after filtering.")
tx_seqs <- as.character(ORFik::txSeqsFromFa(tx, fasta_file))
utr5_lengths <- sum(width(leaders)); cds_lengths <- sum(width(cds))
names(utr5_lengths) <- names(cds_lengths) <- names(tx)

message("Importing unique-mapper OFST and mapping 5' ends...")
reads <- set_unknown_circular_to_false(fimport(ofst_file))
mcols(reads)$L <- as.integer(readWidths(reads, after.softclips = TRUE))
if (!"score" %in% names(mcols(reads))) mcols(reads)$score <- 1
reads_5p <- resize(as(reads, "GRanges"), width = 1L, fix = "start")
mapped <- map_read_ends_to_transcripts(reads_5p, tx)
tx_id <- names(tx)[mapped$transcript_index]
dt_reads <- data.table(
  transcript_id = tx_id, tx_5p_pos = as.integer(start(mapped$mapped)),
  L = as.integer(mcols(mapped$reads)$L), score = as.numeric(mcols(mapped$reads)$score),
  utr5_len = as.integer(utr5_lengths[tx_id]), cds_len = as.integer(cds_lengths[tx_id])
)
dt_reads <- dt_reads[complete.cases(dt_reads)]
dt_reads[, dist_start := tx_5p_pos - utr5_len - 1L]

periodicity <- select_periodic_lengths(dt_reads, candidate_lengths = 25:30,
                                       min_reads = 1000L, min_frame_prop = 0.5)
message("Accepted periodic lengths: ", paste(periodicity$lengths, collapse = ", "))
tis <- infer_tis_offsets(dt_reads, periodicity$lengths,
                         tis_window = c(-40L, -3L), min_peak_reads = 50L)
off_long <- build_frame_offset_map(tis$offsets)
dt_reads[, pos_frame := dist_start %% 3L]
dt_reads <- merge(dt_reads, off_long, by.x = c("L", "pos_frame"),
                  by.y = c("qwidth", "pos_frame"), all = FALSE)
dt_reads[, A_start := tx_5p_pos + d5]
dt_reads[, d3 := L - d5 - 3L]
dt_reads[, cod_idx := (A_start - utr5_len + 2L) %/% 3L]
dt_reads <- dt_reads[d5 >= 0L & d3 >= 0L & cod_idx >= 1L & cod_idx <= cds_len %/% 3L]
if (!nrow(dt_reads)) stop("No reads have valid A-site geometry after QC.")
dt_obs <- dt_reads[, .(count = sum(score)), by = .(transcript_id, cod_idx, d5, d3, utr5_len, cds_len)]

prepared <- list(run_id = run_id, dt_obs = dt_obs, tx_seqs = tx_seqs,
                 off_long = off_long, offsets = tis$offsets,
                 periodicity_qc = periodicity$qc, tis_profile = tis$profile,
                 source_ofst = normalizePath(ofst_file))
prepared_file <- file.path(prepared_dir, paste0(run_id, "_choros_input.rds"))
saveRDS(prepared, prepared_file)
fwrite(tis$offsets, file.path(prepared_dir, paste0(run_id, "_offsets.csv")))
fwrite(periodicity$qc, file.path(prepared_dir, paste0(run_id, "_periodicity_qc.csv")))
fwrite(tis$profile, file.path(prepared_dir, paste0(run_id, "_tis_profile.csv")))
message("Prepared CHOROS input: ", prepared_file)

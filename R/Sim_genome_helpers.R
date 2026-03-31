genome_input_test_controller <- function() {
  with(rlang::caller_env(),{
    # Sanity tests for input
    message("-- Verifying input")
    stopifnot(n > 0)
    stopifnot(all(cds_intron_length%%3 == 0))
    stopifnot(all(cds_length%%3 == 0))
    stopifnot(length(strand) == n)
    stopifnot(length(leader_length) == n)
    stopifnot(length(trailer_length) == n)
    stopifnot(length(seqnames) == length(unique(seqnames)))
    stopifnot(length(gene_names) == length(unique(gene_names)))
    stopifnot(length(tx_names) == length(unique(tx_names)))
    stopifnot(length(seqnames) <= n)
    stopifnot(length(uorfs_can_overlap_cds) == 1)
    stopifnot(uorfs_can_overlap_cds %in% seq(0, 2))
    stopifnot(min(cds_length) >= 9 & min(cds_length) > cds_exons & (min(cds_length) >= uorf_max_length + 3))
    stopifnot(identical(sort(stop_codons), sort(c("TAG", "TGA", "TAA"))))
    })
}

genome_output_test_controller <- function() {
  with(rlang::caller_env(),{
    # Sanity tests for output
    # Check CDS integrity:
    cds_AA <- translate(chromosome_seqs[cds_grl])
    cds_internal_stops <-  alphabetFrequency(cds_AA)[,"*"]
    if (any(cds_internal_stops != 1) & debug_on){
      warning("malformed CDS created, corner case, rerun to avoid this")
      browser()
      print(sum(cds_internal_stops != 1))
      print(cds_AA[cds_internal_stops != 1])
    }
    cds_starts <-  alphabetFrequency(translate(chromosome_seqs[startCodons(cds_grl, TRUE)]))[,"M"]
    stopifnot(all(cds_starts == 1))
  })
}

genome_exon_ranges_controller <- function() {
  with(rlang::caller_env(),{
    exon_cds_ranges <- cds_ranges
    mcols(exon_cds_ranges)[c("phase", "type")] <- DataFrame(phase = rep(as.integer(NA), length(exon_cds_ranges)), type = factor("exon"))
    cds_grl <- split(exon_cds_ranges, exon_cds_ranges$transcript_id)
    leader_ranges <- GRanges(seqnames, IRanges(flank_length + 1, width = leader_length), strand = strand)
    leaders_neg <- !strandBool(leader_ranges)
    ranges(leader_ranges[leaders_neg]) <-
      IRanges(end = (flank_length + leader_length + cds_genomic_length + trailer_length)[leaders_neg],
              width = width(leader_ranges[leaders_neg]))
    trailer_ranges <- GRanges(seqnames, IRanges(flank_length + leader_length + cds_genomic_length + 1, width = trailer_length), strand = strand)
    trailers_neg <- !strandBool(trailer_ranges)
    ranges(trailer_ranges[trailers_neg]) <-
      IRanges(start = flank_length[trailers_neg] + 1,
              width = width(trailer_ranges[trailers_neg]))

    mcols(leader_ranges) <-  gff_metadata(gene_names, transcript_id = tx_names, exon_number = "1")
    mcols(trailer_ranges) <- gff_metadata(gene_id = gene_names, transcript_id = tx_names, exon_number = "1")
    mcols(trailer_ranges) <- DataFrame(source = factor("toy_data"), type = factor("exon"), score = as.numeric(NA), phase = as.integer(NA),
                                       gene_id = gene_names, gene_version = 1,
                                       gene_biotype = "protein_coding", transcript_biotype = "protein_coding",
                                       transcript_id = tx_names, transcript_version = 1, exon_number = as.character(cds_exons + 2),
                                       exon_id = "NA")

    # Make Gene and Tx ranges
    gene_ranges <- GRanges(seqnames, IRanges(flank_length + 1,
                                             flank_length + leader_length + cds_genomic_length + trailer_length),
                           strand = strand, source = factor("toy_data"), type = factor("gene"), score = as.numeric(NA), phase = as.integer(NA),
                           gene_id = gene_names, gene_version = 1,
                           gene_biotype = "protein_coding", transcript_biotype = "protein_coding",
                           transcript_id = tx_names, transcript_version = 1, exon_number = "NA", exon_id = "NA")
    transcript_ranges <- gene_ranges
    transcript_ranges$type = "transcript"
  })
}

genome_exon_flanks_controller <- function() {
  with(rlang::caller_env(),{
    # GEN_CODE <- GENETIC_CODE
    # not_default_stop_code <- !(all(c("TAG", "TGA", "TAA") %in% stop_codons) & length(stop_codons) == 3)
    # if (not_default_stop_code) {
    #   GEN_CODE[GEN_CODE == "*" & !(names(GEN_CODE) %in% stop_codons)] <- "A"
    #   GEN_CODE[(names(GEN_CODE) %in% stop_codons)] <- "*"
    # }
    # Input argument re-definitions
    names(strand) <- tx_names
    cds_intron_length <- cds_intron_length
    total_introns <- cds_exons - 1

    flank_size <- flank_length #flank_length will be updated
    per_gene_length <- leader_length + cds_length + 6 +
      (cds_intron_length*total_introns) + trailer_length
    gene_and_intergenic_len <- flank_length + per_gene_length + flank_length
    flank_length <- rep_len(flank_length, n)
    if (length(seqnames) < n) { #Decide which genes are on same chromosomes
      gene_seqnames <- seq.int(length(seqnames)) #Make sure each chr get at least 1
      gene_seqnames <- sort(c(gene_seqnames,
                              sample(length(seqnames), n - length(seqnames), replace = TRUE)))
      seqnames <- seqnames[gene_seqnames]
      dt <- data.table::data.table(seqnames = gene_seqnames, strand, gene_and_intergenic_len)
      dt[, flank3_gene_end := cumsum(gene_and_intergenic_len), by = seqnames]
      dt[, flank5_gene_start := flank3_gene_end - gene_and_intergenic_len + 1]
      dt[, flank5_gene_end := flank5_gene_start + flank_size - 1]
      flank_length <- dt$flank5_gene_end
    }
    cds_starts <- ifelse(strand == "+",
                         flank_length + leader_length + 1,
                         flank_length + trailer_length + cds_length + 6
                         + (cds_intron_length*total_introns))
  })
}

make_transcriptome_sequences_controller <- function() {
  with(rlang::caller_env(),{
    message("-- Creating Transcriptome annotation")
    cds_string <- create_cds_seq(n, start_codons, stop_codons, cds_length)
    cds_list <- create_cds_ranges(cds_string, cds_starts, seqnames, strand,
                                  cds_exons, total_introns, cds_intron_length,
                                  gene_names, tx_names)
    cds_ranges <- cds_list[[1]]; cds_string_genomic <- cds_list[[2]]; rm(cds_list)
    cds_genomic_length <- width(cds_string_genomic) # Add start and stop codon to total length
    ## Make leader and trailer strings
    leader_string <- sample_dna_sequence(leader_length,
                                         c(A = 0.2, T = 0.2, G = 0.3, C = 0.3))
    trailer_string <- sample_dna_sequence(trailer_length,
                                          c(A = 0.3, T = 0.3, G = 0.2,C =  0.2))
    ## Make full gene string and flip - genes
    #cds_string_genomic[strand == "-"] <- reverseComplement(cds_string_genomic[strand == "-"])
    gene_string <- DNAStringSet(paste0(leader_string, cds_string_genomic, trailer_string))
    names(gene_string) <- tx_names
    gene_string[names(strand)[strand == "-"]] <- reverseComplement(gene_string[names(strand)[strand == "-"]])

  })
}

make_genome_sequences_controller <- function() {
  with(rlang::caller_env(),{
    message("-- Creating Genome annotation")
    # Make full chromosome string
    flank_size_all <- flank_size
    if (length(flank_size_all) != n) flank_size_all <- rep(flank_size_all, length.out = n)
    flank_5 <- sample_dna_sequence(flank_size_all)
    flank_3 <- sample_dna_sequence(flank_size_all)
    chromosome_seqs <- paste0(flank_5,  gene_string, flank_3)
    names(chromosome_seqs) <- seqnames
    chromosome_seqs <- DNAStringSet(chromosome_seqs)
    if (length(names(chromosome_seqs)) != length(unique(names(chromosome_seqs)))) {
      # Merge the genetic regions of same chromosomes
      chromosome_seqs <- split(chromosome_seqs, names(chromosome_seqs))
      chromosome_seqs <- DNAStringSet(lapply(a, function(x) unlist(x, use.names = FALSE)))
      chromosome_seqs <- chromosome_seqs[unique(seqnames)]
    }
  })
}

export_genome_output_controller <- function() {
  with(rlang::caller_env(),{
    # Save all annotation
    message("-- Saving genome and annotation")
    dir.create(out_dir, FALSE)
    #dir.create(file.path(out_dir, "ORFik_optimized"))
    out_prefix <- ORFik:::pasteDir(file.path(out_dir, genome_name))
    out_fasta <- paste0(out_prefix, ".fasta")
    out_gtf <- paste0(out_prefix, ".gtf")
    out_txdb <- paste0(out_gtf, ".db")
    writeXStringSet(chromosome_seqs, out_fasta)
    #rtracklayer::export.2bit(chromosome_seqs, file.path(out_dir, "Homo_sapiens_dummy.twoBit"))
    indexFa(out_fasta)
    rtracklayer::export.gff(con = out_gtf, gtf_ranges)
    simGenome <- c(genome = out_fasta, gtf = out_gtf)
    if (export_txdb) {
      txdb <- txdbmaker::makeTxDbFromGFF(out_gtf,
                                               chrominfo = seqinfo(findFa(out_fasta)),
                                               organism = "Homo sapiens")
      AnnotationDbi::saveDb(txdb, out_txdb)
      simGenome <- c(simGenome, txdb = out_txdb)
    }


    if (max_uorfs > 0) {
      saveRDS(uorf_ranges, file = file.path(out_dir, "true_uORFs.rds"))
      simGenome <- c(simGenome, uorfs = file.path(out_dir, "true_uORFs.rds"))
    }
    message("-- Done")
  })
}



make_gtf_from_ranges <- function(leader_ranges, cds_ranges, trailer_ranges,
                                 exon_cds_ranges, gene_ranges, transcript_ranges,
                                 tx_names, gene_names) {
  ## Reduce exons and and make correct exon ranks
  reduced_exons <- create_reduced_exons(leader_ranges, exon_cds_ranges, trailer_ranges,
                                        tx_names, gene_names)
  ## Sort gtf before save: order is gene, tx, exon, cds, exon, cds, next gene, tx, exon..
  cds_and_exons <- sort(c(reduced_exons, cds_ranges))
  gene_and_tx <- sort(c(gene_ranges, transcript_ranges))
  gtf_ranges <- sort(c(gene_and_tx, cds_and_exons))
  dt <- data.table::data.table(groupings = gtf_ranges$transcript_id, start = start(gtf_ranges),
                               end = end(gtf_ranges), type = gtf_ranges$type, index = seq.int(length(gtf_ranges)))
  data.table::setorder(dt, groupings, start, -end, type)
  gtf_ranges <- gtf_ranges[dt$index]
  names(gtf_ranges) <- NULL
  return(gtf_ranges)
}


sample_dna_sequence <- function(lengths, prob = c(A = 0.25, T = 0.25, G = 0.25, C = 0.25),
                                prob_5p_flank_seq = NULL) {
  seq_nchars <- unique(nchar(names(prob)))
  if (length(seq_nchars) != 1)
    stop("All input prob must have names of same length!")
  allowed_bases <- c("A", "T", "G", "C")
  unique_elements <- names(prob)
  if (seq_nchars > 1) {
    unique_elements <- unique(unlist(strsplit(names(prob), split = "")))
    lengths <- lengths / seq_nchars
    stopifnot(lengths == floor(lengths))
  }
  if (!any(allowed_bases %in% unique_elements))
    stop("prob must have at least 1 element named either: A, T, G, C !")

  if (sum(prob) != 1) stop("base probabilities must sum to 1!")
  if (!is.null(prob_5p_flank_seq)) {
    length_prob5 <- length(prob_5p_flank_seq)
    if (length_prob5 == 0)
      stop("When not NULL, prob_5p_flank_seq must have length > 1")
    if (is.character(prob_5p_flank_seq)) {
      prob_5p_flank_seq_temp <- rep(1/length_prob5, length_prob5)
      names(prob_5p_flank_seq_temp) <- prob_5p_flank_seq
      prob_5p_flank_seq <- prob_5p_flank_seq_temp
    }
    if (sum(prob_5p_flank_seq) != 1) stop("5' probabilities must sum to 1!")
    DNAStringSet(lapply(lengths, function(x)
      DNAString(paste0(sample(names(prob_5p_flank_seq), 1, prob = prob_5p_flank_seq),
                       paste(sample(names(prob), prob = prob, size = x, replace = TRUE),
                             sep = "", collapse = "")))))
  } else {
    DNAStringSet(lapply(lengths, function(x) DNAString(paste(
      sample(names(prob), prob = prob, size = x, replace = TRUE),
      sep = "", collapse = ""))))
  }
}


create_cds_seq <- function(n, start_codons, stop_codons, cds_length,
                           GEN_CODE = GENETIC_CODE,
                           GENETIC_CODE_prob = rep(1, length(GEN_CODE))/length(GEN_CODE)) {
  # Sample CDSs
  code_table <- GENETIC_CODE_prob
  names(code_table) <- names(GEN_CODE)
  code_table <- code_table[!(GEN_CODE %in% "*")]
  code_table <- (code_table / max(code_table) ) / length(code_table)
  cds_string <- sample_dna_sequence(cds_length, prob = code_table,
                                    prob_5p_flank_seq = start_codons)
  return(DNAStringSet(paste(cds_string, sample(stop_codons, length(cds_string), replace = TRUE), sep = "")))
}

gff_metadata <- function(gene_id, transcript_id, exon_number,
                         source = factor("toy_data"), type = factor("exon"),
                         score = as.numeric(NA), phase = as.integer(NA),
                         gene_version = 1, transcript_version = 1,
                         gene_biotype = "protein_coding",
                         transcript_biotype = "protein_coding",
                         exon_id = "NA") {
  DataFrame(source, type, score, phase, gene_id, gene_version,
            gene_biotype, transcript_biotype,
            transcript_id, transcript_version, exon_number,
            exon_id)
}

create_cds_ranges <- function(cds_string, cds_starts, seqnames, strand,
                              cds_exons, total_introns, cds_intron_length,
                              gene_names, tx_names) {
  # Make CDS Ranges with introns
  if (any(cds_exons > 1)) {
    intron_string <- sample_dna_sequence(cds_intron_length - 2,
                                         prob = c(A = 0.2, T = 0.2, G = 0.3, C = 0.3),
                                         prob_5p_flank_seq = "GT")
    cds_split <- floor(width(cds_string) / cds_exons)
    cds_string_genomic <- c("")
    current_start <-  1
    current_exon_start <- cds_starts
    cds_ranges <- GRanges()
    phase <- rep(0, length(cds_string))
    for (i in seq(total_introns)) {
      cds_string_genomic <- paste0(cds_string_genomic,
                                   DNAStringSet(cds_string, start = current_start, width = cds_split),
                                   intron_string)
      new_ranges <- IRanges(current_exon_start, width = cds_split)
      new_ranges[strand == "-"] <- IRanges(end = current_exon_start[strand == "-"],
                                           width = cds_split[strand == "-"])
      new_exons <- GRanges(seqnames, new_ranges, strand = strand,
                           source = factor("toy_data"), type = factor("CDS"), score = as.numeric(NA), phase = phase, # Update phase if needed
                           gene_id = gene_names, gene_version = 1,
                           gene_biotype = "protein_coding", transcript_biotype = "protein_coding",
                           transcript_id = tx_names, transcript_version = 1, exon_number = as.character(i), exon_id = "NA")
      cds_ranges <- c(cds_ranges, new_exons)
      current_start <- current_start + cds_split
      current_exon_start <- current_exon_start + ifelse(strand == "+", 1, -1)*(cds_split + width(intron_string))
      phase <- (phase + width(new_exons)) %% 3
    }
    last_exon <- DNAStringSet(cds_string, start = current_start, end = width(cds_string))
    cds_string_genomic <- DNAStringSet(paste0(cds_string_genomic, last_exon))
    new_ranges <- IRanges(current_exon_start, width = width(last_exon))
    new_ranges[strand == "-"] <- IRanges(end = current_exon_start[strand == "-"],
                                         width = width(last_exon)[strand == "-"])
    cds_ranges <- c(cds_ranges, GRanges(seqnames, new_ranges, strand = strand,
                                        source = factor("toy_data"), type = factor("CDS"), score = as.numeric(NA), phase = phase, # Update phase if needed
                                        gene_id = gene_names, gene_version = 1,
                                        gene_biotype = "protein_coding", transcript_biotype = "protein_coding",
                                        transcript_id = tx_names, transcript_version = 1, exon_number = as.character(i + 1), exon_id = "NA"))

  } else {
    cds_string_genomic <- cds_string
    new_ranges <- IRanges(cds_starts, width = cds_length + 6)
    new_ranges[strand == "-"] <- IRanges(end = cds_starts[strand == "-"],
                                         width = cds_length + 6)
    cds_ranges <- GRanges(seqnames, new_ranges, strand = strand,
                          source = factor("toy_data"), type = factor("CDS"), score = as.numeric(NA), phase = as.integer(0), # Update phase if needed
                          gene_id = gene_names, gene_version = 1,
                          gene_biotype = "protein_coding", transcript_biotype = "protein_coding",
                          transcript_id = tx_names, transcript_version = 1, exon_number = "1", exon_id = "NA")
  }
  return(list(cds_ranges, cds_string_genomic))
}

create_reduced_exons <- function(leader_ranges, exon_cds_ranges, trailer_ranges,
                                 tx_names, gene_names) {
  reduced_exons <- sort((c(leader_ranges, exon_cds_ranges, trailer_ranges)))
  reduced_exons <-reduceKeepAttr(split(reduced_exons, reduced_exons$transcript_id), keep.names = TRUE)
  reduced_exons <- reduced_exons[tx_names]
  dt <- data.table::data.table(groupings = groupings(reduced_exons), ones = 1L)
  dt <- dt[, exon_ranks := cumsum(ones), by = groupings]
  reduced_exons <- unlistGrl(reduced_exons)
  reduced_exons$exon_number <- dt$exon_ranks
  reduced_exons$gene_id <- gene_names[dt$groupings]
  reduced_exons$transcript_id <- tx_names[dt$groupings]
  return(reduced_exons)
}

create_uORFs <- function(leader_string, chromosome_seqs, n,
                         max_uorfs, start_codons, stop_codons,
                         leader_ranges, cds_ranges, cds_grl, seqnames, strand,
                         uorf_max_length, uorfs_can_overlap_cds, uorfs_can_overlap,
                         tx_names, debug_on) {
  message("-- Creating uORF annotation")
  names(leader_string) <- tx_names
  uorf_ranges <- create_uorf_ranges(leader_string, uorf_max_length, max_uorfs,
                                    uorfs_can_overlap_cds, uorfs_can_overlap)
  # Update leader,cds annotation with uORFs
  message("- Making valid uORFs sequences")
  uorf_groupings <- groupings(uorf_ranges)
  temp_all <- create_uorf_seq_template(uorf_ranges, uorf_groupings,
                                       seqnames, strand,
                                       start_codons, stop_codons)
  # Map uORFs to genomic coordinates
  uorf_ranges <- unlist(uorf_ranges, use.names = FALSE)
  names(uorf_ranges) <- uorf_groupings
  uorf_searchspace <- split(c(leader_ranges, cds_ranges), c(leader_ranges$transcript_id, cds_ranges$transcript_id))
  uorf_searchspace <- sortPerGroup(reduce(uorf_searchspace))[leader_ranges$transcript_id]
  uorf_ranges <- pmapFromTranscriptF(uorf_ranges, uorf_searchspace, removeEmpty = TRUE)
  mcols(uorf_ranges) <- NULL
  seqinfo(uorf_ranges) <- seqinfo(chromosome_seqs)

  seqnames_uorfs <- seqnamesPerGroup(uorf_ranges, FALSE)

  # Tile uORFs for overlap tests
  uorf_ranges_index <- uorf_ranges
  names(uorf_ranges_index) <- seq.int(length(uorf_ranges_index))
  tile <- start(tile1(uorf_ranges_index, sort.on.return = FALSE, matchNaming = FALSE))
  # Avoid CDS start codon overlapping parts
  startCodons <- sortPerGroup(reduce(split(cds_ranges, cds_ranges$transcript_id)), quick.rev = TRUE)
  startCodons <- ORFik::startCodons(startCodons[unique(cds_ranges$transcript_id)], TRUE)
  startCodons <- start(tile1(startCodons, sort.on.return = FALSE, matchNaming = FALSE))
  overlapping_list <- !(tile %in% startCodons[uorf_groupings])
  tile <- tile[overlapping_list]
  for (i in seq_along(uorf_ranges)) { # Insert uORF sequence except CDS start codon
    chr <- seqnames_uorfs[i]
    chromosome_seqs[[chr]] <- replaceLetterAt(chromosome_seqs[[chr]],
                                              at = tile[[i]],
                                              letter = temp_all[[i]][overlapping_list[[i]]])
  }
  # Fix 4th base of CDS leading to internal stop codon in uORF
  freq_stop <- alphabetFrequency(translate(chromosome_seqs[uorf_ranges]))[,"*"]
  start_pluss_1_wrong <- which(freq_stop > 1)
  for (i in start_pluss_1_wrong) {
    chr <- seqnames_uorfs[i]
    ir <- ranges(uorf_ranges[i])
    names(ir) <- chr
    internal_stops <- unlist((start(Biostrings::vmatchPattern("*", translate(chromosome_seqs[uorf_ranges[i]]))))*3) -2
    internal_stops <- internal_stops[-length(internal_stops)]
    if (length(names(internal_stops)) != length(seq(length(internal_stops)))) next
    names(internal_stops) <- seq(length(internal_stops))
    pick_base_to_change <- sort(c(internal_stops, internal_stops+1, internal_stops+2))
    pick_base_to_change <- pick_base_to_change[!(pick_base_to_change %in%
                                                   which(if(strandBool(uorf_ranges[i])) {!overlapping_list[[i]]} else {rev(!overlapping_list[[i]])}))]
    pick_base_to_change <- pick_base_to_change[!duplicated(names(pick_base_to_change))]
    chromosome_seqs[ir][[1]] <- if (strandBool(uorf_ranges[i])) {
      replaceLetterAt(chromosome_seqs[uorf_ranges[i]][[1]],
                      pick_base_to_change,
                      DNAString(paste(sample(c("C", "G"), length(pick_base_to_change),
                                             replace = TRUE), collapse = "")))
    } else {
      reverseComplement(replaceLetterAt(chromosome_seqs[uorf_ranges[i]][[1]],
                                        pick_base_to_change,
                                        DNAString(paste(sample(c("C", "G"), length(pick_base_to_change),
                                                               replace = TRUE), collapse = ""))))
    }
  }
  # Fix CDSs that now got new internal stop codons
  new_cds_string <- chromosome_seqs[cds_grl]
  internal_inframe_stops <- 3 + (start(Biostrings::vmatchPattern("*", heads(translate(new_cds_string), -1))) - 1)*3
  if(!all(lengths(internal_inframe_stops) == 0)) {
    stops_matrix <- as.matrix(IntegerList(lapply(seq(n), function(x) seq.int(width(new_cds_string[x])))) %in% internal_inframe_stops)
    temp <- Biostrings::replaceLetterAt(new_cds_string, stops_matrix,
                                        letter = unlist(DNAStringSetList(lapply(rowSums(stops_matrix), function(x) paste(rep("C", x), collapse = "")))))
    a <- cds_grl; names(a) <- NULL
    b <- temp; b[!strandBool(a)] <- reverseComplement(b[!strandBool(a)])
    matching <- chmatch(names(chromosome_seqs), seqnamesPerGroup(a, FALSE))
    chromosome_seqs[ranges(sort(a[matching]))] <- as(unlist(b[matching]), "DNAStringSet")
  }

  # Fix uORF Starts that were overwritten
  uorfs_without_start_sc <- startCodons(uorf_ranges, TRUE)
  total_start <- alphabetFrequency(translate(chromosome_seqs[uorfs_without_start_sc]))[,"M"]
  if (any(total_start == 0)) {
    uorfs_without_start <- uorf_ranges[total_start == 0]
    seqnames_uorfs_ws <- seqnamesPerGroup(uorfs_without_start, FALSE)
    for (i in seq_along(uorfs_without_start)) {
      ir <- ranges(uorfs_without_start[i])
      names(ir) <- seqnames_uorfs_ws[i]
      chromosome_seqs[ir][[1]] <- if (strandBool(uorfs_without_start[i])) {
        replaceLetterAt(chromosome_seqs[uorfs_without_start[i]][[1]],
                        c(1,2,3),
                        DNAString(sample(start_codons, 1)))
      } else {
        reverseComplement(replaceLetterAt(chromosome_seqs[uorfs_without_start[i]][[1]],
                                          c(1,2,3),
                                          DNAString(sample(start_codons, 1))))
      }
    }
  }
  # Fix uORF Stops that were overwritten
  uorfs_without_stop_sc <- stopCodons(uorf_ranges, TRUE)
  total_stops <- alphabetFrequency(translate(chromosome_seqs[uorfs_without_stop_sc]))[,"*"]
  if (any(total_stops == 0)) {
    uorfs_without_stop_sc <- uorfs_without_stop_sc[total_stops == 0]
    seqnames_uorfs_ws <- seqnamesPerGroup(uorfs_without_stop_sc, FALSE)
    for (i in seq_along(uorfs_without_stop_sc)) {
      ir <- ranges(uorfs_without_stop_sc[i])
      names(ir) <- seqnames_uorfs_ws[i]
      chromosome_seqs[ir][[1]] <- if (strandBool(uorfs_without_stop_sc[i])) {
        replaceLetterAt(chromosome_seqs[uorfs_without_stop_sc[i]][[1]],
                        c(1,2,3),
                        DNAString(sample(stop_codons, 1)))
      } else {
        reverseComplement(replaceLetterAt(chromosome_seqs[uorfs_without_stop_sc[i]][[1]],
                                          c(1,2,3),
                                          DNAString(sample(stop_codons, 1))))
      }
    }
  }
  # Fix uORFs that now got new internal stop codons
  print("Fixing uORF internal stop --")
  new_uorf_string <- chromosome_seqs[uorf_ranges]
  internal_inframe_stops <- 3 + (start(Biostrings::vmatchPattern("*", heads(translate(new_uorf_string), -1))) - 1)*3
  if(!all(lengths(internal_inframe_stops) == 0)) {
    hits <- lengths(internal_inframe_stops) > 0
    temp <- rep(ifelse(strandBool(uorf_ranges[hits]), "C", "G"), lengths(internal_inframe_stops)[hits])
    temp <- DNAString(paste(strsplit(temp, split = ""), collapse = ""))
    a <- IRangesList(internal_inframe_stops[hits]); names(a) <- which(hits)
    a <- pmapFromTranscriptF(a, uorf_ranges); names(a) <- NULL
    chromosome_seqs[seqnamesPerGroup(a, FALSE)][ranges(a)] <- as(temp, "DNAStringSet")
  }
  print("uORF internal stop done")
  # uORF sanity tests
  total_stops <- alphabetFrequency(translate(chromosome_seqs[uorf_ranges]))[,"*"]
  start_debug <- FALSE
  if (any(total_stops > 1) & debug_on) {
    warning("Some uORFs have internal stop codons that could not be fixed!")
    message("Starting browser debug mode, rerun if you do not care about the details")
    start_debug <- TRUE
  } else if (any(total_stops == 0) & debug_on) {
    warning("Some uORFs have no stop codon")
    message("Starting browser debug mode, rerun if you do not care about the details")
    start_debug <- TRUE
  }
  if (start_debug) {
    all <- which(total_stops != 1)
    i <- all[1]
    browser()
    print(table(total_stops))
    print(translate(chromosome_seqs[uorf_ranges[i]]))
    print(chromosome_seqs[uorf_ranges[i]])
    print(table(as.integer(table(uorf_groupings))))
    print(strandPerGroup(uorf_ranges[i], FALSE))
  }
  return(list(uorf_ranges, chromosome_seqs))
}

create_uorf_ranges <- function(leader_string, uorf_max_length, max_uorfs,
                               uorfs_can_overlap_cds, uorfs_can_overlap) {
  message("- Making valid uORFs ranges")
  removed_start_area <-
  uorf_ranges <- IRangesList(lapply(nchar(leader_string), function(size) {
    all_uorfs <- IRanges()
    # Not allow to start overlapping CDS ATG, or inframe of CDS
    if (uorfs_can_overlap_cds) {
      allowed_cds_spots <- seq.int(size+7, size + uorf_max_length - 4)
      allowed_cds_spots <- allowed_cds_spots[seq.int(length(allowed_cds_spots)) %% 3 != 0]
    } else allowed_cds_spots <- NULL

    possible_starts <- seq.int(size-3)
    possible_ends <- c(allowed_cds_spots)
    if (uorfs_can_overlap_cds != 2) possible_ends <- c(possible_starts[-seq.int(5)], possible_ends)
    for (i in seq(max_uorfs)) {
      if (length(possible_starts) == 0 | length(possible_ends) == 0) break;
      end <- ifelse(length(possible_ends) == 1, possible_ends, sample(possible_ends, size = 1))
      pstart_all <- seq.int(end-2, max(min(possible_starts), end - uorf_max_length), by = -3)[-1]
      pstart <- pstart_all[pstart_all %in% possible_starts]
      if (i > 1) {
        uorf_upwards <- end > end(all_uorfs)
        if (uorfs_can_overlap) {
          uorf_upwards <- uorf_upwards & ((end - end(all_uorfs)) %% 3 == 0)
        }
        if (any(uorf_upwards)) {
          pstart <- pstart[pstart > max(end(all_uorfs)[uorf_upwards])]
        }
      }
      if (length(pstart) == 0) next
      start <- ifelse(length(pstart) == 1, pstart, sample(pstart, size = 1))
      #overlap_same_frame_cds <- end > size & (((end-(size+1)) %% 3) == 2)
      uorf <- sort(IRanges(start, end))
      names(uorf) <- paste0("uORF", i)
      all_uorfs <- c(all_uorfs, uorf)
      invalid_start <- if (uorfs_can_overlap) { # Start / end restriction
        c((start-2):(start+2), ((end-5):(end)))
      } else seq.int(start-6, end)
      invalid_end <- if (uorfs_can_overlap) {
        c(pstart_all[pstart_all >= start] + 2, # Invalid in frame as stopcodon inside
          invalid_start[-length(invalid_start)], # Invalid around start and stop, allow same stop
          seq(end+1, end+3)) # Invalid after end
      } else c(invalid_start[-c(1,2,3,4,5)], seq(end+1, end+5))
      possible_ends <- possible_ends[!(possible_ends %in% invalid_end)]
      possible_starts <- possible_starts[!(possible_starts %in% invalid_start)]
      possible_ends <- possible_ends[(possible_ends > min(possible_starts) + 5) &
                                       (possible_ends < max(possible_starts) + uorf_max_length -1)]
    }
    #i; all_uorfs; cat("start:"); possible_starts; cat("stop:");possible_ends; end
    all_uorfs
  }))
  return(uorf_ranges)
}

create_uorf_seq_template <- function(uorf_ranges, uorf_groupings,
                                     seqnames, strand,
                                     start_codons, stop_codons) {
  uorf_widths <- width(unlist(uorf_ranges))
  temp_all <- DNAStringSet(lapply(seq_along(uorf_groupings), function(x)
    DNAString(paste0(sample(start_codons, 1),
                     paste(sample(c("A", "T", "G", "C"), size = uorf_widths[x] - 6,
                                  replace = TRUE),
                           sep = "", collapse = "")))
  ))
  internal_inframe_stops <- 1 + (start(Biostrings::vmatchPattern("*", translate(temp_all))) - 1)*3
  for (i in seq_along(uorf_groupings)) {
    stops <- internal_inframe_stops[[i]]
    if (length(stops) > 0) {
      temp_all[[i]] <- Biostrings::replaceLetterAt(temp_all[[i]], stops,
                                                   letter = rep("C", length(stops)))
    }
  }
  temp_all <- DNAStringSet(paste0(temp_all, sample(stop_codons, length(temp_all), replace = TRUE)))
  temp_all[strand[uorf_groupings] == "-"] <- reverseComplement(temp_all[strand[uorf_groupings] == "-"])
  names(temp_all) <- seqnames[uorf_groupings]
  return(temp_all)
}

distribute_reads_to_uORFs <- function(region_counts, assay, uorf_ranges, uorf_prop_mode,
                                      uorf_prop_within_gene) {
  if (is.null(names(region_counts))) names(region_counts) <- rownames(assay)
  uorf_tab <- table(txNames(uorf_ranges))

  prob <- if (uorf_prop_mode == "character") {
    if (uorf_prop_within_gene == "uniform") {
      a <- rep(1, length(uorf_ranges))
      names(a) <- txNames(uorf_ranges)
    } else {
      stop("Not implemented yet")
      lapply(names(uorf_tab), function(x) rmultinom(1, region_counts[x], rep(1, uorf_tab[x])))
    }
    a
  } else uorf_prop_within_gene
  a <-  lapply(names(uorf_tab), function(x) rmultinom(1, region_counts[x], prob[names(prob) == x]))
  a <- unlist(a); names(a) <- txNames(uorf_ranges)
  return(a)
}

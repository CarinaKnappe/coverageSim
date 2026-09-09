# IDEA:
# Sample regions by beta distribution, with given shape parameters.
# Sample then positional coverage by either codon, or location bias. Urn sampling.
# Maybe then also sample noise from beta again

#' Simulate NGS coverage on from count tables on a genome.
#'
#' Simulated read coordinates from count tables given a genome.
#' For each read type bias the way reads are distributed by their properties.
#' Like sampling type, variance, skew, region bias, sequence bias etc.
#' @param simGenome named character vector, the genome to simulate data on.
#' @param count_table a summarized experiment of gene counts with rownames being transcript names.
#' Must contain at least 2 out of 5 named assays: gene, leaders, cds, trailer, uorf
#' @param out_dir directory to save simulated read files, name of individual files decided by naming columns
#' @param exp_name name of ORFik experiment to be created for these read files and annotation
#' @param exp_save_dir path, default ORFik::config()["exp"]. Set to tempdir()
#' or other location to not clutter your real experiments.
#' @param transcripts character vector, default \code{names(count_table)}, default uses all transcripts,
#' alternativly a subset of genes to assign reads to.
#' @param ideal_coverage list, default: list(leader =  list(CAGE = quote(c(1, rep.int(0, x -1)))),
#' cds =          list(RFP = quote(rep.int(c(1, 0, 0), length.out = x))),
#' trailer =      list(PAS = quote(c(rep.int(0, x -1), 1))),
#' uorf =         list(RFP = quote(rep.int(c(1, 0, 0), length.out = x))))
#' @param rnase_bias = list(leader =  list(RFP = 1, RNA = 1, CAGE = quote(1/seq.int(x)), PAS = quote(1/seq.int(x, 1))),
#'                  cds =     list(RFP = c(5, 2, 1), RNA = 1, CAGE = quote(1/seq.int(x)), PAS = quote(1/seq.int(x, 1))),
#'                  trailer = list(RFP = 1, RNA = 1, CAGE = quote(1/seq.int(x)), PAS = quote(1/seq.int(x, 1))),
#'                  uorf = list(RFP = c(5, 2, 1), RNA = 1, CAGE = quote(1/seq.int(x)), PAS = quote(1/seq.int(x, 1))))
#' @param auto_correlation = list(leader = list(RFP = 1, RNA = 1, CAGE = 1, PAS = 1),
#'    cds = list(RFP = quote(abs(cos(seq(0, pi*(x/30), pi/30))) + 0.1), RNA = 1, CAGE = 1, PAS = 1),
#'    trailer = list(RFP = 1, RNA = 1, CAGE = 1, PAS = 1),
#'    uorf = list(RFP = quote(abs(cos(seq(0, pi*(x/30), pi/30))) + 0.1), RNA = 1, CAGE = 1, PAS = 1))
#'    - Alternative: abs(cos(x)) + abs(cos(x*b[i])),
#'     where b is an internal value sampled per gene in interval 3 to 15.
#' @param read_lengths_per a list, default: list(RFP = 27:29, RNA = 100, CAGE = 1, PAS = 1)
#' @param sampling a list, default: list(leader = list(RFP = "MN", RNA = "MN"),
#' cds =    list(RFP = "DMN"),
#' trailer =list(RFP = "MN"),
#' uorf =   list(RFP = "DMN"))
#' @param libClasses = list(RFP = "GRanges", RNA = "GAlignment", CAGE = "GRanges", PAS = "GRanges")
#' @param libFormats The output formats, a named list of characters:
#'  list(RFP = "ofst", RNA = "ofst", CAGE = "ofst", PAS = "ofst"). Alternatives:
#'  "sam" and "bam".
#' @param true_uorf_ranges = "AUTO". Load from uorf string in 'simGenome'.
#' @param seq_bias the sequence bias used for simulation, default:
#' load_seq_bias(), the Ribosome profiling estimates from Amino acid dwell times from
#' P-site shifted reads from library "R2" from the coverageSim paper. Specifically it is:
#' a data.table with correctly named columns:
#'  "seq", "alpha"\cr
#'  the column seq (Amino acid, 1 letter per)
#'  or by codons (3 letter per), and a mean column and dispersion, with relative usage.
#'  In addition iMethionin
#'  (start codon) can be differentiated by setting seq = "#" and stop codon
#'  as seq = "*". Gives ability to scale these two important regions seperatly
#'  from other codons.
#'  If a \code{variable} column is present, it must identify exactly one profile.
#'  Tables returned by \code{load_seq_bias(bias = "all")} must be subset first.
#' @param uorf_prop_within_gene = "uniform". How should counts for uORF be
#' distributed. Alternatives: 'length' (length biased), or user specified
#' list of numeric values for all ORFs.
#' @param validate logical, TRUE, check that experiment was correctly made.
#' Set to false if you know it won't and you will fix it later. This happens
#' if you want to put a new sample with different parameters,
#' in a folder with other libraries you made
#' earlier, so then only validate when you create the ORFik experiment after the
#' last sample is made.
#' @param fragment_mode RFP representation. `"simulated_rpf"` (default) creates
#' complete simulated transcript-derived RPFs; `"legacy_point"` preserves the
#' historical behavior where the simulated signal position is the read start.
#' The former name `"physical"` remains as a deprecated compatibility alias.
#' @param fragment_geometry Named list controlling simulated RPFs. Supports
#' `source` (`"default"`, `"user"`, or `"learned"`), `site_reference`
#' (`"p_site"` or `"a_site"`), a joint `distribution` with fragment_length,
#' site_offset, and probability columns, and `boundary_action`. The default
#' `"renormalize"` conditions probabilities on valid boundary geometry and
#' preserves counts; `"error"` rejects impossible geometry. End-bias source
#' fields are reserved and currently must be `list(source = "none")`.
#' End-bias sources may also be `"user"` or `"learned"` with a table of
#' `kmer`, `weight`, and optionally `fragment_length`; use
#' \code{make_synthetic_end_bias()} to create a simple explicit profile.
#' @param ground_truth FALSE, TRUE, or a directory path. In simulated-RPF mode, TRUE
#' writes one compressed fragment truth table next to each simulated library.
#' @param debug_coverage logical, default FALSE. If TRUE, debug steps of coverage calculation,
#' for parameter errors or just to understand how it all works.
#' @return an \code{\link[ORFik]{experiment}}
#' @import ORFik data.table GenomicRanges
#' @export
#' @examples
#' ## Simple example
#' # 6 genes on 6 chromosomes (No active uORFs)
#' simGenome6 <- simGenome(n = 6, max_uorfs = 0)
#' # Simulate Ribo-seq only
#' gene_count_table <-simCountTables(loadRegion(simGenome6["txdb"], "cds"),
#'  libtypes = "RFP", print_statistics = FALSE)
#' region_count_table <- simCountTablesRegions(gene_count_table,
#'                  regionsToSample = c("leader", "cds", "trailer"))
#' simNGScoverage(simGenome6, region_count_table)
simNGScoverage <- function(simGenome,
                           count_table = simCountTablesRegions(
                             simCountTables(loadRegion(simGenome["txdb"], "cds"))),
                           out_dir = file.path(dirname(simGenome["genome"]), ""),
                           exp_name = "simulated_data",
                           exp_save_dir = ORFik::config()["exp"],
                           transcripts = names(count_table),
                           ideal_coverage = list(leader =  list(CAGE = quote(c(1, rep.int(0, x -1)))),
                                            cds =          list(RFP = quote(rep(c(1, 0, 0), length.out = x))),
                                            trailer =      list(PAS = quote(c(rep.int(0, x -1), 1))),
                                            uorf =         list(RFP = quote(rep(c(1, 0, 0), length.out = x)))),
                           rnase_bias = list(RFP = c(0.5,2,1,10,2,1,0.5), RNA = rnase_models(),
                                             CAGE = rnase_models(), PAS = rnase_models()),
                           auto_correlation = list(cds = list(RFP = shapes(9)),
                                                   uorf =list(RFP = shapes(9))),
                           read_lengths_per = list(RFP = 27:29, RNA = 100,
                                                   CAGE = 1, PAS = 1),
                           sampling = list(leader = list(RFP = "MN", RNA = "MN"),
                                           cds =    list(RFP = "DMN"),
                                           trailer =list(RFP = "MN"),
                                           uorf =   list(RFP = "DMN")),
                           libClasses = list(RFP = "GRanges", RNA = "GAlignment",
                                             CAGE = "GRanges", PAS = "GRanges"),
                           libFormats = list(RFP = "ofst", RNA = "ofst",
                                             CAGE = "ofst", PAS = "ofst"),
                           seq_bias = load_seq_bias(),
                           true_uorf_ranges = "AUTO", uorf_prop_within_gene = "uniform",
                           validate = TRUE,
                           fragment_mode = c("simulated_rpf", "legacy_point", "physical"),
                           fragment_geometry = list(
                             source = "default", site_reference = "a_site",
                             distribution = NULL, boundary_action = "renormalize",
                             five_prime_bias = list(source = "none"),
                             three_prime_bias = list(source = "none")
                           ),
                           ground_truth = FALSE,
                           debug_coverage = FALSE) {
  fragment_mode <- match.arg(fragment_mode)
  validate_sequence_profile(seq_bias)
  if (fragment_mode == "physical") {
    warning("fragment_mode = 'physical' is deprecated; use 'simulated_rpf'")
    fragment_mode <- "simulated_rpf"
  }
  fragment_geometry <- normalize_fragment_geometry(fragment_geometry)
  input_validation_controller()

  # Load annotation
  txdb <- loadTxdb(simGenome["txdb"])
  if (fragment_mode == "simulated_rpf") {
    mrna_ranges <- loadRegion(txdb, "mrna", names.keep = transcripts)
    models <- transcript_models(mrna_ranges, simGenome["genome"])
  }
  loadRegions(txdb, parts = regionsToSample[!(regionsToSample %in% "uorf")],
              envir = environment(), extension = "_ranges", names.keep = transcripts)
  if ("cds" %in% regionsToSample) {
    if (any((widthPerGroup(cds_ranges, FALSE) %% 3) != 0)) {
      warning("Detected CDS ranges that ends on incomplete codon (is not %% 3 == 0 in length")
    }
  }
  # Create sequence table and sequence bias, initiate coverage table
  sequence_table_controller()

  # Set up count tables
  assay <- assay(count_table); assay <- assay - assay
  for (i in seq_along(assayNames(count_table))[-1]) {
    assay <- assay + assay(count_table, i)
  }
  assay_by_chromosome <- assay_by_chromo(assay,
                seqnamesPerGroup(cds_ranges, keep.names = FALSE))

  colnames <- colnames(assay)
  libtypes <- as.character(colData(count_table)$libtype)
  files <- c()
  message("- Sample coverage")
  # For each sample create pdf for each (ORF/region) and sample
  for (s in seq_along(colnames)) {
    assay_column <- colnames[s]
    message("-- ", assay_column)

    libClass <- unlist(libClasses[libtypes[s]], use.names = FALSE)
    # Step 4 (NT coverage distribution)
    dt_final <- nt_coverage_all_regions(count_table[, s], libClass, ideal_coverage,
                                        rnase_bias, auto_correlation, read_lengths_per,
                                        uorf_ranges, uorf_prop_mode,
                                        uorf_prop_within_gene, sampling,
                                        debug_coverage, env = environment())
    # Verify all reads have been distributed correctly
    expected_counts <- data.table::data.table(
      seqnames = assay_by_chromosome$seqnamesPer,
      expected = assay_by_chromosome[, assay_column, with = FALSE][[1]]
    )
    actual_counts <- dt_final[, .(actual = sum(score)), by = "seqnames"]
    count_check <- merge(expected_counts, actual_counts, by = "seqnames", all = TRUE)
    count_check[is.na(count_check)] <- 0
    stopifnot(all(count_check$expected == count_check$actual))
    dt_final <- dt_final[score > 0,]
    simulated_rpf <- fragment_mode == "simulated_rpf" && libtypes[s] == "RFP"
    if (simulated_rpf) {
      fragment_table <- make_simulated_rpf_fragments(
        dt_final,
        transcript_models = models,
        fragment_lengths = unlist(read_lengths_per[libtypes[s]], use.names = FALSE),
        fragment_geometry = fragment_geometry
      )
      stopifnot(sum(fragment_table$score) == sum(dt_final$score))
      gr_final <- simulated_rpf_alignments(
        fragment_table,
        seqinfo = GenomeInfoDb::seqinfo(mrna_ranges)
      )
    } else if (libClass == "GRanges") {
      gr_final <- makeGRangesFromDataFrame(dt_final, keep.extra.columns = TRUE)
    } else {
      gr_final <- ORFik:::getGAlignments(dt_final)
    }
    format <- resolve_export_format(libFormats, libtypes[s])
    file_base <- file.path(out_dir, assay_column)
    if (simulated_rpf) {
      write_fragment_ground_truth(fragment_table, file_base, ground_truth)
    } else if (!identical(ground_truth, FALSE) && !is.null(ground_truth)) {
      warning("Ground truth is currently written only for simulated RPF libraries")
    }
    written_files <- write_simulated_library(
      gr_final,
      file_base = file_base,
      format = format,
      seqinfo = GenomeInfoDb::seqinfo(gr_final)
    )
    files <- c(files, unname(written_files["default"]))
  }
  replicates <- as.character(colData(count_table)$replicate)
  conditions <- as.character(colData(count_table)$condition)
  create.experiment(out_dir, exp_name, txdb = simGenome["txdb"],
                    fa = simGenome["genome"], saveDir = exp_save_dir,
                    organism = "Homo sapiens", author = "Simulated by ORFikSim",
                    libtype = libtypes,
                    condition = conditions,
                    rep = replicates, #[order(colnames(count_table))]
                    files = files)
  message("Saved ORFik experiment with name:")
  message(exp_name)
  return(read.experiment(exp_name, in.dir = exp_save_dir, validate = validate))
}

rnase_models <- function(radius = c(0.5,1,2,10,2,1,0.5),
                         type = "RFP") {
  if (type == "exp") {
    window <- c()
  }
  return(1)
}

nt_coverage_all_regions <- function(count_table_regions, libClass,
                                    ideal_coverage,
                                    rnase_bias, auto_correlation,
                                    read_lengths_per, uorf_ranges,
                                    uorf_prop_mode, uorf_prop_within_gene,
                                    sampling,
                                    debug_coverage, env) {
  regionsToSample <- assayNames(count_table_regions)[-1]
  data.table::rbindlist(lapply(regionsToSample, function(region) {
    assay <- assay(count_table_regions, region)
    region_counts <- as.vector(assay)
    if (sum(region_counts) > 0) {
      lengths <- get(paste0("lengths_", region), envir = env)
      libtype <- as.character(count_table_regions$libtype)

      if (region == "uorf") {
        region_counts <-  distribute_reads_to_uORFs(region_counts, assay,
                                                    uorf_ranges, uorf_prop_mode,
                                                    uorf_prop_within_gene)
      }
      region_ranges <- get(paste0(region, "_ranges"), mode = "S4", envir = env)
      dt_region <- copy(get(paste0("dt_", region), envir = env))
      ideal_cov <- get_value(ideal_coverage, region, libtype) # skew_cov
      if (is.null(ideal_cov)) ideal_cov <- quote(rep(1, x))
      auto_cor <- get_value(auto_correlation, region, libtype) # bump_shape
      read_lengths <- unlist(read_lengths_per[libtype], use.names = FALSE)
      sampling_mode <- get_value(sampling, region, libtype)
      if (is.null(sampling_mode)) sampling_mode <- "MN"
      if (debug_coverage) browser()
      if (sampling_mode == "DMN") {
        region_length_matrix <- get(paste0("region_length_matrix", region),
                                    envir = env)
        alpha_matrix <- get(paste0("seq_alpha_list_", region), envir = env)
        # Calculate bias of coverage
        res <- sim_sequence_bias(ideal_cov, lengths,
                                 alpha_matrix, auto_cor,
                                 rnase_bias[[libtype]])
        res_lengths <- lengths(res)
        res_matrix <- pack_alpha_rows(res, region_length_matrix)
        #i <- 3; 57- sum(alpha_mat_3[i,] == 1e-24); lengths[i]
        n_genes <- length(res)
        sample <- extraDistr::rdirmnom(n = n_genes, size = region_counts,
                                       alpha = res_matrix)
        sample <- flatten_sample_rows(sample, res_lengths)
      } else { #MN
        sample <- lapply(seq_along(region_ranges),
                          function(y, fun, x = lengths[y])
                            a <- as.vector(rmultinom(1, region_counts[y], eval(fun))),
                            fun = ideal_cov)
        sample <- unlist(sample, use.names = FALSE)
      }

      dt_region <- dt_region[, !(colnames(dt_region) %in% c("AA", "genes", "position", "proportion")), with = FALSE]
      dt_region[, score := sample]
      #dt_region[, score := unlist(sample, use.names = FALSE)]
      if (libClass == "GRanges") {
        if (length(read_lengths) == 1) {
          dt_region[, size := rep(read_lengths, .N)]
        } else dt_region[, size := sample(read_lengths, size = .N, replace = TRUE)]
      } else {
        if (length(read_lengths) == 1) {
          dt_region[, cigar := paste0(sample(rep(read_lengths, .N), size = .N, replace = TRUE), "M")]
        } else dt_region[, cigar := paste0(sample(read_lengths, size = .N, replace = TRUE), "M")]
      }
    } else dt_region <- data.table::data.table()
    return(dt_region)
  }))
}

# Some ideas for CIGAR I did not create
# mrna <- loadRegion(df, "mrna")
# overlaps <- findOverlaps(RFP, mrna)
# overlaps <- overlaps[!duplicated(from(overlaps))]
# x_ir <- IRanges(start=start(RFP), width = readWidths(RFP), names = seq(length(RFP)))

#' Convert GAlignment to sam/bam
#' @param x A GRanges or GAlignment object
#' @param path character path to save sam/bam file
#' @param seqinfo a Seqinfo object, default GenomeInfoDb::seqinfo(x)
#' @param MAPQ numeric, full read quality score, default 30.
#' @param sequences character, default "*", empty, this saves space if not needed.
#' @param per_base_quality character, per base quality score, default "*"
#' @param RNEXT character, default "*", see bam format manual if relevant
#' @param PNEXT numeric, default 0, see bam format manual if relevant
#' @param TLEN numeric, default 0, see bam format information if relevant
#' @param FLAGS numeric flags, default is vector of ifelse(strandBool(x), 0, 16),
#' giving +/- strand flag.
#' @param make_bam logical, default TRUE
#' @return a named character vector of path to SAM and BAM (if bam was created)
samFromGAlignment <- function(x, path, seqinfo = GenomeInfoDb::seqinfo(x),
                              MAPQ = 30, sequences = "*", per_base_quality = "*",
                              RNEXT = "*", PNEXT = 0, TLEN = 0, FLAGS = ifelse(strandBool(x), 0, 16),
                              make_bam = TRUE) {
  stopifnot(is(seqinfo, "Seqinfo"))
  if (is(x, "GRanges")) x <- GAlignments(seqnames = seqnames(x), pos = start(x),
                                         cigar = paste0(readWidths(x),"M"), strand = strand(x),
                                         score = mcols(x)$score)
  if (identical(sequences, "*") && !is.null(S4Vectors::mcols(x)$reference_sequence)) {
    # Reference/genome-plus-strand-oriented sequence (see fragment_geometry.R),
    # already reverse-complemented for minus-strand alignments where needed.
    # Preferred over `sequence` (transcript-oriented) whenever present.
    sequences <- as.character(S4Vectors::mcols(x)$reference_sequence)
  } else if (identical(sequences, "*") && !is.null(S4Vectors::mcols(x)$sequence)) {
    sequences <- as.character(S4Vectors::mcols(x)$sequence)
  }
  chr_header <- paste("@SQ", paste0("SN:", seqnames(seqinfo)),
                      paste0("LN:",seqlengths(seqinfo)), sep = "\t")

  reads <- paste(FLAGS, seqnames(x), start(x), MAPQ, cigar(x),
                 RNEXT, PNEXT, TLEN, sequences, per_base_quality, "NH:i:1", sep = "\t")
  if (!is.null(mcols(x)$score)) reads <- reads[rep.int(seq(length(reads)), mcols(x)$score)]
  reads <- paste(paste0("read", seq(length(reads))), reads, sep = "\t")

  writeLines(c(chr_header, reads), path)
  res <- c(SAM = path)
  if (make_bam) {
    res <- c(res, BAM = Rsamtools::asBam(path))
  }
  message("Done")
  return(res)
}

normalize_lib_formats <- function(libFormats) {
  if (is.character(libFormats)) {
    return(as.list(libFormats))
  }
  libFormats
}

resolve_export_format <- function(libFormats, libtype) {
  libFormats <- normalize_lib_formats(libFormats)
  format <- unlist(libFormats[libtype], use.names = FALSE)
  if (!length(format)) {
    stop("Missing libFormats entry for libtype: ", libtype)
  }
  if (length(format) != 1) {
    stop("Exactly one output format must be provided per libtype")
  }
  format
}

validate_lib_formats <- function(libFormats, allowed = c("ofst", "sam", "bam")) {
  libFormats <- normalize_lib_formats(libFormats)
  format_values <- unlist(libFormats, use.names = FALSE)
  stopifnot(length(format_values) > 0)
  if (!all(format_values %in% allowed)) {
    stop(
      "Unsupported libFormats value. Allowed values are: ",
      paste(allowed, collapse = ", ")
    )
  }
}

write_ofst_library <- function(x, path) {
  export.ofst(x, file = path)
  c(default = path, ofst = path)
}

write_sam_library <- function(x, path, seqinfo = GenomeInfoDb::seqinfo(x)) {
  sam_path <- sub("\\.[^.]+$", ".sam", path)
  if (!grepl("\\.sam$", sam_path)) {
    sam_path <- paste0(sam_path, ".sam")
  }
  samFromGAlignment(x, path = sam_path, seqinfo = seqinfo, make_bam = FALSE)
  c(default = sam_path, sam = sam_path)
}

# Write the simulated reads as BAM, via a SAM file that is removed afterwards.
#
# BAM is compressed and indexed, so it is what downstream tools expect, but it is
# written by first producing plain-text SAM and converting. The intermediate file
# is deleted once the BAM exists, so a run leaves only the finished library.
write_bam_library <- function(x, path, seqinfo = GenomeInfoDb::seqinfo(x)) {
  bam_path <- sub("\\.[^.]+$", ".bam", path)
  if (!grepl("\\.bam$", bam_path)) {
    bam_path <- paste0(bam_path, ".bam")
  }
  sam_path <- sub("\\.bam$", ".sam", bam_path)
  samFromGAlignment(x, path = sam_path, seqinfo = seqinfo, make_bam = TRUE)
  if (file.exists(sam_path)) {
    unlink(sam_path)
  }
  c(default = bam_path, bam = bam_path)
}

write_simulated_library <- function(x, file_base, format,
                                    seqinfo = GenomeInfoDb::seqinfo(x)) {
  switch(
    format,
    ofst = write_ofst_library(x, paste0(file_base, ".ofst")),
    sam = write_sam_library(x, paste0(file_base, ".sam"), seqinfo = seqinfo),
    bam = write_bam_library(x, paste0(file_base, ".bam"), seqinfo = seqinfo),
    stop("Unsupported format: ", format)
  )
}



list_to_mat <- function(lengths, rnase_length) {
  length_max <- max(lengths)
  region_length_matrix <- lapply(lengths, function(x)
    c(rep(TRUE, x), rep(TRUE, rnase_length),
      rep(FALSE, length_max - x)))
  if (length(unique(lengths(region_length_matrix))) != 1)
    stop("list to mat failed to create square matrix,
         report on github!")
  matrix(unlist(region_length_matrix), nrow = length(lengths),
         byrow = TRUE)
}

matrix_row_indices <- function(row_lengths) {
  cbind(
    rep.int(seq_along(row_lengths), row_lengths),
    sequence(row_lengths)
  )
}

# Lay the per-gene weight vectors into one rectangular matrix.
#
# Genes differ in length, but the sampler wants a matrix, so every row is as wide
# as the longest gene and shorter genes are padded on the right. The padding is a
# vanishingly small positive number rather than zero, because these entries are
# used as Dirichlet concentrations and zero is not a valid one; the value is
# small enough that padded positions receive no reads in practice.
pack_alpha_rows <- function(alpha_rows, region_length_matrix, pad_value = 1e-24) {
  if (!length(alpha_rows)) {
    return(matrix(pad_value, nrow = 0L, ncol = ncol(region_length_matrix)))
  }

  row_lengths <- lengths(alpha_rows)
  alpha_matrix <- matrix(
    pad_value,
    nrow = length(alpha_rows),
    ncol = ncol(region_length_matrix)
  )
  alpha_matrix[matrix_row_indices(row_lengths)] <- unlist(alpha_rows, use.names = FALSE)
  alpha_matrix
}

flatten_sample_rows <- function(sample_matrix, row_lengths) {
  if (!length(row_lengths)) {
    return(numeric())
  }
  sample_matrix[matrix_row_indices(row_lengths)]
}

validate_dmn_alpha_scale <- function(scale) {
  if (!is.numeric(scale) || length(scale) != 1L || is.na(scale) ||
      !is.finite(scale) || scale <= 0) {
    stop("dmn_alpha_scale must be one finite number greater than zero", call. = FALSE)
  }
  invisible(scale)
}

# rnase_bias is keyed by library type only, while every other per-library
# argument of simNGScoverage() (ideal_coverage, auto_correlation, sampling) and
# simCountTablesRegions()'s region_proportion are nested by region. Passing a
# region-nested rnase_bias is therefore an easy mistake, and a silent one:
# rnase_bias[["RFP"]] is then NULL, which switches the RNase kernel off for
# every region instead of failing -- collapsing the simulated reading-frame
# distribution from roughly 5:1:2 to 100:0:0. Reject that shape explicitly, and
# reject non-numeric entries too, since a quoted expression is never evaluated
# for this argument (unlike ideal_coverage's).
validate_rnase_bias <- function(rnase_bias) {
  if (is.null(rnase_bias)) return(invisible(NULL))
  if (!is.list(rnase_bias) || is.null(names(rnase_bias)) ||
      any(!nzchar(names(rnase_bias)))) {
    stop("rnase_bias must be a named list with one entry per library type",
         call. = FALSE)
  }
  # Duplicate names are another silent way to lose the kernel: rnase_bias[["RFP"]]
  # returns the FIRST match, so list(RFP = NULL, RFP = c(1, 2, 1)) disables
  # smearing while looking like it enables it.
  duplicated_names <- unique(names(rnase_bias)[duplicated(names(rnase_bias))])
  if (length(duplicated_names)) {
    stop("rnase_bias must name each library type once, but got duplicates: ",
         paste(duplicated_names, collapse = ", "),
         ". Only the first entry per name is ever used.", call. = FALSE)
  }
  nested <- names(rnase_bias)[vapply(rnase_bias, is.list, logical(1))]
  if (length(nested)) {
    stop("rnase_bias ",
         if (length(nested) > 1L) "entries '" else "entry '",
         paste(nested, collapse = "', '"), "' ",
         if (length(nested) > 1L) "contain" else "contains",
         " a list, so this looks like a per-region list. rnase_bias is ",
         "the one argument here that is grouped by library type only, not by ",
         "region. Written per region it would be read as 'no kernel for RFP', ",
         "which switches RNase smearing off everywhere: the simulated coverage ",
         "then sits only on exact codon positions instead of spreading over all ",
         "three reading-frame positions. Use one kernel per library type, e.g. ",
         "list(RFP = c(0.5, 2, 1, 10, 2, 1, 0.5), RNA = 1).", call. = FALSE)
  }
  # Odd length: the kernel is applied center-aligned, and rnase_kernel_reach()
  # extends each gene by floor(length/2) rows per side. An even-length kernel has
  # no center, so sim_sequence_bias() returns one alpha value fewer than the row
  # count -- e.g. a 9 nt region with c(1, 1) yields 10 alphas against 11 rows.
  # An all-zero kernel (including the scalar 0, which looks like a natural way to
  # say "off") makes every smoothed weight 0, so the mean-preserving rescale
  # divides by zero and the whole alpha vector becomes NaN; that surfaces later as
  # a misleading "dmn_alpha_scale produced invalid Dirichlet alpha values".
  # Use 1 or NULL to disable smearing instead.
  # Each message says what is wrong, why it matters and what to do instead: the
  # reader of a validation error is by definition someone who does not know this
  # file, so "applied center-aligned" or "reach" would help nobody.
  problem <- vapply(names(rnase_bias), function(name) {
    kernel <- rnase_bias[[name]]
    label <- paste0("rnase_bias$", name)
    if (is.null(kernel)) return(NA_character_)
    if (!is.numeric(kernel) || !length(kernel)) {
      return(paste0(
        label, " is not a numeric vector. This argument takes plain numbers, one",
        " weight per position; unlike ideal_coverage and auto_correlation a",
        " quote(...) expression is never evaluated here, so it would silently",
        " switch RNase smearing off."
      ))
    }
    if (!all(is.finite(kernel))) {
      return(paste0(label, " contains NA, NaN or Inf. Every weight must be a",
                    " finite number."))
    }
    if (length(kernel) %% 2L != 1L) {
      return(paste0(
        label, " has ", length(kernel), " weights, which is an even number. The",
        " weights describe how far RNase digestion spreads a read's signal onto",
        " neighbouring positions, and they are applied centred on the read's own",
        " position -- so one weight has to land on that position, which needs an",
        " odd count. The weights themselves need not be symmetric: the default is",
        " not, and setting the weights on one side to zero spreads the signal in",
        " one direction only, which is a perfectly reasonable thing to ask for.",
        " What an even count breaks is the centring, so the simulated coverage",
        " would come out shifted by one nucleotide."
      ))
    }
    if (any(kernel < 0)) {
      return(paste0(
        label, " contains a negative weight. A weight is the share of a read's",
        " signal that lands on a neighbouring position, so it cannot be negative."
      ))
    }
    if (!any(kernel > 0)) {
      return(paste0(
        label, " has no positive weight, so every position would receive a share",
        " of zero and the coverage could not be rescaled (you would get NaN",
        " further down). To switch RNase smearing off for a library type, use 1",
        " or NULL rather than 0."
      ))
    }
    NA_character_
  }, character(1), USE.NAMES = FALSE)
  invalid <- !is.na(problem)
  if (any(invalid)) {
    stop(paste(problem[invalid], collapse = " "), call. = FALSE)
  }
  invisible(NULL)
}

resolve_dmn_alpha_scale <- function(scale, seq_bias) {
  if (!is.null(scale)) {
    validate_dmn_alpha_scale(scale)
    return(scale)
  }
  learned <- if (!is.null(seq_bias) && "dmn_alpha_scale" %in% names(seq_bias)) {
    unique(seq_bias$dmn_alpha_scale[!is.na(seq_bias$dmn_alpha_scale)])
  } else numeric()
  if (length(learned) > 1L) {
    stop("seq_bias contains multiple dmn_alpha_scale values", call. = FALSE)
  }
  if (!length(learned)) return(1)
  validate_dmn_alpha_scale(learned)
  learned
}

# "AUTO" (simNGScoverage()'s seq_bias default) picks the bundled table whose
# `shift` matches fragment_geometry$site_reference, so the codon-bias table
# and the simulated-RPF fragment placement describe the same ribosome site by
# default (see simNGScoverage()'s seq_bias/fragment_geometry docs). site_reference
# only affects simulated_rpf's physical fragment placement, so "AUTO" keeps the
# historical P-site table for every other fragment_mode. Matches the
# true_uorf_ranges = "AUTO" sentinel convention already used by
# simNGScoverage(); any other seq_bias value -- including NULL, which
# disables codon bias -- is always used as-is, and (unlike a missing()-based
# check) this still works when a wrapper function redeclares its own
# seq_bias = "AUTO" default and forwards it explicitly.
resolve_seq_bias <- function(seq_bias, fragment_mode, site_reference) {
  if (!identical(seq_bias, "AUTO")) return(seq_bias)
  # site_reference only governs simulated_rpf's physical fragment placement;
  # every other fragment_mode keeps the historical P-site table.
  shift <- if (fragment_mode == "simulated_rpf") site_reference else "p_site"
  load_seq_bias(shift = gsub("_", "-", shift, fixed = TRUE))
}

scale_dmn_alpha <- function(alpha_rows, scale) {
  validate_dmn_alpha_scale(scale)
  scaled <- lapply(alpha_rows, function(alpha) alpha * scale)
  valid <- vapply(scaled, function(alpha) {
    length(alpha) > 0L && all(is.finite(alpha)) && all(alpha > 0)
  }, logical(1))
  if (!all(valid)) {
    stop("dmn_alpha_scale produced invalid Dirichlet alpha values", call. = FALSE)
  }
  scaled
}

is_list_or_null <- function(...) {
  args <- list(...)
  mc <- match.call(expand.dots = FALSE)
  for (i in seq_along(args)) {
    l <- args[[i]]
    if(!(is.list(l) | is.null(l))) {
      stop("Argument: '", mc$...[[i]], "' is neither list or null!")
    }
  }
}

get_value <- function(x, region, libtype) {
  res <- try(x[[region]][[libtype]], silent = TRUE)
  if (is(res, "try-error")) res <- NULL
  return(res)
}

assay_by_chromo <- function(assay, seqnamesPer) {
  # as.matrix(): as.data.table() on a single-column DelayedMatrix silently
  # drops its column name (its as.array(x, drop = TRUE) collapses to a bare
  # vector, so data.table falls back to naming the column after that
  # deparsed expression instead) -- only reproducible with exactly one
  # column, since a wider DelayedMatrix keeps its real dimensions and names.
  # A plain base matrix (or a DelayedMatrix coerced to one) isn't affected.
  assay_by_chromosome <- as.data.table(as.matrix(assay))
  assay_by_chromosome <- if (nrow(assay_by_chromosome) == 1) {
    # Column name must match the `else` branch's own grouping column
    # (seqnamesPer, not seqnamesPerGroup) -- the caller
    # (assay_by_chromosome$seqnamesPer in simNGScoverage()) currently only
    # resolves correctly here by data.table's $ partial-name matching.
    cbind(seqnamesPer = seqnamesPer, assay_by_chromosome)
  } else {
    assay_by_chromosome[, lapply(.SD, sum, na.rm=TRUE), by=.(seqnamesPer)]
  }
}

# Check the caller's arguments before any work starts.
#
# Everything here is cheap to test and expensive to get wrong: a simulation that
# fails after an hour because a region name was misspelled costs far more than
# the check. The function deliberately runs inside its caller's environment, so
# it sees that call's arguments directly and can also fill in values the caller
# then uses -- uorf_ranges and regionsToSample are set here, not just verified.
input_validation_controller <- function() {
  with(rlang::caller_env(), {
    message("- Validating input")
    # Sanity test of input
    stopifnot(dir.exists(out_dir))
    stopifnot(is(simGenome, "character") & (c("genome", "gtf", "txdb") %in% names(simGenome)))
    stopifnot(all(transcripts %in% names(count_table)))
    is_list_or_null(ideal_coverage, rnase_bias, auto_correlation,
                    sampling, read_lengths_per)
    stopifnot(is(count_table, "SummarizedExperiment"))
    if (!is.null(seq_bias))
      stopifnot(c("seqs", "alpha") %in% colnames(seq_bias))
    validate_lib_formats(libFormats)
    all_allowed_regions <- c("leader", "cds", "trailer", "uorf")
    regionsToSample <- assayNames(count_table)[-1]
    stopifnot(all(regionsToSample %in% all_allowed_regions))

    if ("uorf" %in% regionsToSample) {
      if (is.character(true_uorf_ranges)) {
        uorf_ranges <- readRDS(file.path(dirname(simGenome["genome"]), "true_uORFs.rds"))
      } else {
        uorf_ranges <- true_uorf_ranges
      }
      if (any((widthPerGroup(uorf_ranges, FALSE) %% 3) != 0)) {
        warning("Detected uORF ranges that ends on incomplete codon (is not %% 3 == 0 in length")
      }
      uorf_prop_mode <- ifelse(is(uorf_prop_within_gene, "character"),
                               "character", "numeric")
      if (uorf_prop_mode == "numeric") {
        stopifnot(length(uorf_ranges) == length(uorf_prop_within_gene))
        stopifnot(all(txNames(uorf_ranges) == names(uorf_prop_within_gene)))
      } else stopifnot(uorf_prop_within_gene %in% c("uniform", "length"))
    }
  })
}

sequence_table_controller <- function() {
  with(rlang::caller_env(),{
    message("- Setting up sequence tables")
    # Create tiling for ranges
    for (region in regionsToSample) {
      region_ranges <- get(paste0(region, "_ranges"), mode = "S4")
      lengths <- widthPerGroup(region_ranges, FALSE)
      assign(paste0("lengths_", region), lengths)
      if (is(region_ranges, "GRanges")) {
        tile <- tile(region_ranges, width = 1)
        tile <- sortPerGroup(tile, quick.rev = TRUE)
      } else {
        region_ranges_noname <- region_ranges
        names(region_ranges_noname) <- NULL
        tile <- tile1(region_ranges_noname, sort.on.return = TRUE)
      }
      dt_range <- data.table::data.table(seqnames = as.character(unlist(seqnames(tile), use.names = FALSE)),
                                         start = unlist(start(tile), use.names = FALSE),
                                         end = unlist(start(tile), use.names = FALSE),
                                         strand = as.character(unlist(strand(tile), use.names = FALSE)))
      tile_groups <- groupings(tile)
      transcript_ids <- txNames(region_ranges)
      dt_range[, transcript_id := transcript_ids[tile_groups]]
      dt_range[, region_position := seq_len(.N), by = tile_groups]
      alpha_matrix <- NULL
      add_sequence_bias <- region %in% c("cds", "uorf") & !is.null(seq_bias)
      # Only DMN sampling (via sim_sequence_bias()'s RNase-kernel smoothing,
      # which grows each gene's alpha vector by the kernel's extra reach)
      # consumes the RNase-extended rows added below. MN sampling never
      # applies rnase_bias at all and only ever produces one score per
      # original position, so extending dt_range for a region no libtype
      # samples with DMN would only add rows nothing fills in, breaking the
      # score-assignment length match later in nt_coverage_all_regions().
      # Resolve each libtype actually present in count_table's per-sample
      # sampling mode the same way nt_coverage_all_regions() itself does
      # (get_value() falling back to "MN") -- not just the raw sampling[[region]]
      # list -- so a libtype named in `sampling` that isn't actually part of
      # this experiment is ignored, and a real libtype silently defaulting to
      # MN (because it's missing from `sampling`) is still accounted for.
      # Also skip a libtype whose samples all have zero counts *for this
      # region*: nt_coverage_all_regions() itself never samples those (its
      # `if (sum(region_counts) > 0)` guard), so they can never hit the
      # length-mismatch this check exists to catch, e.g. an RNA library with
      # reads only in the leader while cds is being set up here.
      # as.matrix(): assay() preserves whatever matrix-like class the assay
      # was stored as (e.g. a DelayedMatrix), which base::colSums() (the
      # unqualified name in scope here) does not accept.
      region_assay <- as.matrix(assay(count_table, region))
      sample_libtypes <- as.character(colData(count_table)$libtype)
      region_libtypes <- unique(sample_libtypes[colSums(region_assay) > 0])
      region_modes <- vapply(region_libtypes, function(lt) {
        mode <- get_value(sampling, region, lt)
        if (is.null(mode)) "MN" else mode
      }, character(1))
      region_uses_dmn <- any(region_modes == "DMN")
      if (region_uses_dmn) {
        rnase_extra <- max(0, length(rnase_bias[["RFP"]]) - 1)
        region_length_matrix <- list_to_mat(lengths, rnase_extra)
        assign(paste0("region_length_matrix", region), region_length_matrix)
      }

      if (add_sequence_bias) {
        dt_range[, genes := groupings(tile)]
        alpha_matrix <- add_sequence_bias(simGenome, dt_range, seq_bias,
                                          region_ranges, lengths, region)
      }
      needs_rnase_extension <- !is.null(rnase_bias[["RFP"]]) &&
        (region %in% c("cds", "uorf")) && region_uses_dmn &&
        rnase_kernel_reach(rnase_bias) > 0L
      # A non-DMN libtype sharing this region's (possibly RNase-extended)
      # dt_range never applies rnase_bias itself, so its MN branch in
      # nt_coverage_all_regions() zero-pads its own sample vector by this
      # same reach at both ends of every gene (matching the prefix/original/
      # suffix row layout append_rnase_to_dt()/append_rnase_to_simulated_rpf_
      # table() both produce) instead of contributing invented smeared reads
      # to the flanking positions it doesn't model.
      assign(
        paste0("region_rnase_reach_", region),
        if (needs_rnase_extension) rnase_kernel_reach(rnase_bias) else 0L
      )
      if (needs_rnase_extension) {
        if (is.null(dt_range$genes)) {
          dt_range[, genes := groupings(tile)]
          dt_range[, position := seq_len(.N), by = genes]
        }
        dt_range <- if (exists("fragment_mode") && fragment_mode == "simulated_rpf") {
          append_rnase_to_simulated_rpf_table(dt_range, rnase_bias, models)
        } else {
          append_rnase_to_dt(dt_range, lengths, rnase_bias)
        }
      }
      dt_range[, signal_position := start]


      assign(paste0("seq_alpha_list_", region), alpha_matrix)
      assign(paste0("dt_", region), dt_range)
    }})
}

# Insist on exactly one weight profile before simulating.
#
# load_seq_bias(bias = "all") returns all ten measured libraries stacked in one
# table, which is useful for comparing them in a figure but meaningless as a
# simulation input: there is no way to tell which library a given position should
# follow. Rejecting it here turns what would otherwise be a silently wrong
# simulation into a message that names the alternative.
validate_sequence_profile <- function(seq_bias) {
  if (is.null(seq_bias) || is.null(seq_bias$variable)) return(invisible(NULL))
  profiles <- unique(as.character(seq_bias$variable))
  if (length(profiles) != 1L || anyNA(profiles) || !nzchar(profiles)) {
    stop(
      "seq_bias must contain exactly one named sequence profile. ",
      "Select one with load_seq_bias(bias = 'median') or subset your table; ",
      "bias = 'all' is for inspection, not simulation.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

# Look up a weight for every codon (or amino acid) in a region's actual sequence.
#
# Translation does not run at a constant speed. How long a ribosome dwells over a
# codon depends on how readily the matching tRNA is available, so each codon
# carries its own weight; `tAI` is that lookup table, one positive weight per
# motif. This function reads the region's real sequence out of the genome,
# replaces each motif by its weight, and hands back one weight vector per gene,
# ready to be turned into per-nucleotide weights.
#
# Four positions are treated as motifs of their own, marked by the symbols #, %,
# & and *: the start codon, the one after it, the one before the stop, and the
# stop itself. Initiation and termination are slower and differently regulated
# than elongation, so those positions need weights that are not tied to which
# codon happens to sit there. The symbols only appear if the supplied table
# defines them, so a table without them simply treats those positions normally.
#
# The checks here are deliberately strict and fail rather than guess: a motif
# present in the sequence but missing from the table would otherwise silently
# shift every later weight onto the wrong position.
add_sequence_bias <- function(simGenome, dt_range, tAI, region_ranges, lengths, region) {
  validate_sequence_profile(tAI)
  tAI <- data.table::copy(tAI)
  if (!is.null(tAI$variable)) {
    tAI$variable <- NULL
  }
  if (is.factor(tAI$seqs)) tAI[, seqs := as.character(seqs)]
  if (anyNA(tAI$seqs) || anyDuplicated(tAI$seqs) || anyNA(tAI$alpha) ||
      any(!is.finite(tAI$alpha)) || any(tAI$alpha <= 0)) {
    stop("seq_bias needs exactly one finite, positive alpha value per motif",
         call. = FALSE)
  }
  unique_seqs <- tAI$seqs
  by_AA <- all(nchar(unique_seqs) == 1)
  by_codon <- any(nchar(unique_seqs) == 3)
  if (by_AA) {
    special_symbols <- c("#", "%", "&", "*") %in% unique_seqs
    as_seq <- "AA"
  } else if (by_codon) {
    special_symbols <- c("###", "%%%", "&&&", "***") %in% unique_seqs
    as_seq <- "codon"
  } else stop("Malformed format of tAI table of codon/AA scores")
  message("-- Biasing ", as_seq, " estimators from: ", region)

  seqs <- translate_orf_seq(region_ranges, simGenome["genome"], is.sorted = TRUE,
                            as = as_seq,
                            start.as.hash = special_symbols[1], startp1.as.per = special_symbols[2],
                            stopm1.as.amp = special_symbols[3], return.as.list = TRUE)
  missing_motifs <- setdiff(unique(unlist(seqs, use.names = FALSE)), tAI$seqs)
  if (length(missing_motifs)) {
    stop("seq_bias does not contain every ", as_seq, " found in the ", region,
         " sequences (missing: ", paste(head(missing_motifs, 10), collapse = ", "),
         "). Alpha values would be assigned to the wrong positions.", call. = FALSE)
  }
  dt_range[, position := seq_len(.N), by = genes]
  tAI_short <- tAI[, c("seqs", "alpha")]
  tAI_merged <- data.table::merge.data.table(data.table(seqs = unlist(seqs, use.names = FALSE)),
                                             tAI_short, by = "seqs", sort = FALSE)
  # Create the alpha matrix for region
  seq_alpha <- tAI_merged$alpha[!is.na(tAI_merged$alpha)]
  if (length(seq_alpha) != sum(lengths) / 3) {
    stop("seq_bias contains missing alpha values for motifs in the ", region,
         " sequences.", call. = FALSE)
  }
  # One weight per codon, but dt_range has one row per nucleotide, so every third
  # row names the gene that codon belongs to.
  seq_alpha <- split(seq_alpha, dt_range$genes[c(T, F, F)])
  # seq_lengths <- lengths / 3
  # if (any(seq_lengths != as.integer(lengths/3))) stop("Mismatch of seqlength and divisor")
  # seq_length_max <- max(seq_lengths)
  # seq_alpha_max <- lapply(seq_alpha, function(x)
  #   c(x, rep(1e-24, length.out = seq_length_max - length(x))))
  # alpha_matrix <- t(matrix(unlist(seq_alpha_max, F, F), nrow = seq_length_max))
  seq_alpha
}


apply_autocorrelation_kernel <- function(signal, kernel) {
  if (!is.numeric(kernel) || !length(kernel) || length(kernel) %% 2L != 1L ||
      any(!is.finite(kernel) | kernel < 0) || !any(kernel > 0)) {
    stop("Numeric auto_correlation must be an odd, non-negative kernel")
  }
  radius <- (length(kernel) - 1L) %/% 2L
  lags <- seq.int(-radius, radius)
  vapply(seq_along(signal), function(i) {
    positions <- i + lags
    keep <- positions >= 1L & positions <= length(signal)
    sum(signal[positions[keep]] * kernel[keep]) / sum(kernel[keep])
  }, numeric(1))
}

# Build the per-nucleotide weights that decide where reads land inside a region.
#
# The biology this walks through, in order. A ribosome does not sit evenly along a
# transcript: it lingers at some codons and hurries past others, so each codon
# attracts a different amount of signal. That is what `alpha_matrix` holds -- one
# weight per codon, larger meaning more signal. Four things then happen to those
# weights, and each has a physical reason:
#
#   1. Neighbouring codons influence each other. A ribosome covers roughly 30
#      nucleotides, so what is happening a few codons away is not independent of
#      the A-site. `seq_acf` spreads each codon's weight onto its neighbours.
#   2. The ribosome steps one whole codon at a time, so signal piles up on one
#      nucleotide of each triplet and not on the other two. The weights move from
#      codon level to nucleotide level as (weight, 0, 0), the ideal triplet
#      translocation pattern. This is what makes the 3-nt periodicity of
#      Ribo-seq.
#   3. RNase does not cut at exactly the same place every time. A read's signal
#      therefore lands slightly before or after where the ribosome actually was,
#      spread over a few neighbouring positions. `rnase_acf` is that spread.
#   4. All of the above only redistributes signal; it must not create or destroy
#      any. The last step rescales each region back to its original mean.
#
# Returns one numeric vector per region, used as the Dirichlet concentration
# (alpha) when the region's read budget is distributed over its positions.
#
# A warning about `seq_acf`, because its meaning depends on its TYPE:
#   * a numeric vector is used directly as kernel weights. A single number is
#     therefore a kernel of width one, which averages each position with itself
#     and changes nothing -- `seq_acf = 9` does NOT mean "9 codons each side".
#   * an unevaluated expression (what `shapes(9)` returns) is evaluated below,
#     and only that branch also applies the rescaling in step 1b.
sim_sequence_bias <- function(ideal_coverage, lengths, alpha_matrix,
                              seq_acf = 9, rnase_acf =
                                c(0.5,1,2,6,2,1,0.5)) {
  if (!is.null(alpha_matrix)){
    # res <- as.data.table(t(alpha_matrix))
    res <- alpha_matrix
    alpha_means <- unlist(lapply(alpha_matrix, mean), use.names = FALSE)
  } else {
    res <- lapply(lengths, function(x) eval(ideal_coverage))
    alpha_means <- rep(mean(res[[1]]), length(res))
  }
  # Step 1: smear each codon's weight onto its neighbouring codons.
  if (is.numeric(seq_acf)) {
    res <- lapply(res, apply_autocorrelation_kernel, kernel = seq_acf)
  } else if (!is.null(seq_acf)) { # Higher order auto correlation
    # Step 1b: sharpen the contrast between codons before smearing, so that
    # smearing does not flatten the differences away entirely. Two properties of
    # this step are worth knowing before changing anything here. It raises the
    # weights to a fixed power and then divides by a quantile of themselves, so
    # its effect on the spread depends on the gene; and the quantile is drawn at
    # random on every call, so repeated runs do not produce identical alphas
    # unless the random seed is fixed. Both are deliberate and long-standing.
    if (!is.null(alpha_matrix)) { # Rescale alpha values
      scalers <- unlist(lapply(res, function(x) {
        codon_extreme <- max(x) / median(x)
        codon_variance <- sd(x) / median(x)
        30*(codon_variance / codon_extreme)
      }), use.names = FALSE)

      # Scale signal before auto correlation smearing
      res <- lapply(res, function(x) {
        x**1.25
      })
      quantile <- sample(seq(7, 9), length(res), TRUE, prob = c(0.4, 0.5, 0.1))/10
      # Sample number of extreme peaks
      # quants <- unlist(alpha_matrix[, lapply(.SD, function(x) quantile(x, 0.8))],
      #                  use.names = FALSE)
      res <- lapply(seq_along(res), function(x) {
        (res[[x]] /
          (quantile(res[[x]], quantile[x])))**(scalers[x])
      })
    }
    if (any(c("sin","cos") %in% as.character(quote(1)))) stop("Implement me")
    res <- lapply(res, function(alpha_vec)
      eval(seq_acf))
  }
  # Step 2: move from one weight per codon to one weight per nucleotide. The
  # ribosome translocates a whole codon at a time, so all of a codon's weight is
  # placed on its first nucleotide and the other two are left at zero. Reading
  # across a transcript this gives the (1, 0, 0) pattern that produces Ribo-seq's
  # 3-nt periodicity.
  if (!is.null(alpha_matrix)) { # Codon to NT level
    res <- lapply(res, function(x) {
      weights_per_nucleotide <- rep(x, each = 3)
      weights_per_nucleotide[seq(weights_per_nucleotide) %% 3 %in% c(2,0)] = 0
      weights_per_nucleotide
    })
  }


  # Step 3: spread each position's weight onto its immediate neighbours, because
  # the enzyme cuts slightly before or after the ribosome's true position. The
  # weights are applied centred, which is why an even number of them has no
  # middle and would shift the whole profile by one nucleotide (see
  # validate_rnase_bias()). The zero padding lets the first and last real
  # positions be smeared like any other instead of being cut short.
  if (!is.null(rnase_acf)) { # Lower order auto correlation
    edge_padding <- rep(0, length(rnase_acf) - 1)
    res <- lapply(res, function(x) {
      x <- c(edge_padding, x, edge_padding)
      frollapply_compat(x,
                        window = length(rnase_acf),
                        FUN = function(i) sum(i*rnase_acf),
                        align = "center", fill = NA)
    })
  }
  # Step 4: restore the scale. Smearing moved signal around; rescaling to the
  # original mean makes sure it neither created nor destroyed any. Exact zeros
  # are then nudged to a tiny positive number, because a Dirichlet concentration
  # of zero is not a valid parameter -- it would assert that a position can never
  # receive a read, which is stronger than "very unlikely" and is not what a
  # weight of zero is meant to say here.
  # Cleanup
  res <- lapply(seq_along(res), function(i) {
    codon_ac_rnase_alphas <- res[[i]]
    codon_ac_rnase_alphas <- codon_ac_rnase_alphas[!is.na(codon_ac_rnase_alphas)]
    codon_ac_rnase_alphas <- codon_ac_rnase_alphas * (alpha_means[i] / mean(codon_ac_rnase_alphas))
    codon_ac_rnase_alphas[codon_ac_rnase_alphas == 0] <- 1e-24
    return(codon_ac_rnase_alphas)
  })

  return(res)
}

#' Fetch internal sequence bias tables
#'
#' Translation does not run at a constant speed: a ribosome dwells longer over
#' some codons than others. These tables give one weight per sequence motif,
#' measured from ten real human Ribo-seq libraries, and are used as the
#' Dirichlet concentration that decides how strongly a position attracts reads.
#'
#' The tables are read from files named
#' \code{paste0(type, "_bias_", shift, "_estimates_human.csv")} in \code{dir}.
#' @param type character, default "AA". How finely the sequence is described.
#' \code{"codon"} gives one weight per triplet; \code{"AA"} pools the triplets
#' that code for the same amino acid, which is steadier when a library is small
#' but cannot distinguish synonymous codons.
#' @param shift character, default "p-site". Alternative: "a-site". Which site
#' inside the ribosome the weights were referenced to. The P-site holds the
#' growing peptide chain and the A-site receives the incoming tRNA; they sit one
#' codon apart, so a table referenced to one site is displaced by three
#' nucleotides if used as though it were referenced to the other. Match this to
#' wherever the simulation places its reads.
#' @param dir Directory with sequence biases, default is internal path
#' predefined estimators: system.file(package = "coverageSim", "extdata")
#' @param bias The default, \code{"median"}, calculates the motif-wise median
#' of the alpha values in libraries R1--R10. Alternatively, use one profile
#' name (e.g. "R2"), one of the aliases "start_codon" (R2), "stop_codon"
#' (R1), and "similar" (R10), or \code{"all"} to load every profile for
#' inspection. Select one profile before passing a table loaded with
#' \code{bias = "all"} to \code{simNGScoverage()}.
#' @return a data.table of bias per sequence motif
#' @export
#' @examples
#' load_seq_bias()
#' load_seq_bias(type = "codon")
#' load_seq_bias(bias = "all")
# The three aliases name the feature each library shows most strongly: R2 a
# pronounced start-codon peak, R1 a pronounced stop-codon one, R10 the least
# distinctive of the ten. "all" is rejected downstream by
# validate_sequence_profile(), since a simulation cannot follow ten profiles at
# once.
load_seq_bias <- function(type = "AA", shift = "p-site",
                          dir = system.file(package = "coverageSim", "extdata"),
                          bias = "median") {
  stopifnot(type %in% c("AA", "codon"))
  stopifnot(shift %in% c("p-site", "a-site"))
  if (!is.character(bias) || length(bias) != 1L || is.na(bias) || !nzchar(bias)) {
    stop("bias must be one profile name, alias, or 'all'", call. = FALSE)
  }
  shift <- gsub("-", "_", shift)
  file <- paste0(type, "_bias_", shift, "_estimates_human.csv")
  dt <- fread(file.path(dir, file))
  if (bias == "all") return(dt)
  if (bias == "median") return(median_sequence_profile(dt))
  aliases <- c(start_codon = "R2", stop_codon = "R1", similar = "R10")
  profile <- if (bias %in% names(aliases)) unname(aliases[[bias]]) else bias
  if (!profile %in% dt$variable) {
    stop(
      "Unknown bias profile: ", bias, ". Choose one of: ",
      paste(c("median", names(aliases), unique(dt$variable), "all"), collapse = ", "),
      call. = FALSE
    )
  }
  dt[variable == profile, ]
}

# Combine the ten measured libraries into one representative weight table.
#
# The bundled codon and amino-acid weights were measured from ten real human
# Ribo-seq libraries, R1 to R10. Any single one of them carries that experiment's
# own quirks -- one has a pronounced start-codon peak, another an unusual stop
# signal. Taking the median per motif keeps what the libraries agree on and drops
# what only one of them shows, which is what makes it a sensible default.
#
# The median is taken motif by motif, so the result is not any one library but a
# consensus. All ten must be present and describe the same motifs, or the
# comparison would be between different things; that is what the checks enforce.
median_sequence_profile <- function(profiles) {
  required_profiles <- paste0("R", seq_len(10L))
  selected <- profiles[variable %in% required_profiles]
  if (!setequal(unique(selected$variable), required_profiles)) {
    stop("Median sequence bias requires profiles R1 through R10", call. = FALSE)
  }
  if (anyNA(selected$seqs) || anyNA(selected$alpha) ||
      any(!is.finite(selected$alpha)) || any(selected$alpha <= 0)) {
    stop("Sequence-bias alpha values must be finite and positive", call. = FALSE)
  }
  support <- selected[, .N, by = .(variable, seqs)]
  if (any(support$N != 1L)) {
    stop("Each sequence motif must occur once in every profile", call. = FALSE)
  }
  motif_support <- support[, .N, by = seqs]
  if (any(motif_support$N != length(required_profiles))) {
    stop("Profiles R1 through R10 must contain the same sequence motifs", call. = FALSE)
  }
  motif_order <- unique(selected$seqs)
  result <- selected[, .(alpha = stats::median(alpha)), by = seqs]
  result <- result[match(motif_order, seqs)]
  result[, variable := "median"]
  data.table::setcolorder(result, c("variable", "seqs", "alpha"))
  result[]
}

# Give each gene a margin at both ends so RNase smearing has somewhere to go.
#
# The smearing spreads a position's signal onto its neighbours. Near a gene's
# first and last nucleotide some of those neighbours lie outside the gene, and
# without room for them the outermost positions would be smeared with less signal
# than the rest and end up artificially quiet. Adding as many extra positions at
# each end as the kernel reaches lets every real position be treated alike.
append_rnase_to_dt <- function(dt_range, lengths, rnase_bias) {
  gene_split_sites_end <- cumsum(lengths)
  gene_split_sites_start <- c(1, (gene_split_sites_end + 1)[-length(lengths)])
  gene_split_sites <- sort(c(gene_split_sites_start, gene_split_sites_end))
  rnase_reach <- rnase_kernel_reach(rnase_bias)
  gene_split_sites <- rep.int(gene_split_sites, rnase_reach)
  new_index_map <- sort(c(seq.int(nrow(dt_range)), gene_split_sites))
  dt_range <- dt_range[new_index_map, ]
  # Extend each gene by the RNase reach at its 5' and 3' end, in transcript
  # direction: genomic coordinates run backwards along minus-strand genes.
  extension <- seq.int(rnase_reach, 0L)
  direction <- function(strand) data.table::fifelse(strand == "-", -1L, 1L)
  dt_range[position == 1, start := start - direction(strand) * extension]
  dt_range[dt_range[, .I[position == max(position)], by = genes]$V1,
           start := start + direction(strand) * rev(extension)]
  dt_range[, end := start]
  return(dt_range)
}

#' Get shape function
#'
#' Returns a quoted expression describing the codon-level autocorrelation
#' smoothing applied along a gene: a rolling window (\code{autocor_window()},
#' width \code{2*i + 1}, its \code{max.lag} argument) that locally correlates
#' neighboring codon-bias values, modeling real tRNA/wobble-position sharing
#' between nearby codons. Larger \code{i} smooths over a wider
#' neighborhood, which matters most for sparsely sampled genes: at low read
#' depth it noticeably increases how much reads cluster into a few
#' positions rather than spreading out (higher skew and peak-to-median
#' ratio), while at high read depth the effect is negligible.
#' @param i integer, default 9 (the same default \code{simNGScoverage()}
#'  itself uses). Half-width of the local autocorrelation smoothing window
#'  (the full window is \code{2*i + 1} codons). \code{i = 0} disables only
#'  this smoothing step; codon bias and RNase-kernel processing elsewhere
#'  in the pipeline still apply, so this alone does not produce a fully
#'  uniform coverage model.
#' @return a quote object to be evaluated inside the coverage function, or
#'  NULL when \code{i = 0}.
#' @export
shapes <- function(i = 9) {
  if (i == 0) return(NULL)
  bquote(autocor_window(alpha_vec, .(i), padding.rm = T))
}

# Sequence lengths for SAM/BAM headers: use the reference FASTA where the
# annotation carries none (a TxDb built from a GTF has NA lengths).
# Recover chromosome lengths that the annotation did not carry.
#
# An alignment file has to declare how long each reference sequence is, but a
# TxDb built from a GTF often does not know: a GTF describes features, not the
# genome they sit on. The FASTA index does know, so the missing lengths are
# filled in from there. Names the FASTA does not have are left as they are rather
# than guessed.
fill_missing_seqlengths <- function(seqinfo, fasta_file) {
  lengths <- GenomeInfoDb::seqlengths(seqinfo)
  missing <- is.na(lengths) | lengths <= 0L
  if (!any(missing)) return(seqinfo)
  index <- Rsamtools::scanFaIndex(fasta_file)
  reference <- stats::setNames(GenomicRanges::width(index),
                               as.character(GenomicRanges::seqnames(index)))
  fill <- names(lengths)[missing]
  known <- fill %in% names(reference)
  lengths[fill[known]] <- reference[fill[known]]
  GenomeInfoDb::Seqinfo(GenomeInfoDb::seqnames(seqinfo), unname(lengths),
                        GenomeInfoDb::isCircular(seqinfo),
                        GenomeInfoDb::genome(seqinfo))
}

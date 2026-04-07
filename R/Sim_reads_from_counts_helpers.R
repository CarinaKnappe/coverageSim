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
  assay_by_chromosome <- as.data.table(assay)
  assay_by_chromosome <- if (nrow(assay_by_chromosome) == 1) {
    cbind(seqnamesPerGroup = seqnamesPer, assay_by_chromosome)
  } else {
    assay_by_chromosome[, lapply(.SD, sum, na.rm=TRUE), by=.(seqnamesPer)]
  }
}

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
    if (!all(unlist(libFormats) == "ofst")) stop("Only 'ofst' is supported as NGS file format currently!")
    all_allowed_regions <- c("leader", "cds", "trailer", "uorf")
    regionsToSample <- assayNames(count_table)[-1]
    stopifnot(all(regionsToSample %in% all_allowed_regions))

    if ("uorf" %in% regionsToSample) {
      if (is.character(true_uorf_ranges)) {
        uorf_ranges <- readRDS(file.path(dirname(simGenome["genome"]), "true_uORFs.rds"))
      } else {
        uorf_ranges <- true_uorf_ranges
      }
      uorf_gene_grouping <- chmatch(txNames(uorf_ranges),
                                    rownames(assay(count_table, "uorf")))
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
      alpha_matrix <- NULL
      add_sequence_bias <- region %in% c("cds", "uorf") & !is.null(seq_bias)
      if (unlist(sampling[[region]])[1] == "DMN") {
        rnase_extra <- max(0, length(rnase_bias[["RFP"]]) - 1)
        region_length_matrix <- list_to_mat(lengths, rnase_extra)
        assign(paste0("region_length_matrix", region), region_length_matrix)
      }

      if (add_sequence_bias) {
        dt_range[, genes := groupings(tile)]
        alpha_matrix <- add_sequence_bias(simGenome, dt_range, seq_bias,
                                          region_ranges, lengths, region)
      }
      if (!is.null(rnase_bias[["RFP"]]) & (region %in% c("cds", "uorf"))) {
        if (is.null(dt_range$genes)) {
          dt_range[, genes := groupings(tile)]
          dt_range[, position := seq_len(.N), by = genes]
        }
        dt_range <- append_rnase_to_dt(dt_range, lengths, rnase_bias)
      }


      assign(paste0("seq_alpha_list_", region), alpha_matrix)
      assign(paste0("dt_", region), dt_range)
    }})
}

nt_coverage_all_regions_old <- function(count_table_regions, libClass, assay,
                                    ideal_coverage,
                                    rnase_bias, auto_correlation,
                                    read_lengths_per, uorf_ranges,
                                    uorf_prop_mode, debug_coverage) {
  data.table::rbindlist(lapply(regionsToSample, function(region) {
    region_counts <- dt[, region]

    if (sum(region_counts) > 0) {
      lengths <- get(paste0("lengths_", region))

      if (region == "uorf") {
        region_counts <-  distribute_reads_to_uORFs(region_counts, assay,
                                                    uorf_ranges, uorf_prop_mode)
      }
      region_ranges <- get(paste0(region, "_ranges"), mode = "S4")
      dt_region <- copy(get(paste0("dt_", region)))
      ideal_cov <- ideal_coverage[[region]][[libtype]] # skew_cov
      auto_cor <- auto_correlation[[region]][[libtype]] # bump_shape
      read_lengths <- unlist(read_lengths_per[libtype], use.names = FALSE)
      if (debug_coverage) browser()
      if (!is(ideal_coverage, "call")) {
        # Initial division of counts per codon

        alpha_matrix <- get(paste0("seq_alpha_matrix_", region))
        seq_lengths <- lengths/3
        res <- if (!is.null(alpha_matrix)) {
          split_prop <- if (region == "uorf") {groupings(seasonality)} else dt_region$genes
          # Bias by codon
          seq_bias <- extraDistr::rdirmnom(n = nrow(alpha_matrix),
                                           size = region_counts,
                                           alpha = alpha_matrix)
          seq_bias_t <- t(seq_bias)[t(alpha_matrix) != 1e-24]

          split(seq_bias_skewed, split_prop)
        } else if (length(skew_cov) != 1) {
          skew_alpha <- if (is(skew_cov, "data.frame")) {
            stopifnot(nrow(skew_cov) == 3)
            skew_cov$alpha
          } else {
            stopifnot(length(skew_cov) == 3 & is.numeric(skew_cov))
            skew_cov
          }
          seq_bias_skewed <- extraDistr::rdirmnom(n = sum(seq_lengths),
                                                  size = 100,
                                                  alpha = skew_alpha)
          seq_bias_skewed <- as.vector(t(seq_bias_skewed))
          split(seq_bias_skewed, split_prop)
        } else NULL


        use_trend <- is.call(bump_shape)
        if (use_trend) {
          b <- sample(seq(3,15, by = 0.1), size = length(lengths), replace = TRUE)
          pi_s <- sample(pi*seq(5), prob = 1/seq(5), size = length(lengths), replace = TRUE)
          periods <-  lapply(seq_along(lengths), function(i) seq(0, pi_s[i], pi_s[i]/lengths[i]))
          trend <- lapply(seq_along(lengths),
                          function(i, l = lengths[i], x = periods[[i]])
                            rep_len(eval(bump_shape), l))
          trend <- lapply(trend, function(shape) shape / max(shape))
          res <- if (is.null(res)) {
            trend
          } else split(unlist(res, use.names = F) * unlist(trend, use.names=F),
                       groupings(res))
        }
        cov <- lapply(seq_along(region_ranges),
                      function(x) a <- as.vector(rmultinom(1, region_counts[x], res[[x]] + 1e-24)))
      } else {
        cov <- lapply(seq_along(region_ranges),
                      function(y, res, x = lengths[y]) a <- as.vector(rmultinom(1, region_counts[y], eval(res))),
                      res = res)
      }
      dt_region <- dt_region[, !(colnames(dt_region) %in% c("AA", "genes", "position", "proportion")), with = FALSE]
      dt_region[, score := unlist(cov, use.names = FALSE)]
      if (libClass == "GRanges") {
        if (length(read_lengths) == 1) {
          dt_region[, size := sample(rep(read_lengths, .N), size = .N, replace = TRUE)]
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

add_sequence_bias <- function(simGenome, dt_range, tAI, region_ranges, lengths, region) {
  if (!is.null(tAI$variable)) {
    tAI <- tAI[variable == unique(variable)[1],]
    stopifnot(nrow(tAI) > 0)
    tAI$variable <- NULL
  }
  if (is.factor(tAI$seqs)) tAI[, seqs := as.character(seqs)]
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
  dt_range[, position := seq_len(.N), by = genes]
  tAI_short <- tAI[, c("seqs", "alpha")]
  tAI_merged <- data.table::merge.data.table(data.table(seqs = unlist(seqs, use.names = FALSE)),
                                             tAI_short, by = "seqs", sort = FALSE)
  # Create the alpha matrix for region
  seq_alpha <- tAI_merged$alpha[!is.na(tAI_merged$alpha)]
  seq_alpha <- split(seq_alpha, dt_range$genes[c(T, F, F)])
  # seq_lengths <- lengths / 3
  # if (any(seq_lengths != as.integer(lengths/3))) stop("Mismatch of seqlength and divisor")
  # seq_length_max <- max(seq_lengths)
  # seq_alpha_max <- lapply(seq_alpha, function(x)
  #   c(x, rep(1e-24, length.out = seq_length_max - length(x))))
  # alpha_matrix <- t(matrix(unlist(seq_alpha_max, F, F), nrow = seq_length_max))
  seq_alpha
}


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
  if (!is.null(seq_acf)) { # Higher order auto correlation
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
  if (!is.null(alpha_matrix)) { # Codon to NT level
    res <- lapply(res, function(x) {
      codon_ac_rnase_alphas_ideal <- rep(x, each = 3)
      codon_ac_rnase_alphas_ideal[seq(codon_ac_rnase_alphas_ideal) %% 3 %in% c(2,0)] = 0
      codon_ac_rnase_alphas_ideal
    })
  }


  if (!is.null(rnase_acf)) { # Lower order auto correlation
    rnase2 <- rep(0, length(rnase_acf) - 1)
    res <- lapply(res, function(x) {
      x <- c(rnase2, x, rnase2)
      frollapply(x, FUN = function(i) sum(i*rnase_acf),
              N = length(rnase_acf), align = "center", fill = NA)
    })
  }
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
#' Should be Direclet alpha dispersion value tables per sequence motif.
#' Stored in files called
#' \code{paste0(type, "_bias_", shift, "_estimates_human.csv")} in the
#' 'dir' folder.
#' @param type character, default: "AA". Sequence motif type, alternatives:
#' codon
#' @param shift character, default "p-site". Alternative: "a-site"
#' @param dir Directory with sequence biases, default is internal path
#' predefined estimators: system.file(package = "coverageSim", "extdata")
#' @param bias character, default: "start_codon". Alternative: "stop_codon",
#' "similar" (uniform like), or "all" (to load all).
#' @return a data.table of bias per sequence motif
#' @export
#' @examples
#' load_seq_bias()
#' load_seq_bias(type = "codon")
#' load_seq_bias(bias = "all")
load_seq_bias <- function(type = "AA", shift = "p-site",
                          dir = system.file(package = "coverageSim", "extdata"),
                          bias = "start_codon") {
  stopifnot(type %in% c("AA", "codon"))
  stopifnot(shift %in% c("p-site", "a-site"))
  shift <- gsub("-", "_", shift)
  file <- paste0(type, "_bias_", shift, "_estimates_human.csv")
  dt <- fread(file.path(dir, file))
  if (bias == "start_codon") {
    dt <- dt[variable == "R2",]
  } else if (bias == "stop_codon") {
    dt <- dt[variable == "R1",]
  } else if (bias == "similar") {
    dt <- dt[variable == "R10",]
  } else if (bias == "all") {
    # Return full
  } else stop("bias must be either of: start_codon, stop_codon or similar!")
  return(dt)
}

append_rnase_to_dt <- function(dt_range, lengths, rnase_bias) {
  gene_split_sites_end <- cumsum(lengths)
  gene_split_sites_start <- c(1, (gene_split_sites_end + 1)[-length(lengths)])
  gene_split_sites <- sort(c(gene_split_sites_start, gene_split_sites_end))
  rnase_reach <- floor(length(rnase_bias[["RFP"]]) / 2)
  gene_split_sites <- rep.int(gene_split_sites, rnase_reach)
  new_index_map <- sort(c(seq.int(nrow(dt_range)), gene_split_sites))
  dt_range <- dt_range[new_index_map, ]
  dt_range[position == 1, start := start - seq.int(3,0)]
  dt_range[dt_range[, .I[position == max(position)], by=genes]$V1, start := start + seq.int(0,3)]
  dt_range[, end := start]
  return(dt_range)
}

#' Get shape function
#'
#' Either uniform 0 (null model), cosine signals 1-4 or
#' rolling windows with auto correlation max.lag 'i'
#' @param i integer index of default shapes to get,
#'  default is 3. Which returns a a function using 3 cosings.
#' @param Z integer, default 0.1. For cosine shapes,
#'  the zero inflation scaler,  the avoid 0 probability positions
#' @param v numeric, default 1. For cosine shapes, scales the first
#' order cosine by the power of 'v'
#' @return a quote object to be evaluated inside the coverage function
#' @export
shapes <- function(i = 3, Z = 0.1, v = 1) {
  if (i == 0) {
    res <- NULL
  } else if (i == 1) {
    #quote(abs(cos(seq(0, pi*(x/30), pi/30))) + 0.1)
    res <- bquote(abs(cos(x))^.(v) + .(Z))
  } else if (i == 2) {
    res <- bquote((abs(cos(x)) + abs(cos(x*b[i])))^.(v) + .(Z))
  } else if (i == 3) {
    res <- bquote((abs(cos(x)) + abs(cos(x*b[i])) + abs(cos(x*(b[i]^2))))^.(v) + .(Z))
  } else if (i == 4) {
    res <-bquote(abs(cos(seq(0, pi*(l/30), pi/30))) + .(Z))
  } else {
    res <- bquote(autocor_window(alpha_vec, .(i), padding.rm = T))
  }
  return(res)
}

reset_matrix <- function(res_matrix, alpha_matrix, lengths) {
  alpha_mat_3 <- alpha_matrix[, rep(seq(ncol(alpha_matrix)), each=3)]
  flank <- (ncol(res_matrix) - ncol(alpha_mat_3))/2
  if (flank < 0) stop("flank < 0, report this error on github!")
  if (flank > 0) {
    alpha_mat_3 <- cbind(alpha_mat_3[, rep(1, flank)],
                         alpha_mat_3, alpha_mat_3[, rep(ncol(alpha_mat_3), flank)])
    all_lengths <- rep(lengths, each = flank) + seq(flank) + flank
    for (i in seq(nrow(alpha_mat_3))) {
      index <- i*3 -3 + 1
      alpha_mat_3[i, all_lengths[seq(index, index+2)]] <- rep(1, flank)
    }
  }
  return(alpha_mat_3)
}

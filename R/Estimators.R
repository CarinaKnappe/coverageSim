
#' Ribo seq frame estimator
#' @param dt data.table of coverage of counts per column
#' @param normalize_to numeric, default 5.
#' @param only_used_codons logical, default TRUE. If FALSE, use all codons.
#' @param relative logical, default TRUE. Relative usage, if FALSE total usage.
#' @return data.table of estimators
#' @export
frame_usage <- function(dt, normalize_to = 5, only_used_codons = TRUE, relative = TRUE) {
  all_frame_usage <- suppressWarnings(melt(dt))
  all_frame_usage[, frame := as.factor(rep.int(seq.int(3L), nrow(all_frame_usage) / 3))]
  all_frame_usage[, codon_sum := frollsum(x = value, n = 3, align = "left")]
  all_frame_usage[, codon_sum := rep(codon_sum[frame == 1], each = 3)]

  if (only_used_codons) {
    all_frame_usage <- all_frame_usage[codon_sum > 0,]
  }
  if (relative) {
    all_frame_usage[, value := (value / codon_sum)]
  }

  all_frame_usage <- all_frame_usage[, .(sum = sum(value, na.rm = TRUE), var = var(value, na.rm = TRUE),
                                         N = .N, g_zero = sum(value > 0)),
                                     by = .(variable, frame)]

  all_frame_usage[, `:=`(mean = sum/N, sd = sqrt(var))]
  all_frame_usage[, alpha := dirichlet_params(mean, sd), by = variable]
  all_frame_usage[, dispersion := (mean^2) / (var - mean)]
  all_frame_usage[, g_zero_rel := round((g_zero/N)*normalize_to, 1), by = .(variable)]
  all_frame_usage[, frame_usage := round((sum/max(sum))*normalize_to, 1), by = .(variable)]
  all_frame_usage[, frame_usage_int := round((sum/max(sum))*normalize_to), by = .(variable)]
  #all_frame_usage[, median(frame_usage), by = frame]
  all_frame_usage[]
  return(all_frame_usage)
}

# dierlich MOM estimater
dirichlet_params <- function(p.mean, sigma){
  n.params <- length(p.mean)
  if(n.params != length(sigma)){
    stop("Length of mean different from length of sigma")
  }
  # Compute second moment
  p.2 <- sigma^2 + p.mean^2
  # Initialize alpa vector
  alpha <- numeric(n.params)
  for (i in 1:(n.params-1)){
    alpha[i] <- (p.mean[1] - p.2[1])*p.mean[i]/(p.2[1] - p.mean[1]^2)
  }
  alpha[n.params] <- (p.mean[1] - p.2[1])*(1-sum(p.mean[-n.params]))/(p.2[1] - p.mean[1]^2)
  return(alpha)
}

dispersion_test <- function(x, round = 5, silent = TRUE)
{
  res <- 1-2 * abs((1 - pchisq((sum((x - mean(x))^2)/mean(x)), length(x) - 1))-0.5)
  if (!silent) {
    cat("Dispersion test of count data:\n",
        length(x), " data points.\n",
        "Mean: ",mean(x),"\n",
        "Variance: ",var(x),"\n",
        "Probability of being drawn from Poisson distribution: ",
        round(res, round),"\n", sep = "")
  }
  invisible(res)
}


GENETIC_CODE_ORFik <- function(as.dt = FALSE, with.charge = FALSE,
                               code = Biostrings::GENETIC_CODE) {
  x <- Biostrings::GENETIC_CODE
  stops <- names(x)[x == "*"]
  x_stop <- x[(names(x) %in% stops)]

  x <- c("###" = "#", "%%%" = "%", x[!(names(x) %in% stops)], "&&&" = "&", x_stop)
  if (as.dt) x <- data.table(AA = x, codon = names(x))
  return(x)
}

## Amino Acid content check
translate_orf_seq <- function(cds, faFile, is.sorted = TRUE,
                              as = "AA", start.as.hash = FALSE,
                              stopm1.as.amp = FALSE, startp1.as.per = FALSE,
                              return.as.list = FALSE,
                              genetic.code = GENETIC_CODE) {
  stopifnot(all(widthPerGroup(cds) %% 3 == 0))
  stopifnot(as %in% c("AA", "codon"))
  seqs <- txSeqsFromFa(cds, faFile, is.sorted = TRUE)
  if (as == "AA") {
    hash <- "#"; amp <- "&"; per <- "%"
    end <- end_amp <- start <- 1; m_width <- 2; ms_width <- 3
    seqs <- as.character(translate(seqs, genetic.code = genetic.code))
  } else {
    hash <- "###"; amp <- "&&&"; per <- "%%%"
    end <- end_amp <- 3; start <- 5; m_width <- 6; ms_width <- 9
    seqs <- as.character(seqs)
  }

  seq_width <- nchar(seqs)

  if (start.as.hash) substring(seqs, 1, end) <- hash
  if (stopm1.as.amp) {
    lt2 <- seq_width > m_width
    if (any(lt2)) {
      substring(seqs[lt2], seq_width[lt2] - start, seq_width[lt2] - end_amp) <- amp
    }
  }
  if (startp1.as.per) {
    lt3 <- seq_width > ms_width
    if (any(lt3)) {
      substring(seqs[lt3], end + 1, start + 1) <- per
    }
  }

  if (as == "codon"){
    #subseq(seqs[lt2], width(seqs[lt2]) - 2, width(seqs[lt2])) <- "***"
    seqs <- lapply(seqs, function(seq) {
      stringr::str_sub(string = seq,
                       start = seq(1, nchar(seq) - 2, by = 3),
                       end = seq(3, nchar(seq), by = 3))
    })
    if (!return.as.list) seqs <- unlist(seqs, use.names = FALSE)
  }
  if (return.as.list && as == "AA") {
    seqs <- unlist(strsplit(seqs, split = ""))
  }
  return(seqs)
}

seq_usage <- function(dt, seqs, genes, input.dt.length = 1, output.seq.length = 3,
                      seqs.order.table = NULL, dispersion_method = "MME") {
  stopifnot(is(dt, "data.table"))
  if ("genes" %in% colnames(dt)) {
    stop("dt can not contain column called 'genes'",
         "it should be given seperatly in argument 'genes'!")
  }
  if (length(genes) != nrow(dt)) stop("Length of dt and genes must match!")
  if (length(seqs) != (nrow(dt) / output.seq.length))
    stop("Size of 'seqs' must be nrow(dt) / output.seq.length")
  codon_sums <- copy(dt)
  n.samples <- ncol(dt)
  if (n.samples > 1) {
    codon_sums <- suppressWarnings(melt.data.table(codon_sums,
                                                   value.name = "score"))
  } else {
    stopifnot("score" %in% colnames(codon_sums))
    codon_sums[, variable := "lib1"]
  }
  codon_sums[, genes := rep(genes, length.out = .N)]
  we_must_collapse <- input.dt.length != output.seq.length
  if (we_must_collapse) {
    codon_sums[, codon_sum := frollsum(x = score, n = output.seq.length, align = "left")]
    # Keep only first per seq group
    codon_sums <- codon_sums[rep(c(T, rep(FALSE, output.seq.length - 1)), length.out = .N),]
  }

  codon_sums[, seqs := rep(seqs, length.out = .N)]

  # Normalize
  codon_sums[, `:=`(gene_sum = sum(score, na.rm = TRUE)), by = .(variable, genes)]
  codon_sums[, `:=`(N_AA_of_type_per_gene = .N), by = .(variable, genes, seqs)]
  codon_sums[, `:=`(as_prob_normalized = score / gene_sum / N_AA_of_type_per_gene )]
  codon_sums[, `:=`(as_prob_normalized = as_prob_normalized / sum(as_prob_normalized)),
             by = .(variable, genes)]
  codon_sums <- codon_sums[, .(score = sum(score), as_prob_normalized = sum(as_prob_normalized),
                           N_total = .N), by = .(variable, genes, seqs)]
  seq_scores <- codon_sums[, .(sum = sum(score, na.rm = TRUE),
                              sum_txNorm = sum(as_prob_normalized, na.rm = TRUE),
                              var = var(score, na.rm = TRUE),
                              var_txNorm = var(as_prob_normalized, na.rm = TRUE),
                              N = .N, N.total = sum(N_total)), by = .(variable, seqs)]

  seq_scores[, mean := sum / N]
  seq_scores[!is.finite(mean), mean := 0]
  seq_scores[, mean_txNorm := sum_txNorm / N]
  seq_scores[!is.finite(mean_txNorm), mean_txNorm := 0]
  if (dispersion_method == "MLE") {
    stop("MLE dispersion not implemented yet!")
  } else if (dispersion_method == "MME") {
    seq_scores[, dispersion := (mean^2) / (var - mean)]
    seq_scores[!is.finite(dispersion) | dispersion < 0, dispersion := 0.001]
    seq_scores[, dispersion_txNorm := (mean_txNorm^2) / (var_txNorm - mean_txNorm)]
    seq_scores[!is.finite(dispersion_txNorm) | dispersion_txNorm < 0, dispersion_txNorm := 0.001]
  }
  seq_scores[, mean_percentage := (mean / sum(mean))*100, by = variable]
  seq_scores[, mean_txNorm_prob := (mean_txNorm / sum(mean_txNorm)), by = variable]
  seq_scores[, mean_txNorm_percentage := mean_txNorm_prob*100]
  seq_scores[, alpha := dirichlet_params(mean_txNorm_prob, sqrt(var_txNorm)),
             by = variable]
  if (is.null(seqs.order.table)) {
    setorderv(seq_scores, c("variable","mean_txNorm_percentage"), order = c(1,-1))
  } else {
    seq_scores[, seqs := factor(seqs, levels = seqs.order.table, ordered = TRUE)]
    setorderv(seq_scores, c("variable","seqs"), order = c(1,1))
  }

  seq_scores[]
  return(seq_scores)
}

AA_score <- function(grl, reads, seqs, weight = "score", is.sorted = TRUE,
                     algorithm = ifelse(all(widthPerGroup(grl) %% 3 == 0), "fast", "exact"),
                     dispersion_method = "MME") {
  stopifnot(dispersion_method %in% c("MME", "MLE"))
  dt <- coveragePerTiling(grl, reads, as.data.table = T, weight = weight, is.sorted = is.sorted)
  return(seq_usage(dt[, .(score = count)], seqs, dt[, genes], dispersion_method))
}

#' Find proportion used of uORFs vs gene regions
#' @param gene character name of gene
#' @param uORFs GRangesList of uORF genomic coordinates
#' @param merged_gr the coverage track as GRanges
#' @export
#' @return a list of length 2: data.table of gene tiling coverage and sum coverage.
overlap_props <- function(gene, uORFs, merged_gr) {
  mrna <- loadRegion(df, "mrna", names.keep = gene)
  count_mRNA <- countOverlapsW(mrna, merged_gr, "score")
  count_leader <- countOverlapsW(loadRegion(df, "leaders", names.keep = gene), merged_gr, "score")
  count_cds <- countOverlapsW(loadRegion(df, "cds", names.keep = gene), merged_gr, "score")
  count_uORFs <- countOverlapsW(uORFs, merged_gr, "score")
  dt_gene <- data.table(counts = c(count_mRNA, count_leader, count_cds, count_uORFs),
                        region = c("MRNA", "LEADERS", "CDS", paste0("UORF", seq(length(uORFs)))))
  count_table <- matrix(count_mRNA)

  colData <- DataFrame(libtype = as.factor("RFP"),
                       condition = as.factor("WT"),
                       replicate = as.factor(1))
  colData$SAMPLE <- paste(colData$libtype, colData$condition, colData$replicate, sep = "_")
  colnames(count_table) <- colData$SAMPLE
  rownames(colData) <- colnames(count_table)
  count_table <- SummarizedExperiment(assays = count_table, rowRanges = mrna, colData = colData)
  return(list(dt_gene, count_table))
}

coverage_all_cds_all_samples <- function(df, cds, prefix, lib.type = "pshifted") {
  if (!is.character(prefix)) stop("Prefix must be character vector!")
  stopifnot(length(cds) > 0)
  message(name(df))
  message("-- Number of transcripts: ", length(cds))
  dt <- data.table()
  for (i in seq(nrow(df))) {
    file <- filepath(df, type = lib.type)[i]
    print(basename(file))
    dt_samp <- coveragePerTiling(cds, reads = fimport(file),
                                 is.sorted = TRUE, as.data.table = T)
    dt <- cbind(dt, dt_samp$count)
  }
  colnames(dt) <- paste0(prefix, seq(1, ncol(dt)))
  return(dt)
}

# All vs all comparison
coverage_cor <- function(..., method = "spearman", na.rm = TRUE,
                         rm.selfmatch = TRUE, is.listed = FALSE) {
  stopifnot(is.logical(na.rm))
  if (is.listed) {
    dt <- as.data.table(...)
  } else dt <- cbind(...)
  print("-- Fast calculating correlation..")
  dt_melt <- suppressWarnings(melt(get_upper_tri(round(cor(dt, method = method), 2)),
                                   variable.name = c("Sample1", "Sample2"), value.name = "Cor"))
  print("-- Cor calculated")
  setDT(dt_melt)
  if (na.rm) dt_melt <- dt_melt[!is.na(Cor),]
  if (rm.selfmatch) dt_melt <- dt_melt[Var1 != Var2,]
  dt_melt[, Var1_index := as.integer(substr(Var1, 2, 3))]
  dt_melt[, Var2_index := as.integer(substr(Var2, 2, 3))]
  dt_melt[, Var1_type := substr(Var1, 1, 1)]
  dt_melt[, Var2_type := substr(Var2, 1, 1)]
  dt_melt[, Comparison := paste(Var1_type, "v", Var2_type)]
  dt_melt[]
  return(dt_melt)
}

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

cor_upper_tri <- function(dt, method = c("pearson", "spearman")[1],
                          decimals = 2, melt = TRUE) {
  cor <- round(cor(dt, method = method), decimals)
  cor <- get_upper_tri(dt)
}

auto_correlation <- function(dt, dist = 6, by.codon = TRUE, codon.vs.nt = FALSE,
                             method = "spearman", fast.cor = TRUE) {
  if (dist > 50) stop("Max dist allowed is 50!")
  dt_per_codon <-
  if (by.codon) {
    if (codon.vs.nt) {
      dt[rep(seq(1, .N, by = 3), each = 3),]
    } else dt[(seq(nrow(dt))-1) %% 3 == 0,]
  } else dt
  # Set names to R
  colnames(dt_per_codon) <- paste0("R", seq(ncol(dt_per_codon)))
  codons_dist <- dist
  selection <- seq(codons_dist)
  # Remove F, to avoid match with logical F
  if (dist > 17) selection <- c(selection[-18], codons_dist + 1)
  if (dist > 5) selection <- c(selection[-6], max(selection) + 1)
  names_codon_dist <- c(LETTERS, letters)[selection]
  shift_list <- list()
  print("-- Calculating auto correlation")
  for (c in seq(codons_dist)) {
    to_r <- seq(c)
    dt_dist_shifts <-
    if (codon.vs.nt) {
      rbind(dt[-to_r,], dt[rep(1, c)*nrow(dt),])
    } else rbind(dt_per_codon[-to_r,], dt_per_codon[rep(1, c)*nrow(dt_per_codon),])

    colnames(dt_dist_shifts) <- paste0(names_codon_dist[c], seq(ncol(dt_dist_shifts)))
    shift_list <- c(shift_list, dt_dist_shifts)
  }
  distance_index <- rep(seq(codons_dist), each = ncol(dt_per_codon))
  all_cor <-
  if (fast.cor) {
    cor_temp <- coverage_cor(cbind(dt_per_codon, setDT(shift_list)), method = method,
                             na.rm = T, rm.selfmatch = T)
    cor_temp[Var1_type != Var2_type & Var1_type == "R" & Var1_index == Var2_index,]$Cor
  } else {
    unlist(lapply(unique(distance_index), function(i) {
      cat(i, ",")
      round(diag(cor(dt_per_codon, setDT(shift_list[distance_index == i]),
                     method = method)),
            2)
    }))
  }
  dt_all_c_c_comp <- data.table(Cor = all_cor)
  dt_all_c_c_comp[, distance := as.character(distance_index)]
  return(dt_all_c_c_comp)
}

auto_correlation_fast <- function(dt, dist = 6, by.codon = TRUE, codon.vs.nt = FALSE,
                                  genes = NULL, fun = acf) {
  dt <-
  if (by.codon) {
    if (codon.vs.nt) {
      dt[rep(seq(1, .N, by = 3), each = 3),]
    } else dt[(seq(nrow(dt))-1) %% 3 == 0,]
  } else dt

  ress <- lapply(dt, function(x) fun(x, lag.max = dist, plot = FALSE)$acf[-1])
  ress_table <- data.table(Cor = unlist(ress),
                           distance = rep(seq(1, dist), length.out = length(unlist(ress))))
  return(ress_table)
}

#' Get gene auto correlation
#' @param dt a data.table of counts
#' @param dist numeric, default 6. Distance in nt or codons to check
#' @param by.codon logical, default TRUE. Else by nt
#' @param codon.vs.nt logical, default FALSE Else convert codon to nt space
#' (1,0,0)
#' @param genes numeric vector, genes, index of which gene this is from
#' @param fun which auto correlation function, default stats::acf
#' @param mean logical, default FALSE Else get only mean acf per position.
#' @return a data.table of counts
#' @export
auto_correlation_genes <- function(dt, dist = 6, by.codon = TRUE, codon.vs.nt = FALSE,
                                   genes, fun = acf, mean = FALSE) {
  iterator <- seq_along(dt)
  if (!is.null(genes)) dt[, genes := genes]
  dt <-
    if (by.codon) {
      if (codon.vs.nt) {
        dt[rep(seq(1, .N, by = 3), each = 3),]
      } else dt[(seq(nrow(dt))-1) %% 3 == 0,]
    } else dt

  res <- melt.data.table(dt, id.vars = "genes", variable.name = "sample", value.name = "count", variable.factor = T)
  res[, sample := as.integer(sample)]
  is_acf <- identical(fun, acf)
  res <- if (is_acf) {
    res[, .(Cor = fun(count, lag.max = dist, plot = FALSE)$acf[-1]), by = .(sample, genes)]
  } else res[, .(Cor = as.vector(fun(count, lag.max = dist, plot = FALSE)$acf)), by = .(sample, genes)]

  res[, distance := seq.int(.N), by = .(sample, genes)]
  if (mean == TRUE) res <- res[, .(Cor = mean(Cor, na.rm = T)), by = distance]
  return(res)
}

auto_correlation_genes_all <- function(dt_all_list_named, dist = 100, by.codon = T, genes,
                                       codon.vs.nt = T, fun = acf, mean = F) {
  ac_list <- lapply(dt_all_list_named, function(x)
    auto_correlation_genes(x, dist = dist, by.codon = by.codon, genes = genes,
                           codon.vs.nt = codon.vs.nt, fun = fun, mean = mean))
  return(ac_list_to_dt(ac_list))
}

#' Auto correlation
#' @param vec vector of coverage
#' @param max.lag integer, the max lag for correlation window
#' @param fill NA
#' @param na.rm logical, FALSE
#' @param padding.rm logical or integer, if integer, keep this amount
#' of padding.
#' @param penalty 1.5, for auto correlation window, the power scaler for
#' penalty. Higher value, lower correlation far away.
#' @return the input vector with applied rolling function
autocor_window <- function(vec, max.lag, fill = NA, na.rm = FALSE,
                           padding.rm = FALSE, penalty = 1.5) {
  window <- max.lag*2 + 1
  split_2 <- ceiling(window/2)
  split_2_low <- floor(window/2)
  roll_function <- c(seq(split_2 , 2), 1, seq(2, split_2))^penalty
  padding_left <- padding_right <- rep(0, split_2_low*2)
  roll <- frollapply_compat(c(padding_left, vec, padding_right),
                            window = window,
                            FUN = function(x) mean(x/roll_function, na.rm =T),
                            align = "center", fill = NA)
  if (padding.rm) {
    padd_to_keep <- ifelse(padding.rm > 1, padding.rm, 0)
    left_pad_index <-  seq_along(padding_left - padd_to_keep)
    right_pad_index <-  seq(length(roll) - length(padding_right) +
                              1 + padd_to_keep, length(roll))
    roll <- roll[-c(left_pad_index, right_pad_index)]
  } else if (na.rm) roll <- roll[!is.na(roll)]
  return(roll)
}

ac_list_to_dt <- function(ac_list) {
  dt <- rbindlist(ac_list)
  dt[, id := rep(names(ac_list), each = nrow(ac_list[[1]]))]
  dt[, distance := factor(distance, levels = unique(distance), ordered = TRUE)]
  dt[, col := factor(as.integer(distance) %% 3)]
  if (!is.null(dt$genes)) dt[, mean_Cor := mean(Cor, na.rm = T), by = genes]
  dt <- dt[!is.na(Cor),]
  dt[]
  return(dt)
}

#' Positional boxplot of auto correlation
#' @param ac_dt a data.table
#' @param breaks.by 9 (x-axis breaks, default 9 (9 codons = 1 ribosome))
#' @param autocor_name "codon", else "nt"
#' @param plot logical, default TRUE, If FALSE, don't plot, only return
#' @return a ggplot object
#' @export
ac_boxplot <- function(ac_dt, breaks.by = 9, autocor_name = "codon", plot = T) {
  gg <- ggplot(ac_dt, aes(y = Cor, x = distance, fill = col)) + geom_boxplot(outlier.size = 0.1, ) +
    xlab(paste("Upstream", autocor_name)) + theme_classic() + facet_wrap(~ id) +
    scale_x_discrete(breaks = breaks.by) + geom_vline(xintercept =
      seq(breaks.by, nrow(ac_dt), by = breaks.by), colour = "gray", size = 0.1, linetype = "dashed")
  if (plot) plot(gg)
  gg
}

period_detector <- function(acf, quantile_min = 0.8, xlims = c(1, 25), line_center = 9) {
  # mean.filter = quantile(acf$mean_Cor, quantile_min)
  # merge_dt <- acf[mean_Cor > mean.filter,]
  merge_dt <- acf[acf[, .I[mean_Cor > quantile(mean_Cor, quantile_min)], by = sample]$V1]
  merge_dt <- merge_dt[, .(Cor = mean(Cor, na.rm = T)), by = .(sample, distance, id)]
  merge_dt <- merge_dt[, (spec.pgram(x = Cor, plot = F)[c(1,2)]), by = sample]
  ggplot(merge_dt, aes(1 / freq, spec)) + geom_line() + theme_classic() + facet_wrap(~ sample, scales = "free_y") +
    xlim(xlims) + geom_vline(xintercept = line_center, colour = "gray", size = 0.1, linetype = "dashed")
}

skewness <- function (x, na.rm = FALSE) {
    if (is.matrix(x))
      apply(x, 2, skewness, na.rm = na.rm)
    else if (is.vector(x)) {
      if (na.rm) x <- x[!is.na(x)]
      n <- length(x)
      (sum((x-mean(x))^3)/n)/(sum((x-mean(x))^2)/n)^(3/2)
    }
    else if (is.data.frame(x))
      sapply(x, skewness, na.rm = na.rm)
    else skewness(as.vector(x), na.rm = na.rm)
}

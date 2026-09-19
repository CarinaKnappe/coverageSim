
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

# NOTE: overlap_props() was removed here (2026-09-18). It referenced an
# undefined global `df` (an ORFik experiment object was clearly intended,
# but never passed as a parameter), so it errored on every call unless a
# `df` happened to exist by coincidence in the calling environment. It was
# unused anywhere else in this package, in the test suite, or in the
# rust_helpers/choros scripts, and was pre-existing article/figure code
# (not part of the simulation pipeline) rather than a Carina change. See
# claude/code-review-carina-dev-2026-09.md in the CoverageSimulator Claude
# project for details. Run `devtools::document()` to regenerate NAMESPACE
# and remove man/overlap_props.Rd accordingly.

# NOTE (2026-09-19): Estimators.R was split into "pipeline-relevant" vs.
# "article/figure code" (see claude/code-review-carina-dev-2026-09.md,
# section "Estimators.R"). Removed here as unused, unexported, and with no
# remaining callers in this package, the test suite, or the choros/rust
# scripts: dispersion_test(), GENETIC_CODE_ORFik() (note: despite the
# similar name it was never actually used by translate_orf_seq() above --
# that function's default genetic.code argument resolves to
# Biostrings::GENETIC_CODE, not this one), seq_usage(), AA_score() (its only
# caller), coverage_all_cds_all_samples(), get_lower_tri(), cor_upper_tri()
# (itself dead/buggy code -- computed `cor` but never used it),
# auto_correlation_genes_all() and its only caller-less helper
# ac_list_to_dt(), period_detector(), and skewness(). Also removed in a
# second pass the same day: auto_correlation(), auto_correlation_fast(), and
# their only dependencies coverage_cor() and get_upper_tri() -- none of the
# four were exported or called anywhere in the package/tests/choros/rust
# scripts. They are unrelated to the actual autocorrelation bias model used
# by the simulation pipeline, which lives entirely in
# Sim_reads_from_counts_helpers.R / learn_end_bias.R / Sim_reads_from_counts.R
# as the `auto_correlation` parameter/kernel -- a different code path,
# untouched by this. None of the removed functions had roxygen documentation
# or @export, so no NAMESPACE/man/ regeneration is needed. Full removed code
# is preserved in git history (`git log -p -- R/Estimators.R`) if any of it
# is needed later for the separate article-code repo.

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

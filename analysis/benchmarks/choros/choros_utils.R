library(data.table)
library(stringi)
library(ggplot2)

resolve_covsim_config <- function(base_dir, output_dir = "", run_id = "RFP_WT_1") {
  if (!nzchar(base_dir)) {
    stop("Set CHOROS_COVSIM_BASE to a coverageSim dataset directory.")
  }
  base_dir <- normalizePath(path.expand(base_dir), mustWork = TRUE)
  if (!nzchar(output_dir)) output_dir <- file.path(base_dir, "choros_output")
  list(
    base_dir = base_dir,
    output_dir = path.expand(output_dir),
    run_id = run_id,
    bam_file = file.path(base_dir, "reads", paste0(run_id, ".bam")),
    fasta_file = file.path(base_dir, "genome", "human_flavoured_sim.fasta"),
    txdb_file = file.path(base_dir, "genome", "human_flavoured_sim.gtf.db")
  )
}

set_unknown_circular_to_false <- function(x) {
  info <- GenomeInfoDb::seqinfo(x)
  circular <- GenomeInfoDb::isCircular(info)
  circular[] <- FALSE
  GenomeInfoDb::isCircular(info) <- circular
  GenomeInfoDb::seqinfo(x) <- info
  x
}


align_regions_by_name <- function(regions, transcripts, keep_names) {
  keep <- keep_names[
    keep_names %in% names(regions) & keep_names %in% names(transcripts)
  ]
  regions <- regions[keep]
  transcripts <- transcripts[keep]
  if (!identical(names(regions), names(transcripts))) {
    stop("CDS and transcript regions could not be aligned by name.")
  }
  list(regions = regions, transcripts = transcripts)
}

map_read_ends_to_transcripts <- function(reads, transcripts) {
  shared <- intersect(
    GenomeInfoDb::seqlevels(reads), GenomeInfoDb::seqlevels(transcripts)
  )
  reads <- GenomeInfoDb::keepSeqlevels(
    reads, shared, pruning.mode = "coarse"
  )
  mapped <- GenomicFeatures::mapToTranscripts(reads, transcripts)
  if (!length(mapped)) stop("No read/transcript overlaps found.")
  read_index <- S4Vectors::mcols(mapped)$xHits
  transcript_index <- S4Vectors::mcols(mapped)$transcriptsHits
  list(
    mapped = mapped, reads = reads[read_index],
    transcript_index = transcript_index
  )
}

select_periodic_lengths <- function(mapped_reads, candidate_lengths = 25:30,
                                    min_reads = 1000L, min_frame_prop = 0.5) {
  required <- c("L", "dist_start", "score")
  missing <- setdiff(required, names(mapped_reads))
  if (length(missing)) stop("mapped_reads is missing: ", paste(missing, collapse = ", "))
  qc <- mapped_reads[L %in% candidate_lengths,
    .(count = sum(score)), by = .(L, frame = dist_start %% 3L)]
  if (!nrow(qc)) stop("No reads have a candidate RPF length.")
  qc[, total := sum(count), by = L]
  qc[, frame_prop := count / total]
  qc[, max_frame_prop := max(frame_prop), by = L]
  valid <- qc[total >= min_reads & max_frame_prop >= min_frame_prop,
              sort(unique(L))]
  if (!length(valid)) stop("No candidate read length passes the periodicity filters.")
  list(lengths = valid, qc = qc[order(L, frame)])
}

infer_tis_offsets <- function(mapped_reads, read_lengths,
                              tis_window = c(-40L, -3L), min_peak_reads = 50L) {
  required <- c("L", "dist_start", "score")
  missing <- setdiff(required, names(mapped_reads))
  if (length(missing)) stop("mapped_reads is missing: ", paste(missing, collapse = ", "))
  profile <- mapped_reads[
    L %in% read_lengths & dist_start >= tis_window[1] & dist_start <= tis_window[2],
    .(count = sum(score)), by = .(L, dist_start)]
  if (!nrow(profile)) stop("No reads occur in the requested TIS window.")
  setorder(profile, L, -count, dist_start)
  peaks <- profile[, .SD[1L], by = L]
  peaks[, d5_base := -dist_start + 3L]
  peaks[, valid := count >= min_peak_reads & d5_base >= 0L & d5_base <= L - 3L]
  failed <- sort(unique(c(setdiff(read_lengths, peaks$L), peaks[valid == FALSE, L])))
  if (length(failed)) stop("No reliable, valid TIS offset for length(s): ",
                           paste(failed, collapse = ", "))
  list(offsets = peaks[order(L)], profile = profile[order(L, dist_start)])
}

build_frame_offset_map <- function(offsets) {
  result <- rbindlist(lapply(seq_len(nrow(offsets)), function(i) {
    read_width <- as.integer(offsets$L[i])
    d5_base <- -as.integer(offsets$dist_start[i]) + 3L
    data.table(qwidth = read_width, pos_frame = 0:2,
      d5 = vapply(0:2, function(frame) {
        candidates <- (d5_base - 2L):(d5_base + 2L)
        valid <- candidates[candidates %% 3L == (3L - frame) %% 3L &
                            candidates >= 0L & candidates <= read_width - 3L]
        if (!length(valid)) stop("Cannot construct offsets for length ", read_width)
        valid[which.min(abs(valid - d5_base))]
      }, integer(1)))
  }))
  unique(result[order(qwidth, pos_frame)])
}

detect_ribosome_shifts_compat <- function(
    footprints, txdb, tx_names, tx, accepted_lengths = 1:100,
    top_tx = 10L, first_n = 150L, min_reads = 1000L,
    min_reads_tis = 50L, strict_fft = TRUE, verbose = TRUE
) {
  txdb <- ORFik::loadTxdb(txdb)
  cds <- ORFik::loadRegion(txdb, part = "cds", names.keep = tx_names)
  footprints <- ORFik::fimport(footprints, cds)
  matched <- ORFik:::validSeqlevels(cds, footprints)
  cds <- GenomeInfoDb::keepSeqlevels(cds, matched, pruning.mode = "coarse")
  footprints <- GenomeInfoDb::keepSeqlevels(
    footprints, matched, pruning.mode = "coarse"
  )
  if (!length(cds) || !length(footprints)) {
    stop("TxDb and footprints have no matching sequence names.")
  }

  cds <- set_unknown_circular_to_false(cds)
  footprints <- set_unknown_circular_to_false(footprints)
  tx <- set_unknown_circular_to_false(tx)
  aligned <- align_regions_by_name(cds, tx, tx_names)
  cds <- aligned$regions
  tx <- aligned$transcripts
  footprints <- ORFik:::convertToOneBasedRanges(
    footprints, addSizeColumn = TRUE, addScoreColumn = TRUE,
    along.reference = TRUE
  )
  footprints <- set_unknown_circular_to_false(footprints)

  length_counts <- data.table(
    score = footprints$score, size = footprints$size
  )[, .(counts = sum(score)), by = size]
  length_counts <- length_counts[
    counts >= min_reads & size %in% accepted_lengths
  ]
  if (!nrow(length_counts)) stop("No accepted read lengths have enough reads.")

  cds <- cds[ORFik:::countOverlapsW(cds, footprints, "score") > 0]
  aligned <- align_regions_by_name(cds, tx, names(cds))
  cds <- aligned$regions
  tx <- aligned$transcripts
  periodicity <- ORFik::windowPerReadLength(
    cds, tx, footprints, pShifted = FALSE, upstream = 0,
    downstream = first_n - 1L, zeroPosition = 0,
    scoring = "transcriptNormalized",
    acceptedLengths = length_counts$size,
    drop.zero.dt = TRUE, append.zeroes = TRUE
  )
  periodicity <- periodicity[, .(
    periodic = ORFik:::isPeriodic(
      score, unique(fraction), verbose = verbose, strict.fft = strict_fft
    )
  ), by = fraction]
  valid_lengths <- periodicity[periodic == TRUE, fraction]
  if (!length(valid_lengths)) stop("No periodic read lengths were detected.")
  if (verbose) {
    message("Periodic read lengths: ", paste(valid_lengths, collapse = ", "))
  }

  top_tx <- ORFik:::percentage_to_ratio(top_tx, cds)
  tis <- ORFik::windowPerReadLength(
    cds, tx, footprints, pShifted = FALSE, upstream = 30,
    downstream = 29, acceptedLengths = valid_lengths, scoring = NULL
  )
  if (!"count" %in% names(tis)) {
    if (!"score" %in% names(tis)) {
      stop("ORFik TIS coverage contains neither count nor score.")
    }
    tis[, count := score]
  }
  tis[, sum.count := sum(count), by = genes]
  tis <- tis[sum.count >= stats::quantile(sum.count, top_tx)]
  tis <- ORFik::coverageScorings(tis, scoring = "sum")
  tis[, frac.score := sum(score), by = fraction]
  tis <- tis[frac.score > min_reads_tis]
  if (!nrow(tis)) {
    stop("Not enough reads remain in the TIS region for offset detection.")
  }
  tis[, .(
    offsets_start = ORFik:::changePointAnalysis(
      score, info = unique(fraction), verbose = verbose
    )
  ), by = fraction]
}

compute_gc_dt <- function(dt, tx_seqs) {
  dt[, gc := {
    tid <- as.character(transcript_id)
    rpf_5 <- stringi::stri_sub(tx_seqs[tid], A_start - d5, A_start - 7L)
    rpf_3 <- stringi::stri_sub(tx_seqs[tid], A_start + 3L, A_start + 2L + d3)
    rpf <- paste0(rpf_5, rpf_3)
    (stringi::stri_count_fixed(rpf, "G") +
       stringi::stri_count_fixed(rpf, "C")) / (d5 + d3 + 3L)
  }]
  invisible(dt)
}

build_model_dt <- function(dt_obs, tx_seqs, off_long, excl_codons = 20L) {
  required <- c("transcript_id", "cod_idx", "d5", "d3",
                "utr5_len", "cds_len", "count")
  missing <- setdiff(required, names(dt_obs))
  if (length(missing)) stop("dt_obs is missing: ", paste(missing, collapse = ", "))
  if (!nrow(dt_obs)) stop("dt_obs is empty.")

  geom_pairs <- unique(off_long[, .(d5, d3 = qwidth - d5 - 3L)])[d3 >= 0L]
  if (!nrow(geom_pairs)) stop("No valid d5/d3 geometry pairs.")

  tx_info <- unique(dt_obs[, .(transcript_id, utr5_len, cds_len)])
  model_list <- lapply(seq_len(nrow(tx_info)), function(i) {
    cds_codons <- tx_info$cds_len[i] %/% 3L
    if (cds_codons <= 2L * excl_codons) return(NULL)
    data.table(
      transcript_id = tx_info$transcript_id[i],
      cod_idx = seq.int(excl_codons + 1L, cds_codons - excl_codons),
      utr5_len = tx_info$utr5_len[i],
      cds_len = tx_info$cds_len[i]
    )
  })
  model_list <- model_list[!vapply(model_list, is.null, logical(1))]
  if (!length(model_list)) stop("No transcripts are long enough for modeling.")

  model_dt <- rbindlist(model_list)
  model_dt[, join_dummy := 1L]
  gp <- copy(geom_pairs)
  gp[, join_dummy := 1L]
  model_dt <- gp[model_dt, on = "join_dummy", allow.cartesian = TRUE]
  model_dt[, join_dummy := NULL]

  model_dt[, A_start := utr5_len + 3L * (cod_idx - 1L) + 1L]
  tid <- as.character(model_dt$transcript_id)
  model_dt[, f5 := stri_sub(tx_seqs[tid], A_start - d5, A_start - d5 + 2L)]
  model_dt[, f3 := stri_sub(tx_seqs[tid], A_start + d3, A_start + 2L + d3)]
  model_dt[, A := stri_sub(tx_seqs[tid], A_start, A_start + 2L)]
  model_dt[, P := stri_sub(tx_seqs[tid], A_start - 3L, A_start - 1L)]
  model_dt[, E := stri_sub(tx_seqs[tid], A_start - 6L, A_start - 4L)]

  model_dt <- dt_obs[
    model_dt,
    on = c("transcript_id", "cod_idx", "d5", "d3"),
    .(
      transcript_id, cod_idx, d5, d3,
      A_start = i.A_start, utr5_len = i.utr5_len, cds_len = i.cds_len,
      f5 = i.f5, f3 = i.f3, A = i.A, P = i.P, E = i.E,
      count = fifelse(is.na(x.count), 0L, as.integer(x.count))
    )
  ]
  model_dt <- model_dt[
    nchar(f5) == 3L & nchar(f3) == 3L &
      nchar(A) == 3L & nchar(P) == 3L & nchar(E) == 3L
  ]
  if (!nrow(model_dt)) stop("The CHOROS model dataset is empty.")
  model_dt
}

fit_choros <- function(model_dt, tx_seqs, n_top_tx = 250L, min_prop = 0.9) {
  tx_density <- model_dt[
    count > 0L, .(density = mean(count)), by = transcript_id
  ][order(-density)]
  if (!nrow(tx_density)) return(NULL)

  keep_tx <- tx_density[seq_len(min(as.integer(n_top_tx), .N)), transcript_id]
  model_dt <- model_dt[transcript_id %in% keep_tx]

  geom_cnt <- model_dt[, .(n = sum(count)), by = .(d5, d3)][order(-n)]
  total <- sum(geom_cnt$n)
  if (!nrow(geom_cnt) || !is.finite(total) || total <= 0) return(NULL)
  geom_cnt[, cum_prop := cumsum(n) / total]
  first_over <- which(geom_cnt$cum_prop >= min_prop)[1]
  if (is.na(first_over)) first_over <- nrow(geom_cnt)
  keep_geom <- geom_cnt[seq_len(first_over), .(d5, d3)]
  model_dt <- model_dt[keep_geom, on = c("d5", "d3"), nomatch = 0L]
  if (!nrow(model_dt) || max(model_dt$count) <= 0) return(NULL)

  compute_gc_dt(model_dt, tx_seqs)
  factor_cols <- c("transcript_id", "A", "P", "E", "f5", "f3", "d5", "d3")
  model_dt[, (factor_cols) := lapply(.SD, as.factor), .SDcols = factor_cols]
  model_dt <- droplevels(model_dt)

  base_terms <- c("A", "P", "E", "gc")
  full_terms <- c("A", "P", "E", "gc")
  if (nlevels(model_dt$d5) > 1L) {
    base_terms <- c(base_terms, "d5")
    full_terms <- c(full_terms, "d5*f5")
  } else {
    full_terms <- c(full_terms, "f5")
  }
  if (nlevels(model_dt$d3) > 1L) {
    base_terms <- c(base_terms, "d3")
    full_terms <- c(full_terms, "d3*f3")
  } else {
    full_terms <- c(full_terms, "f3")
  }

  fml_base <- as.formula(paste(
    "count ~", paste(base_terms, collapse = " + "), "| transcript_id"
  ))
  fml_full <- as.formula(paste(
    "count ~", paste(full_terms, collapse = " + "), "| transcript_id"
  ))

  message("Fitting base model...")
  fit_base <- fixest::fenegbin(fml_base, data = model_dt, nthreads = 0.5, warn = FALSE)
  message("Fitting full sequence-bias model...")
  fit_full <- fixest::fenegbin(fml_full, data = model_dt, nthreads = 0.5, warn = FALSE)
  if (is.null(fit_base) || is.null(fit_full)) return(NULL)

  bic_base <- BIC(fit_base)
  bic_full <- BIC(fit_full)
  message(sprintf("BIC base: %.1f; BIC full: %.1f", bic_base, bic_full))
  fit <- if (bic_base <= bic_full) fit_base else fit_full

  list(
    model_dt = model_dt,
    fit = fit,
    BIC_base = bic_base,
    BIC_full = bic_full
  )
}

correct_bias_choros <- function(model_dt, fit) {
  b <- coef(fit)
  b[is.na(b)] <- 0
  names(b) <- gsub("::", "", names(b))
  get_coef <- function(name) {
    value <- b[name]
    if (!length(value) || is.na(value)) 0 else as.numeric(value)
  }
  correction_table <- function(seq_levels, digest_levels, seq_prefix, digest_prefix) {
    result <- matrix(
      NA_real_, length(seq_levels), length(digest_levels),
      dimnames = list(seq_levels, digest_levels)
    )
    for (sequence_level in seq_levels) {
      for (digest_level in digest_levels) {
        marginal <- get_coef(paste0(seq_prefix, sequence_level))
        interaction <- get_coef(paste0(
          digest_prefix, digest_level, ":", seq_prefix, sequence_level
        ))
        if (interaction == 0) {
          interaction <- get_coef(paste0(
            seq_prefix, sequence_level, ":", digest_prefix, digest_level
          ))
        }
        result[sequence_level, digest_level] <- exp(marginal + interaction)
      }
    }
    result
  }

  corr5 <- correction_table(
    levels(model_dt$f5), levels(model_dt$d5), "f5", "d5"
  )
  corr3 <- correction_table(
    levels(model_dt$f3), levels(model_dt$d3), "f3", "d3"
  )
  model_dt[, cf5 := corr5[cbind(as.character(f5), as.character(d5))]]
  model_dt[, cf3 := corr3[cbind(as.character(f3), as.character(d3))]]
  # Same non-finite guard and [0.2, 5] clamp as apply_choros_bias_correction()
  # applies to the full observation table: without it, this function's
  # bias_after/absolute_bias evaluation metrics would grade a different,
  # unclamped correction than the one actually written to the output counts.
  model_dt[!is.finite(cf5), cf5 := 1]
  model_dt[!is.finite(cf3), cf3 := 1]
  model_dt[, cf5 := pmin(pmax(cf5, 0.2), 5)]
  model_dt[, cf3 := pmin(pmax(cf3, 0.2), 5)]
  model_dt[, corrected := count / cf5 / cf3]

  if (sum(model_dt$corrected, na.rm = TRUE) > 0) {
    b_gc <- get_coef("gc")
    mean_gc <- model_dt[
      , sum(corrected * gc, na.rm = TRUE) / sum(corrected, na.rm = TRUE)
    ]
    model_dt[, corrected := exp(
      log(pmax(corrected, 1e-10)) + b_gc * mean_gc - b_gc * gc
    )]
    model_dt[, corrected := corrected * sum(count) / sum(corrected)]
  }
  invisible(model_dt)
}

eval_bias_norm <- function(model_dt, tx_seqs, col_name,
                           trunc5 = 20L, trunc3 = 20L) {
  if (!col_name %in% names(model_dt)) stop("Column not found: ", col_name)
  codon_offsets <- -6:6
  codon_names <- c(
    "n6", "n5", "n4", "n3", "E", "P", "A_cod",
    "p1", "p2", "p3", "p4", "p5", "p6"
  )
  agg <- copy(model_dt[, .(
    transcript_id, cod_idx, A_start, utr5_len, cds_len,
    value = get(col_name)
  )])
  agg <- agg[, .(count = as.double(sum(value))), by = .(
    transcript_id, cod_idx, A_start, utr5_len, cds_len
  )]
  agg[, cds_codons := cds_len %/% 3L]
  agg <- agg[cod_idx > trunc5 & cod_idx <= cds_codons - trunc3]
  agg <- agg[(A_start - 18L) >= 1L &
               (A_start + 20L) <= utr5_len + cds_len + 20L]
  positive_tx <- agg[, .(total = sum(count)), by = transcript_id][
    total > 0, transcript_id
  ]
  agg <- agg[transcript_id %in% positive_tx]
  if (!nrow(agg)) stop("No observations remain for bias evaluation.")
  agg[, count := count / mean(count), by = transcript_id]

  tid <- as.character(agg$transcript_id)
  for (i in seq_along(codon_offsets)) {
    agg[, (codon_names[i]) := stri_sub(
      tx_seqs[tid],
      A_start + 3L * codon_offsets[i],
      A_start + 3L * codon_offsets[i] + 2L
    )]
  }
  for (column in codon_names) {
    agg <- agg[nchar(get(column)) == 3L]
    agg[, (column) := as.character(get(column))]
  }
  if (!nrow(agg)) stop("No complete codon windows remain for evaluation.")

  fml <- as.formula(paste("count ~", paste(codon_names, collapse = " + ")))
  coefficients <- coef(lm(fml, data = agg))
  vapply(codon_names, function(column) {
    index <- grep(paste0("^", column), names(coefficients))
    sum(coefficients[index]^2, na.rm = TRUE)
  }, numeric(1))
}

compute_absolute_bias_score <- function(model_dt) {
  required <- c("cf5", "cf3", "count")
  if (!all(required %in% names(model_dt))) {
    stop("model_dt must contain cf5, cf3, and count.")
  }
  observed <- model_dt[
    count > 0 & is.finite(cf5) & is.finite(cf3) & cf5 > 0 & cf3 > 0
  ]
  if (!nrow(observed)) stop("No positive observations for bias scoring.")
  weighted_sd <- function(x, weights) {
    mean_weighted <- sum(x * weights) / sum(weights)
    sqrt(sum(weights * (x - mean_weighted)^2) / sum(weights))
  }
  bias5 <- weighted_sd(log2(observed$cf5), observed$count)
  bias3 <- weighted_sd(log2(observed$cf3), observed$count)
  list(
    bias_5p_score = bias5,
    bias_3p_score = bias3,
    total_bias_score = sqrt(bias5^2 + bias3^2)
  )
}

summarize_end_profile <- function(dt, distance_column,
                                  window = c(-60L, 60L)) {
  if (!distance_column %in% names(dt)) {
    stop("Missing distance column: ", distance_column)
  }
  if (length(window) != 2L || any(!is.finite(window)) || window[1] > window[2]) {
    stop("window must contain two ordered, finite positions")
  }
  positions <- data.table(position = seq.int(window[1], window[2]))
  observed <- dt[
    get(distance_column) >= window[1] & get(distance_column) <= window[2],
    .(count = sum(count)), by = .(position = get(distance_column))
  ]
  profile <- observed[positions, on = "position"]
  profile[is.na(count), count := 0]
  profile[]
}

plot_end_profile <- function(dt, distance_column, title,
                             window = c(-60L, 60L)) {
  profile <- summarize_end_profile(dt, distance_column, window)
  ggplot(profile, aes(position, count)) +
    geom_col(fill = "grey45") +
    geom_vline(xintercept = 0L, linetype = "dashed", colour = "black") +
    scale_x_continuous(
      limits = window,
      breaks = seq.int(window[1], window[2], by = 20L),
      expand = expansion(mult = c(0, 0))
    ) +
    theme_classic() +
    labs(title = title, x = "Distance (nt)", y = "Count")
}

plot_diagnostics <- function(dt_obs) {
  dt <- copy(dt_obs)
  dt[, L := d5 + d3 + 3L]
  if (!"A_start" %in% names(dt)) {
    dt[, A_start := utr5_len + 3L * (cod_idx - 1L) + 1L]
  }
  dt[, tx_5p_pos := A_start - d5]
  dt[, pos_frame := (tx_5p_pos - utr5_len - 1L) %% 3L]
  dt[, dist_start := tx_5p_pos - utr5_len - 1L]
  dt[, dist_stop := tx_5p_pos + L - 1L - utr5_len - cds_len]

  p1 <- ggplot(
    dt[L >= 20L & L <= 35L, .(count = sum(count)), by = .(L, pos_frame)],
    aes(L, count, fill = factor(pos_frame))
  ) + geom_col() + theme_classic() +
    labs(title = "Fragment lengths by frame", fill = "5' frame")

  p2 <- plot_end_profile(
    dt, "dist_start", "5' end relative to start codon"
  )

  p3 <- plot_end_profile(
    dt, "dist_stop", "3' end relative to stop codon"
  )

  geometry <- dt[, .(total = sum(count)), by = .(d5, d3)]
  p4 <- ggplot(geometry, aes(factor(d5), factor(d3), fill = total)) +
    geom_tile(color = "black") +
    geom_text(aes(label = paste0(round(total / sum(total) * 100, 1), "%"))) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    theme_classic() +
    labs(title = "d5 versus d3 geometry", x = "d5", y = "d3")

  (p1 + p4) / (p2 + p3)
}

plot_bias_eval <- function(bias_before, bias_after,
                           title = "CHOROS bias evaluation") {
  bias <- data.table(
    position = -6:6,
    name = c(
      "n6", "n5", "n4", "n3", "E", "P", "A_cod",
      "p1", "p2", "p3", "p4", "p5", "p6"
    ),
    raw = as.numeric(bias_before),
    corrected = as.numeric(bias_after)
  )
  bias[, type := "other"]
  bias[name %in% c("n5", "n4", "p3", "p4"), type := "bias"]
  bias[name == "E", type := "E"]
  bias[name == "P", type := "P"]
  bias[name == "A_cod", type := "A"]
  colors <- c(
    E = "#E41A1C", P = "#377EB8", A = "#4DAF4A",
    bias = "#984EA3", other = "grey60"
  )
  maximum <- max(c(bias$raw, bias$corrected, 1e-6), na.rm = TRUE)
  make_plot <- function(value, plot_title) {
    ggplot(bias, aes(position, .data[[value]], fill = type)) +
      geom_col() + theme_bw() + scale_fill_manual(values = colors) +
      scale_x_continuous(breaks = -6:6) +
      coord_cartesian(ylim = c(0, maximum)) +
      labs(title = plot_title, y = expression(Sigma(beta^2))) +
      theme(legend.position = "none")
  }
  make_plot("raw", "Bias before correction") +
    make_plot("corrected", "Bias after correction") +
    patchwork::plot_annotation(title = title)
}

apply_choros_bias_correction <- function(dt, tx_seqs, fit) {
  result <- copy(dt)
  if (!"A_start" %in% names(result)) {
    result[, A_start := utr5_len + 3L * (cod_idx - 1L) + 1L]
  }
  tid <- as.character(result$transcript_id)
  result[, f5 := stri_sub(tx_seqs[tid], A_start - d5, A_start - d5 + 2L)]
  result[, f3 := stri_sub(tx_seqs[tid], A_start + d3, A_start + 2L + d3)]
  compute_gc_dt(result, tx_seqs)

  coefficients <- coef(fit)
  coefficients[is.na(coefficients)] <- 0
  names(coefficients) <- gsub("::", "", names(coefficients))
  get_values <- function(coefficient_names) {
    values <- numeric(length(coefficient_names))
    matched <- coefficient_names %in% names(coefficients)
    values[matched] <- coefficients[coefficient_names[matched]]
    values
  }
  correction_values <- function(sequence, digest, seq_prefix, digest_prefix) {
    marginal <- get_values(paste0(seq_prefix, sequence))
    order1 <- get_values(paste0(
      digest_prefix, digest, ":", seq_prefix, sequence
    ))
    order2 <- get_values(paste0(
      seq_prefix, sequence, ":", digest_prefix, digest
    ))
    exp(marginal + ifelse(order1 != 0, order1, order2))
  }

  result[, cf5 := correction_values(
    as.character(f5), as.character(d5), "f5", "d5"
  )]
  result[, cf3 := correction_values(
    as.character(f3), as.character(d3), "f3", "d3"
  )]
  result[!is.finite(cf5), cf5 := 1]
  result[!is.finite(cf3), cf3 := 1]
  result[, cf5 := pmin(pmax(cf5, 0.2), 5)]
  result[, cf3 := pmin(pmax(cf3, 0.2), 5)]

  b_gc <- if ("gc" %in% names(coefficients)) coefficients[["gc"]] else 0
  mean_gc <- mean(result$gc, na.rm = TRUE)
  result[, c_gc := exp(b_gc * gc - b_gc * mean_gc)]
  result[!is.finite(c_gc), c_gc := 1]
  result[, corrected := count / (cf5 * cf3 * c_gc)]
  if (sum(result$corrected, na.rm = TRUE) > 0) {
    result[, corrected := corrected *
             sum(count, na.rm = TRUE) / sum(corrected, na.rm = TRUE)]
  }
  result
}

plot_metagenes <- function(dt_corrected) {
  dt <- copy(dt_corrected[, .(
    transcript_id, cod_idx, cds_len, count, corrected
  )])
  dt[, dist_start := cod_idx - 1L]
  dt[, dist_stop := cod_idx - cds_len %/% 3L]
  tis <- dt[dist_start >= -10L & dist_start <= 40L, .(
    Raw = sum(count, na.rm = TRUE),
    Corrected = sum(corrected, na.rm = TRUE)
  ), by = .(Position = dist_start)]
  tts <- dt[dist_stop >= -40L & dist_stop <= 10L, .(
    Raw = sum(count, na.rm = TRUE),
    Corrected = sum(corrected, na.rm = TRUE)
  ), by = .(Position = dist_stop)]
  tis <- melt(tis, id.vars = "Position", variable.name = "Type")
  tts <- melt(tts, id.vars = "Position", variable.name = "Type")
  colors <- c(Raw = "grey60", Corrected = "#E41A1C")
  p_tis <- ggplot(tis, aes(Position, value, color = Type)) +
    geom_line(linewidth = 1) + theme_classic() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = colors) +
    labs(title = "TIS metagene", y = "Total count")
  p_tts <- ggplot(tts, aes(Position, value, color = Type)) +
    geom_line(linewidth = 1) + theme_classic() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = colors) +
    labs(title = "TTS metagene", y = "Total count")
  p_tis + p_tts + patchwork::plot_layout(guides = "collect")
}

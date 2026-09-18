test_that("load_seq_bias returns the default AA p-site bias table", {
  dt <- load_seq_bias()
  expect_s3_class(dt, "data.table")
  expect_equal(unique(dt$variable), "median")
  expect_true(all(c("seqs", "alpha") %in% colnames(dt)))
})

test_that("DMN alpha scaling changes variance but preserves expected coverage", {
  alpha <- replicate(4000L, c(2, 6), simplify = FALSE)
  low <- do.call(rbind, scale_dmn_alpha(alpha, 0.25))
  high <- do.call(rbind, scale_dmn_alpha(alpha, 4))

  expect_equal(low[1, ] / sum(low[1, ]), c(0.25, 0.75))
  expect_equal(high[1, ] / sum(high[1, ]), c(0.25, 0.75))

  set.seed(2026)
  low_counts <- extraDistr::rdirmnom(4000L, 1000L, low)
  set.seed(2026)
  high_counts <- extraDistr::rdirmnom(4000L, 1000L, high)
  low_fraction <- low_counts[, 1] / rowSums(low_counts)
  high_fraction <- high_counts[, 1] / rowSums(high_counts)

  expect_equal(mean(low_fraction), 0.25, tolerance = 0.015)
  expect_equal(mean(high_fraction), 0.25, tolerance = 0.015)
  expect_gt(stats::var(low_fraction), 5 * stats::var(high_fraction))
  expect_equal(rowSums(low_counts), rep(1000, 4000L))
  expect_equal(rowSums(high_counts), rep(1000, 4000L))
  for (invalid in list(0, -1, NA_real_, Inf, c(1, 2), "1")) {
    expect_error(scale_dmn_alpha(alpha[1], invalid), "dmn_alpha_scale")
  }
})

test_that("DMN moment estimator recovers a known concentration", {
  set.seed(20260919)
  sites <- 60L
  reads <- 1200L
  true_scale <- 0.4
  alpha <- rep(3 * true_scale, sites)
  counts <- extraDistr::rdirmnom(400L, reads, alpha)
  estimates <- vapply(seq_len(nrow(counts)), function(i) {
    dmn_alpha_moment(
      counts[i, ], expected = rep(reads / sites, sites),
      nt_positions = 3L * sites
    )$dmn_alpha_scale
  }, numeric(1))

  expect_equal(stats::median(estimates), true_scale, tolerance = 0.08)
  expect_true(all(estimates > 0))
})

test_that("DMN alpha defaults to one and learned profiles override it", {
  profile <- data.table::data.table(
    variable = "learned", seqs = c("AAA", "CCC"), alpha = c(1, 2),
    dmn_alpha_scale = 0.025
  )
  expect_equal(resolve_dmn_alpha_scale(NULL, NULL), 1)
  expect_equal(resolve_dmn_alpha_scale(NULL, profile), 0.025)
  expect_equal(resolve_dmn_alpha_scale(2, profile), 2)
  profile[2, dmn_alpha_scale := 0.05]
  expect_error(resolve_dmn_alpha_scale(NULL, profile), "multiple")
})

test_that("numeric autocorrelation kernels smooth locally and validate input", {
  signal <- c(0, 0, 9, 0, 0)
  smoothed <- apply_autocorrelation_kernel(signal, c(1, 2, 1))
  expect_equal(smoothed, c(0, 2.25, 4.5, 2.25, 0))
  expect_equal(apply_autocorrelation_kernel(signal, 1), signal)
  expect_error(apply_autocorrelation_kernel(signal, c(1, 1)), "odd")
  expect_error(apply_autocorrelation_kernel(signal, c(1, -1, 1)), "non-negative")
})

test_that("load_seq_bias supports alternate built-in bias tables", {
  aa_stop <- load_seq_bias(bias = "stop_codon")
  codon_start <- load_seq_bias(
    type = "codon", shift = "a-site", bias = "start_codon"
  )
  codon_similar <- load_seq_bias(
    type = "codon", shift = "a-site", bias = "similar"
  )
  codon_all <- load_seq_bias(type = "codon", shift = "a-site", bias = "all")

  expect_equal(unique(aa_stop$variable), "R1")
  expect_equal(unique(codon_start$variable), "R2")
  expect_equal(unique(codon_similar$variable), "R10")
  expect_true(length(unique(codon_all$variable)) > 1)
})

test_that("load_seq_bias rejects unsupported arguments", {
  expect_error(load_seq_bias(type = "peptide"))
  expect_error(load_seq_bias(shift = "e-site"))
  expect_error(load_seq_bias(bias = "unknown"))
})

test_that("pack_alpha_rows and flatten_sample_rows preserve row-wise ragged layout", {
  alpha_rows <- list(
    c(1, 2, 3),
    c(4, 5),
    c(6, 7, 8, 9)
  )
  region_length_matrix <- rbind(
    c(TRUE, TRUE, TRUE, FALSE),
    c(TRUE, TRUE, FALSE, FALSE),
    c(TRUE, TRUE, TRUE, TRUE)
  )

  packed <- pack_alpha_rows(alpha_rows, region_length_matrix)
  expected_packed <- t(region_length_matrix)
  expected_packed[expected_packed] <- unlist(alpha_rows, use.names = FALSE)
  expected_packed <- t(expected_packed)
  expected_packed[region_length_matrix == FALSE] <- 1e-24

  expect_equal(packed, expected_packed)

  sampled <- matrix(
    c(
      10, 11, 12, 13,
      20, 21, 22, 23,
      30, 31, 32, 33
    ),
    nrow = 3,
    byrow = TRUE
  )
  expected_flat <- t(sampled)[t(region_length_matrix)]

  expect_equal(flatten_sample_rows(sampled, lengths(alpha_rows)), expected_flat)
})

test_that("replace_letters_by_width_group matches per-sequence replacement on mixed widths", {
  sequences <- Biostrings::DNAStringSet(c(
    "ATGTAAATG",
    "ATGCCCTAA",
    "ATGATG",
    "ATGATGATGATG"
  ))
  positions <- IRanges::IntegerList(
    c(4L, 7L),
    integer(),
    1L,
    c(2L, 6L, 11L)
  )

  expected <- Biostrings::DNAStringSet(lapply(seq_along(sequences), function(i) {
    at <- as.integer(positions[[i]])
    if (!length(at)) {
      return(sequences[[i]])
    }
    Biostrings::replaceLetterAt(
      sequences[[i]],
      at = at,
      letter = Biostrings::DNAString(paste(rep.int("C", length(at)), collapse = ""))
    )
  }))

  actual <- replace_letters_by_width_group(sequences, positions)

  expect_equal(as.character(actual), as.character(expected))
})

test_that("uorf debug warnings explain how to avoid interactive debug mode", {
  expect_warning(
    uorf_debug_warning("Some uORFs have internal stop codons that could not be fixed!"),
    regexp = "set debug_on = FALSE"
  )
  expect_warning(
    uorf_debug_warning("Some uORFs have no stop codon."),
    regexp = "longer leaders, fewer uORFs, or disabled overlaps"
  )
})

test_that("distribute_reads_to_uORFs preserves uorf_ranges order when expanding transcript counts", {
  region_counts <- c(txB = 7, txA = 5)
  assay <- matrix(region_counts, ncol = 1, dimnames = list(names(region_counts), "sample"))
  uorf_ranges <- GenomicRanges::GRanges(
    seqnames = c("chr2", "chr1", "chr2"),
    ranges = IRanges::IRanges(start = c(10, 20, 30), width = 3),
    strand = c("+", "+", "+")
  )
  names(uorf_ranges) <- c("txB", "txA", "txB")

  set.seed(1)
  allocated <- distribute_reads_to_uORFs(
    region_counts = region_counts,
    assay = assay,
    uorf_ranges = uorf_ranges,
    uorf_prop_mode = "character",
    uorf_prop_within_gene = "uniform"
  )

  expect_equal(names(allocated), ORFik::txNames(uorf_ranges))
  expect_equal(sum(allocated[names(allocated) == "txA"]), unname(region_counts["txA"]))
  expect_equal(sum(allocated[names(allocated) == "txB"]), unname(region_counts["txB"]))
})

test_that("libFormats validation accepts supported formats and rejects others", {
  expect_no_error(validate_lib_formats(list(RFP = "ofst", RNA = "bam", CAGE = "sam")))
  expect_equal(resolve_export_format(list(RFP = "bam"), "RFP"), "bam")
  expect_error(validate_lib_formats(list(RFP = "cram")))
  expect_error(resolve_export_format(list(RFP = c("ofst", "bam")), "RFP"))
})

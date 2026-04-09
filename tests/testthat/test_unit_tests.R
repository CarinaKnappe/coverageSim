test_that("load_seq_bias returns the default AA p-site bias table", {
  dt <- load_seq_bias()
  expect_s3_class(dt, "data.table")
  expect_equal(unique(dt$variable), "R2")
  expect_true(all(c("seqs", "alpha") %in% colnames(dt)))
})

test_that("load_seq_bias supports alternate built-in bias tables", {
  aa_stop <- load_seq_bias(bias = "stop_codon")
  codon_all <- load_seq_bias(type = "codon", shift = "a-site", bias = "all")

  expect_equal(unique(aa_stop$variable), "R1")
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

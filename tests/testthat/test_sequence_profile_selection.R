test_that("named built-in profiles retain their original weights", {
  aliases <- c(start_codon = "R2", stop_codon = "R1", similar = "R10")
  for (type in c("AA", "codon")) {
    for (shift in c("p-site", "a-site")) {
      all_profiles <- load_seq_bias(type = type, shift = shift, bias = "all")
      for (profile in unique(all_profiles$variable)) {
        selected <- load_seq_bias(type = type, shift = shift, bias = profile)
        expect_equal(selected, all_profiles[variable == profile])
      }
      for (alias in names(aliases)) {
        expect_equal(
          load_seq_bias(type = type, shift = shift, bias = alias),
          load_seq_bias(type = type, shift = shift, bias = aliases[[alias]])
        )
      }
    }
  }
})

test_that("ambiguous profiles fail before simulation reads or writes files", {
  profiles <- load_seq_bias(type = "codon", shift = "a-site", bias = "all")
  expect_error(simNGScoverage(seq_bias = profiles), "exactly one named sequence profile")
  expect_error(
    simNGScoverage(seq_bias = profiles[nrow(profiles):1L]),
    "exactly one named sequence profile"
  )
  for (invalid in list(character(), NA_character_, "", c("R1", "R2"))) {
    expect_error(load_seq_bias(bias = invalid), "bias must be one")
  }
  expect_error(load_seq_bias(bias = "R999"), "Unknown bias profile")
})

test_that("sequence weights come from the chosen profile without changing input", {
  fasta <- tempfile(fileext = ".fa")
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(c(chr1 = "ATGAAATAA")), fasta
  )
  Rsamtools::indexFa(fasta)
  cds <- GenomicRanges::GRangesList(
    tx1 = GenomicRanges::GRanges("chr1", IRanges::IRanges(1L, 9L), "+")
  )
  profiles <- data.table::data.table(
    variable = rep(c("first", "chosen"), each = 3L),
    seqs = rep(c("ATG", "AAA", "TAA"), 2L),
    alpha = c(1, 2, 3, 10, 20, 30)
  )
  weights <- function(table) {
    add_sequence_bias(
      c(genome = fasta), data.table::data.table(genes = rep(1L, 9L)),
      table, cds, 9L, "cds"
    )
  }
  expect_error(weights(profiles), "exactly one named sequence profile")
  chosen <- profiles[variable == "chosen"]
  original <- data.table::copy(chosen)
  expect_equal(unname(unlist(weights(chosen))), c(10, 20, 30))
  expect_equal(chosen, original)
  expect_equal(unname(unlist(weights(profiles[variable == "first"]))), c(1, 2, 3))
  expect_equal(unname(unlist(weights(chosen[, .(seqs, alpha)]))), c(10, 20, 30))
  for (invalid in list(NA_character_, "")) {
    malformed <- data.table::copy(chosen)
    malformed[, variable := invalid]
    expect_error(weights(malformed), "exactly one named sequence profile")
  }
})

make_artifact_records <- function() {
  fixture <- end_selection_fixture()
  fixture$signal[, score := c(10L, 10L)]
  fragments <- make_simulated_rpf_fragments(
    fixture$signal, fixture$models, 8L,
    list(source = "user", site_offset = 2L)
  )
  set.seed(901)
  simulate_alignment_artifacts(fragments, list(
    duplication_rate = 1, duplicate_copies = 2L,
    multimapping_rate = 1, secondary_alignments = 1L
  ))
}

test_that("PCR and multimapping artifacts create explicit SAM records", {
  records <- make_artifact_records()
  diagnostics <- attr(records, "diagnostics")
  expect_equal(diagnostics[["biological_molecules"]], 20)
  expect_equal(diagnostics[["pcr_duplicate_records"]], 40)
  expect_equal(diagnostics[["multimapping_molecules"]], 20)
  expect_equal(diagnostics[["secondary_alignment_records"]], 20)
  expect_equal(nrow(records), 80)
  expect_equal(sum(bitwAnd(records$flag, 1024L) != 0L), 40)
  expect_equal(sum(bitwAnd(records$flag, 256L) != 0L), 20)
  expect_equal(sum(records$is_duplicate), 40)
  expect_equal(sum(records$is_secondary), 20)
  expect_true(all(records[is_secondary == TRUE, nh] == 2L))
  expect_true(all(records[is_secondary == TRUE, qname] %in%
                  records[is_secondary == FALSE & is_duplicate == FALSE, qname]))
})

test_that("every record of a multimapper carries that molecule's own sequence", {
  # A multimapper is one molecule the aligner could place in more than one spot.
  # The coordinates differ between its records, but the sequence that was read
  # does not, so a tool that compares the records of one QNAME must see one read.
  # The secondary record borrows its coordinates from another fragment, and it
  # must not borrow that fragment's sequence along with them.
  records <- make_artifact_records()
  reverse_complement <- function(x) {
    as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))
  }
  groups <- records[, .(
    reads = uniqueN(sequence),
    # SEQ is written relative to the genome, so it may legitimately differ
    # between records that sit on opposite strands -- but only in that way.
    oriented_correctly = all(
      reference_sequence == ifelse(strand == "+", sequence,
                                   reverse_complement(sequence))
    )
  ), by = qname]
  multi <- groups[records[, .N, by = qname][N > 1L], on = "qname"]
  expect_gt(nrow(multi), 0L)
  expect_true(all(multi$reads == 1L))
  expect_true(all(groups$oriented_correctly))
  # And the sequence length still matches the CIGAR, which is what makes the
  # record readable at all.
  expect_true(all(nchar(records$reference_sequence) == records$fragment_length))
})

test_that("artifact BAM preserves duplicate flags, secondary flags, and NH tags", {
  records <- make_artifact_records()
  model <- end_selection_fixture()$models$tx
  seqinfo <- GenomeInfoDb::seqinfo(model$exons)
  GenomeInfoDb::seqlengths(seqinfo) <- 90L
  base <- tempfile("artifact-bam-")
  paths <- write_artifact_library(records, base, "bam", seqinfo)
  raw <- Rsamtools::scanBam(
    paths[["default"]],
    param = Rsamtools::ScanBamParam(what = c("qname", "flag"), tag = "NH")
  )[[1]]
  expect_length(raw$flag, 80)
  expect_equal(sum(bitwAnd(raw$flag, 1024L) != 0L), 40)
  expect_equal(sum(bitwAnd(raw$flag, 256L) != 0L), 20)
  expect_equal(sum(raw$tag$NH > 1L), 40)
})

test_that("technical artifact defaults are clean and settings are validated", {
  defaults <- normalize_technical_artifacts(NULL)
  expect_false(has_active_technical_artifacts(defaults))
  expect_equal(defaults$duplication_rate, 0)
  expect_equal(defaults$multimapping_rate, 0)
  expect_error(normalize_technical_artifacts(list(duplication_rate = 1.1)),
               "duplication_rate")
  expect_error(normalize_technical_artifacts(list(secondary_alignments = 0)),
               "secondary_alignments")
  expect_error(normalize_technical_artifacts(list(unknown = 1)), "Unknown")
})

test_that("a secondary alignment on the opposite strand stores the reverse complement", {
  # The case that matters for orientation is a molecule whose alternative
  # location sits on the other strand. SAM reports SEQ relative to the genome,
  # so the two records of that one molecule must hold reverse complements of
  # each other -- the molecule was read once, and only the strand it was placed
  # on differs. A plus-strand-only fixture cannot show this, because there the
  # two sequence columns are identical and leaving the orientation out would
  # look correct.
  forward <- "AAAACCCCGGGGTTTTAAAACCCCGGGGTT"
  reverse <- as.character(Biostrings::reverseComplement(
    Biostrings::DNAString(forward)
  ))
  other <- "GGGGTTTTAAAACCCCGGGGTTTTAAAACC"
  fragments <- data.table::data.table(
    seqnames = "chr1",
    start = c(100L, 500L),
    strand = c("+", "-"),
    cigar = "30M",
    fragment_length = 30L,
    score = 1L,
    sequence = c(forward, other),
    reference_sequence = c(
      forward,
      as.character(Biostrings::reverseComplement(Biostrings::DNAString(other)))
    )
  )
  set.seed(3)
  records <- simulate_alignment_artifacts(fragments, list(
    multimapping_rate = 1, secondary_alignments = 1L
  ))
  # Each molecule is now reported twice, once per strand.
  expect_equal(nrow(records), 4L)
  expect_equal(records[, .N, by = qname][, unique(N)], 2L)

  plus_molecule <- records[sequence == forward]
  expect_equal(nrow(plus_molecule), 2L)
  expect_setequal(plus_molecule$strand, c("+", "-"))
  # Read once, so one sequence in read direction ...
  expect_equal(plus_molecule[, uniqueN(sequence)], 1L)
  # ... but two different ones as written to SAM, and they are reverse
  # complements rather than two unrelated strings.
  expect_equal(plus_molecule[strand == "+", reference_sequence], forward)
  expect_equal(plus_molecule[strand == "-", reference_sequence], reverse)
  expect_true(all(nchar(records$reference_sequence) == 30L))
})

make_fragment_fixture <- function() {
  chromosome <- paste(rep(c("A", "C", "G", "T"), 125L), collapse = "")
  fasta <- tempfile(fileext = ".fa")
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(c(chr1 = chromosome)),
    filepath = fasta
  )
  Rsamtools::indexFa(fasta)

  transcripts <- GenomicRanges::GRangesList(
    tx_plus = GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(c(101L, 201L), width = 10L), "+",
      exon_rank = 1:2
    ),
    tx_minus = GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(c(401L, 301L), width = 10L), "-",
      exon_rank = 1:2
    )
  )
  GenomeInfoDb::seqlengths(transcripts) <- c(chr1 = 500L)
  list(
    fasta = fasta,
    transcripts = transcripts,
    models = transcript_models(transcripts, fasta)
  )
}

make_fragment_signals <- function(scores = c(2L, 3L)) {
  data.table::data.table(
    seqnames = "chr1",
    start = c(110L, 401L),
    end = c(110L, 401L),
    strand = c("+", "-"),
    score = scores,
    transcript_id = c("tx_plus", "tx_minus"),
    signal_position = c(110L, 401L)
  )
}

test_that("simulated RPFs preserve site geometry, sequence, and strand", {
  fixture <- make_fragment_fixture()
  fragments <- make_simulated_rpf_fragments(
    make_fragment_signals(), fixture$models, 8L,
    list(source = "user", site_offset = 2L, boundary_action = "error")
  )

  expect_equal(fragments$cigar, c("3M90N5M", "5M90N3M"))
  expect_equal(fragments$fragment_length, c(8L, 8L))
  expect_equal(nchar(fragments$sequence), c(8L, 8L))
  expect_equal(fragments$site_offset, c(2L, 2L))
  expect_equal(fragments$ribosome_site, fragments$signal_position)
  expect_false(any(fragments$ribosome_site == fragments$fragment_start))
  expect_equal(
    fragments$ribosome_site[fragments$strand == "+"] -
      fragments$five_prime_end[fragments$strand == "+"],
    2L
  )
  expect_equal(
    fragments$five_prime_end[fragments$strand == "-"] -
      fragments$ribosome_site[fragments$strand == "-"],
    2L
  )

  expected <- vapply(seq_len(nrow(fragments)), function(i) {
    model <- fixture$models[[fragments$transcript_id[i]]]
    site_tx <- genomic_site_to_transcript(model, fragments$signal_position[i])
    substr(model$sequence, site_tx - 2L, site_tx + 5L)
  }, character(1))
  expect_equal(fragments$sequence, unname(expected))
})

test_that("OFST and BAM represent the same simulated RPFs", {
  fixture <- make_fragment_fixture()
  fragments <- make_simulated_rpf_fragments(
    make_fragment_signals(), fixture$models, 8L,
    list(source = "user", site_offset = 2L, boundary_action = "error")
  )
  alignments <- simulated_rpf_alignments(
    fragments, GenomeInfoDb::seqinfo(fixture$transcripts)
  )
  expect_equal(
    suppressWarnings(as.integer(ORFik::readWidths(alignments))),
    c(8L, 8L)
  )
  expect_equal(as.character(GenomicAlignments::cigar(alignments)), fragments$cigar)

  ofst <- tempfile(fileext = ".ofst")
  bam <- tempfile(fileext = ".bam")
  write_ofst_library(alignments, ofst)
  write_bam_library(alignments, bam)
  imported <- ORFik::fimport(ofst)
  scanned <- Rsamtools::scanBam(
    bam,
    param = Rsamtools::ScanBamParam(what = c("pos", "cigar", "seq", "strand"))
  )[[1]]

  expect_equal(as.integer(GenomicRanges::start(imported)), fragments$fragment_start)
  expect_equal(as.character(GenomicAlignments::cigar(imported)), fragments$cigar)
  expect_equal(
    suppressWarnings(as.integer(ORFik::readWidths(imported))),
    fragments$fragment_length
  )
  expect_equal(length(scanned$pos), sum(fragments$score))
  expect_equal(unique(scanned$cigar), fragments$cigar)
  expect_true(all(Biostrings::width(scanned$seq) == 8L))
  expect_setequal(as.character(scanned$seq), fragments$sequence)
})

test_that("legacy point export keeps the signal as alignment start", {
  point <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(100L, width = 28L), "+", score = 1L,
    seqinfo = GenomeInfoDb::Seqinfo("chr1", seqlengths = 500L)
  )
  sam <- tempfile(fileext = ".sam")
  samFromGAlignment(point, sam, make_bam = FALSE)
  fields <- strsplit(readLines(sam)[2], "\t", fixed = TRUE)[[1]]
  expect_equal(as.integer(fields[4]), 100L)
  expect_equal(fields[6], "28M")
  expect_equal(fields[10], "*")
})

test_that("boundary failures report affected counts instead of dropping them", {
  fixture <- make_fragment_fixture()
  edge <- make_fragment_signals(scores = c(7L, 0L))[1]
  edge[, `:=`(start = 101L, end = 101L, signal_position = 101L)]
  expect_error(
    make_simulated_rpf_fragments(
      edge, fixture$models, 8L,
      list(source = "user", site_offset = 2L, boundary_action = "error")
    ),
    "affected count: 7"
  )
})

test_that("ground truth contains neutral simulated-RPF fields", {
  fixture <- make_fragment_fixture()
  fragments <- make_simulated_rpf_fragments(
    make_fragment_signals(), fixture$models, 8L,
    list(source = "user", site_offset = 2L, boundary_action = "error")
  )
  base <- tempfile()
  path <- write_fragment_ground_truth(fragments, base, TRUE)
  truth <- data.table::fread(path)
  expect_equal(sum(truth$score), sum(fragments$score))
  expect_true(all(c(
    "fragment_id", "transcript_id", "ribosome_site", "signal_position",
    "fragment_start", "fragment_end", "fragment_length", "strand",
    "site_reference", "geometry_probability", "sequence", "score"
  ) %in% names(truth)))
})

test_that("joint geometry probabilities split counts without breaking coupling", {
  fixture <- make_fragment_fixture()
  signal <- make_fragment_signals(scores = c(10000L, 0L))[1]
  distribution <- data.frame(
    fragment_length = c(8L, 9L),
    site_offset = c(2L, 4L),
    probability = c(0.8, 0.2)
  )

  set.seed(42)
  fragments <- make_simulated_rpf_fragments(
    signal, fixture$models, fragment_lengths = 8:9,
    fragment_geometry = list(
      source = "user", site_reference = "p_site",
      distribution = distribution, boundary_action = "error"
    )
  )

  expect_equal(sum(fragments$score), 10000L)
  expect_setequal(
    paste(fragments$fragment_length, fragments$site_offset),
    c("8 2", "9 4")
  )
  expect_equal(fragments$score[fragments$fragment_length == 8L] / 10000,
               0.8, tolerance = 0.02)
  expect_equal(fragments$geometry_probability, c(0.8, 0.2))
  expect_true(all(fragments$site_reference == "p_site"))
})

test_that("default A-site offsets are three nucleotides beyond P-site offsets", {
  p_site <- default_fragment_distribution(28:32, "p_site")
  a_site <- default_fragment_distribution(28:32, "a_site")

  expect_equal(a_site$fragment_length, p_site$fragment_length)
  expect_equal(a_site$site_offset, p_site$site_offset + 3L)
  expect_equal(a_site$probability, p_site$probability)
})

test_that("simulated RPF geometry defaults to A-site", {
  geometry <- normalize_fragment_geometry(NULL)
  default_distribution <- default_fragment_distribution(28:30)
  explicit_a_site <- default_fragment_distribution(28:30, "a_site")

  expect_identical(geometry$site_reference, "a_site")
  expect_equal(default_distribution, explicit_a_site)
})

test_that("boundary renormalization preserves all counts", {
  fixture <- make_fragment_fixture()
  edge <- make_fragment_signals(scores = c(101L, 0L))[1]
  edge[, `:=`(start = 103L, end = 103L, signal_position = 103L)]
  distribution <- data.frame(
    fragment_length = c(8L, 8L),
    site_offset = c(2L, 6L),
    probability = c(0.25, 0.75)
  )
  fragments <- make_simulated_rpf_fragments(
    edge, fixture$models, 8L,
    list(source = "learned", distribution = distribution,
         boundary_action = "renormalize")
  )

  expect_equal(sum(fragments$score), 101L)
  expect_equal(fragments$site_offset, 2L)
  expect_equal(fragments$geometry_probability, 1)
})

test_that("former physical function names remain compatibility aliases", {
  fixture <- make_fragment_fixture()
  expect_warning(
    fragments <- make_physical_fragments(
      make_fragment_signals(), fixture$models, 8L,
      list(source = "user", site_offset = 2L, boundary_action = "error")
    ),
    "deprecated"
  )
  expect_warning(
    alignments <- physical_fragment_alignments(
      fragments, GenomeInfoDb::seqinfo(fixture$transcripts)
    ),
    "deprecated"
  )
  expect_s4_class(alignments, "GAlignments")
})

test_that("end-bias API is neutral until an explicit model is implemented", {
  expect_no_error(normalize_fragment_geometry(list(
    source = "default",
    five_prime_bias = list(source = "none"),
    three_prime_bias = list(source = "none")
  )))
  expect_error(
    normalize_fragment_geometry(list(
      source = "default",
      five_prime_bias = list(source = "user")
    )),
    "reserved"
  )
})

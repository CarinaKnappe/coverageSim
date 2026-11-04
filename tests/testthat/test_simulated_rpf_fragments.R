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

  # reference_sequence is what must be written to SAM/BAM SEQ: identical to
  # the transcript-oriented `sequence` on the plus strand, and its reverse
  # complement on the minus strand (SAM always reports SEQ relative to the
  # reference/plus strand, regardless of the alignment's FLAG).
  expect_equal(
    fragments$reference_sequence[fragments$strand == "+"],
    fragments$sequence[fragments$strand == "+"]
  )
  expect_equal(
    fragments$reference_sequence[fragments$strand == "-"],
    as.character(Biostrings::reverseComplement(
      Biostrings::DNAStringSet(fragments$sequence[fragments$strand == "-"])
    ))
  )
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
  # The BAM SEQ field is always relative to the reference (plus) strand, so
  # it must match `reference_sequence`, not the transcript-oriented
  # `sequence` (they differ for the minus-strand fragment in this fixture).
  expect_setequal(as.character(scanned$seq), fragments$reference_sequence)
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

test_that("a library without reads gives an empty fragment table with every column", {
  fixture <- make_fragment_fixture()
  # Positions with signal, but a read budget of zero.
  no_reads <- make_fragment_signals()[, read_count := 0L]
  fragments <- make_simulated_rpf_fragments(
    no_reads, fixture$models, 8L,
    list(source = "user", site_offset = 2L, boundary_action = "error")
  )
  with_reads <- make_simulated_rpf_fragments(
    make_fragment_signals(), fixture$models, 8L,
    list(source = "user", site_offset = 2L, boundary_action = "error")
  )
  expect_equal(nrow(fragments), 0L)
  expect_identical(names(fragments), names(with_reads))
  expect_length(simulated_rpf_alignments(
    fragments, GenomeInfoDb::seqinfo(fixture$transcripts)
  ), 0L)
})

test_that("a signal table without positions gives the same empty fragment table", {
  fixture <- make_fragment_fixture()
  geometry <- list(source = "user", site_offset = 2L, boundary_action = "error")
  no_positions <- make_simulated_rpf_fragments(
    make_fragment_signals()[0], fixture$models, 8L, geometry
  )
  no_reads <- make_simulated_rpf_fragments(
    make_fragment_signals()[, read_count := 0L], fixture$models, 8L, geometry
  )
  expect_equal(nrow(no_positions), 0L)
  expect_identical(lapply(no_positions, class), lapply(no_reads, class))
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

test_that("positions with no room for a whole fragment pass their reads on", {
  fixture <- make_fragment_fixture()
  # tx_plus is 20 nt: genomic 101-110 are transcript positions 1-10 and
  # 201-210 are 11-20. With 14-nt fragments and offsets 9-12, nothing fits at
  # position 1 or 20, while at position 12 the offsets 9-11 fit.
  signals <- make_fragment_signals()[c(1, 1, 1)]
  signals[, `:=`(signal_position = c(101L, 202L, 210L), score = 60L)]
  signals[, `:=`(start = signal_position, end = signal_position)]
  geometry <- list(source = "user", site_offset = 9:12, boundary_action = "renormalize")
  fragments <- make_simulated_rpf_fragments(signals, fixture$models, 14L, geometry)

  # All reads of the transcript end up at the one position with room, with
  # unchanged offsets, and nothing is left out.
  expect_identical(unique(fragments$signal_position), 202L)
  expect_true(all(fragments$site_offset %in% 9:11))
  expect_equal(sum(fragments$score), 180L)
  expect_equal(nrow(attr(fragments, "dropped_reads")), 0L)
  expect_error(
    make_simulated_rpf_fragments(
      signals, fixture$models, 14L, modifyList(geometry, list(boundary_action = "error"))
    ),
    "Cannot create a complete fragment"
  )
})

test_that("reads stay in their region and are left out where no position has room", {
  fixture <- make_fragment_fixture()
  signals <- make_fragment_signals()[c(1, 1, 1)]
  signals[, `:=`(signal_position = c(101L, 202L, 210L), score = c(60L, 30L, 20L),
                 sampling_group = c("leader:1", "cds:1", "cds:1"))]
  signals[, `:=`(start = signal_position, end = signal_position)]
  plain <- list(source = "user", site_offset = 9:12, boundary_action = "renormalize")
  with_end_bias <- c(plain, list(five_prime_bias = list(
    source = "user", table = data.table::data.table(kmer = c("A", "C", "G", "T"),
                                                    weight = c(2, 1, 1, 1))
  )))

  for (geometry in list(plain, with_end_bias)) {
    fragments <- make_simulated_rpf_fragments(signals, fixture$models, 14L, geometry)
    dropped <- attr(fragments, "dropped_reads")
    # The cds keeps its 50 reads at its one position with room; the leader
    # has no such position, so its 60 reads are left out and reported.
    expect_identical(unique(fragments$signal_position), 202L)
    expect_equal(sum(fragments$score), 50L)
    expect_equal(dropped$sampling_group, "leader:1")
    expect_equal(dropped$reads, 60)
  }

  base <- tempfile()
  path <- write_dropped_reads(dropped, base)
  written <- data.table::fread(path)
  expect_equal(written$region, "leader")
  expect_equal(written$reads, 60L)
  expect_equal(written$transcript_id, "tx_plus")
  # A later run into the same place without losses removes the old report.
  expect_null(write_dropped_reads(dropped[0], base))
  expect_false(file.exists(path))

  expect_error(
    make_simulated_rpf_fragments(signals[, sampling_group := NA_character_],
                                 fixture$models, 14L, plain),
    "must not be missing"
  )
})

test_that("a per-position read_count is what moves when there is no room", {
  fixture <- make_fragment_fixture()
  signals <- make_fragment_signals()[c(1, 1)]
  signals[, `:=`(signal_position = c(101L, 202L), score = 1L, read_count = c(7L, 5L))]
  signals[, `:=`(start = signal_position, end = signal_position)]
  fragments <- make_simulated_rpf_fragments(
    signals, fixture$models, 14L,
    list(source = "user", site_offset = 9:12, boundary_action = "renormalize")
  )
  expect_identical(unique(fragments$signal_position), 202L)
  expect_equal(sum(fragments$score), 12L)
})

test_that("a library whose reads were all left out is written as an empty file", {
  fixture <- make_fragment_fixture()
  signals <- make_fragment_signals()[1][, `:=`(signal_position = 101L, start = 101L, end = 101L)]
  fragments <- make_simulated_rpf_fragments(
    signals, fixture$models, 14L,
    list(source = "user", site_offset = 12L, boundary_action = "renormalize")
  )
  expect_equal(nrow(fragments), 0L)
  expect_equal(attr(fragments, "dropped_reads")$reads, 2)
  alignments <- simulated_rpf_alignments(fragments, GenomeInfoDb::seqinfo(fixture$transcripts))
  files <- write_bam_library(alignments, tempfile(fileext = ".bam"))
  expect_length(GenomicAlignments::readGAlignments(files[["default"]]), 0L)
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

test_that("end-bias API accepts neutral defaults and validates explicit tables", {
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
    "requires a kmer/weight table"
  )
})

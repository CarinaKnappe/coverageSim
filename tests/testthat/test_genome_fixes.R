genome_orf_problems <- function(genome) {
  cds <- ORFik::loadRegion(genome["txdb"], "cds")
  protein <- Biostrings::translate(ORFik::txSeqsFromFa(cds, genome["genome"]))
  stops <- Biostrings::alphabetFrequency(protein)[, "*"]
  last <- as.character(Biostrings::subseq(protein, Biostrings::width(protein),
                                          Biostrings::width(protein)))
  first <- as.character(Biostrings::subseq(protein, 1L, 1L))
  problems <- c(cds_internal_stop = sum(stops > 1L),
                cds_bad_ends = sum(last != "*" | first != "M"))
  if ("uorfs" %in% names(genome)) {
    uorf <- readRDS(genome["uorfs"])
    uorf_protein <- Biostrings::translate(ORFik::txSeqsFromFa(uorf, genome["genome"]))
    problems[["uorf_internal_stop"]] <-
      sum(Biostrings::alphabetFrequency(uorf_protein)[, "*"] > 1L)
  }
  problems
}

test_that("uORFs overlapping the CDS never leave a stop codon inside the CDS", {
  for (settings in list(list(overlap_cds = 2L, max_uorfs = 1),
                        list(overlap_cds = 1L, max_uorfs = 3))) {
    for (seed in 1:2) {
      set.seed(seed)
      genome <- suppressWarnings(suppressMessages(simGenome(
        n = 30, out_dir = tempfile("genome-"), max_uorfs = settings$max_uorfs,
        uorfs_can_overlap = TRUE, uorfs_can_overlap_cds = settings$overlap_cds,
        chromosome_layout = "compact", debug_on = FALSE
      )))
      expect_equal(unname(genome_orf_problems(genome)), c(0, 0, 0))
    }
  }
})

test_that("simGenome supports one and several CDS exons", {
  for (exons in 1:3) {
    set.seed(exons)
    genome <- suppressWarnings(suppressMessages(simGenome(
      n = 6, out_dir = tempfile("genome-"), cds_exons = exons,
      cds_length = rep(600, 6), max_uorfs = 0, debug_on = FALSE
    )))
    expect_equal(unname(genome_orf_problems(genome)), c(0, 0))
    cds <- ORFik::loadRegion(genome["txdb"], "cds")
    expect_true(all(lengths(cds) == exons))
    expect_true(all(ORFik::widthPerGroup(cds, FALSE) == 606))
  }
})

test_that("simGenome() accepts a stop_codons subset, simulating stop-codon usage bias", {
  # stop_codons used to be forced to exactly all three standard stop codons
  # (stopifnot(identical(sort(stop_codons), sort(c("TAG","TGA","TAA"))))),
  # making the parameter pointless -- any non-empty subset (still real
  # stop codons, since the integrity checks below rely on
  # Biostrings::translate()'s default genetic code to recognize them) is
  # now accepted.
  for (subset in list("TGA", c("TAA", "TAG"))) {
    set.seed(5)
    genome <- suppressWarnings(suppressMessages(simGenome(
      n = 20, out_dir = tempfile("genome-"), cds_length = rep(300, 20),
      max_uorfs = 0, debug_on = FALSE, stop_codons = subset
    )))
    expect_equal(unname(genome_orf_problems(genome)), c(0, 0))
    cds <- ORFik::loadRegion(genome["txdb"], "cds")
    cds_seqs <- ORFik::txSeqsFromFa(cds, genome["genome"])
    last_codon <- as.character(Biostrings::subseq(
      cds_seqs, Biostrings::width(cds_seqs) - 2L, Biostrings::width(cds_seqs)
    ))
    expect_true(all(last_codon %in% subset))
  }

  # A non-standard "stop codon", a duplicate, or an empty set are all rejected.
  expect_error(simGenome(n = 2, out_dir = tempfile("genome-"), stop_codons = "AAA"))
  expect_error(simGenome(n = 2, out_dir = tempfile("genome-"), stop_codons = c("TGA", "TGA")))
  expect_error(simGenome(n = 2, out_dir = tempfile("genome-"), stop_codons = character(0)))
})

test_that("simGenome() rejects leader and trailer lengths it cannot build", {
  make <- function(leader = rep(140, 3), trailer = rep(105, 3)) {
    suppressWarnings(suppressMessages(simGenome(
      n = 3, out_dir = tempfile("genome-"), max_uorfs = 0, debug_on = FALSE,
      leader_length = leader, trailer_length = trailer
    )))
  }
  leader_error <- "leader_length must be whole numbers of at least 1"
  trailer_error <- "trailer_length must be whole numbers of at least 0"
  expect_error(make(leader = c(140, 0, 140)), leader_error, fixed = TRUE)
  expect_error(make(leader = rep(-5, 3)), leader_error, fixed = TRUE)
  expect_error(make(leader = c(140, NA, 140)), leader_error, fixed = TRUE)
  expect_error(make(leader = rep(100.5, 3)), leader_error, fixed = TRUE)
  expect_error(make(leader = rep("140", 3)), leader_error, fixed = TRUE)
  expect_error(make(leader = rep(TRUE, 3)), leader_error, fixed = TRUE)
  expect_error(make(trailer = c(105, Inf, 105)), trailer_error, fixed = TRUE)
  expect_error(make(trailer = rep(-1, 3)), trailer_error, fixed = TRUE)
  expect_error(make(trailer = rep(10.5, 3)), trailer_error, fixed = TRUE)
})

test_that("simGenome() builds a one-nucleotide leader and an empty trailer", {
  set.seed(3)
  genome <- suppressWarnings(suppressMessages(simGenome(
    n = 3, out_dir = tempfile("genome-"), max_uorfs = 0, debug_on = FALSE,
    export_txdb = TRUE, leader_length = rep(1, 3), trailer_length = rep(0, 3)
  )))
  lengths <- GenomicFeatures::transcriptLengths(
    ORFik::loadTxdb(genome["txdb"]), with.utr5_len = TRUE, with.utr3_len = TRUE
  )
  expect_equal(lengths$utr5_len, rep(1L, 3))
  expect_equal(lengths$utr3_len, rep(0L, 3))
  expect_equal(unname(genome_orf_problems(genome)), c(0, 0))
})

test_that("simGenome() explains every setting it cannot use", {
  make <- function(...) {
    settings <- utils::modifyList(
      list(n = 3, out_dir = tempfile("genome-"), max_uorfs = 0, debug_on = FALSE),
      list(...)
    )
    suppressWarnings(suppressMessages(do.call(simGenome, settings)))
  }
  codons <- "divisible by 3 (whole codons)"
  cases <- list(
    list(list(n = 1.6), "n must be one whole number of at least 1"),
    list(list(n = 0), "n must be one whole number of at least 1"),
    list(list(cds_length = rep(301, 3)), codons),
    list(list(cds_length = c(300, NA, 300)), "cds_length must be whole numbers of at least 9"),
    list(list(cds_length = c(300, 600)), "cds_length must have one value per gene (n = 3), not 2"),
    list(list(cds_exons = 1.6), "cds_exons must be one whole number of at least 1"),
    list(list(cds_exons = 0), "cds_exons must be one whole number of at least 1"),
    list(list(cds_exons = c(2, 3)), "cds_exons must be one number"),
    list(list(cds_intron_length = rep(31, 3)), codons),
    list(list(cds_intron_length = rep(0, 3)), "cds_intron_length must be whole numbers of at least 3"),
    list(list(max_uorfs = 1.6), "max_uorfs must be one whole number of at least 0"),
    list(list(max_uorfs = -1), "max_uorfs must be one whole number of at least 0"),
    list(list(uorf_max_length = 61), codons),
    list(list(uorf_max_length = 0), "uorf_max_length must be one whole number of at least 6"),
    list(list(uorf_max_length = 600), "every cds_length must be at least uorf_max_length + 3"),
    list(list(uorfs_can_overlap = NA), "uorfs_can_overlap must be TRUE or FALSE"),
    list(list(uorfs_can_overlap_cds = 1.5), "uorfs_can_overlap_cds must be 0 (FALSE), 1 (TRUE) or 2"),
    list(list(flank_length = 1.6), "flank_length must be whole numbers of at least 0"),
    list(list(chromosome_count = 1.6), "chromosome_count must be one whole number of at least 1"),
    list(list(start_codons = "XYZ"), "start_codons must be codons of three letters"),
    list(list(stop_codons = "TTT"), "stop_codons must be one or more of TAG, TGA, TAA"),
    list(list(uorfs_can_overlap_cds = "1"), "uorfs_can_overlap_cds must be 0 (FALSE), 1 (TRUE) or 2"),
    list(list(flank_length = c(100, 200)), "flank_length must have one value per gene (n = 3), not 2"),
    list(list(export_txdb = NA), "export_txdb must be TRUE or FALSE"),
    list(list(debug_on = NA), "debug_on must be TRUE or FALSE")
  )
  for (case in cases) {
    expect_error(do.call(make, case[[1]]), case[[2]], fixed = TRUE,
                 info = paste(names(case[[1]]), format(case[[1]][[1]]), collapse = " "))
  }
})

test_that("simGenome() accepts the per-gene and single-value forms it supports", {
  build <- function(...) {
    set.seed(4)
    suppressWarnings(suppressMessages(simGenome(
      n = 3, out_dir = tempfile("genome-"), max_uorfs = 0, debug_on = FALSE, ...
    )))
  }
  expect_equal(unname(genome_orf_problems(build(cds_intron_length = 3))), c(0, 0))
  expect_equal(unname(genome_orf_problems(build(flank_length = c(100, 200, 300)))), c(0, 0))
  # Without a second exon there is no intron, so its length is not used.
  expect_equal(unname(genome_orf_problems(build(cds_exons = 1,
                                                cds_intron_length = rep(0, 3)))), c(0, 0))
})

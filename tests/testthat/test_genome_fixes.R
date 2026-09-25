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

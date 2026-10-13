seq_bias_matching_fixture <- function() {
  set.seed(11)
  genome <- suppressWarnings(suppressMessages(simGenome(
    n = 4, out_dir = tempfile("genome-"), cds_length = rep(300, 4),
    max_uorfs = 0, debug_on = FALSE
  )))
  cds <- ORFik::loadRegion(genome["txdb"], "cds")
  counts <- SummarizedExperiment::SummarizedExperiment(
    assays = list(gene = matrix(500L, length(cds), 1L,
                                dimnames = list(names(cds), "RFP_x"))),
    rowRanges = cds,
    colData = S4Vectors::DataFrame(libtype = factor("RFP"), condition = factor("x"),
                                   replicate = "1", row.names = "RFP_x")
  )
  counts <- suppressMessages(simCountTablesRegions(
    counts, regionsToSample = "cds",
    region_proportion = list(cds = list(RFP = 1)), sampling = c(RFP = "MN")
  ))
  list(genome = genome, counts = counts)
}

run_seq_bias_matching_experiment <- function(fixture, ..., fragment_mode = "simulated_rpf") {
  out_dir <- tempfile("reads-"); exp_dir <- tempfile("exp-")
  dir.create(out_dir); dir.create(exp_dir)
  suppressMessages(simNGScoverage(
    fixture$genome, fixture$counts, out_dir = out_dir,
    exp_name = basename(tempfile("exp-")), exp_save_dir = exp_dir,
    fragment_mode = fragment_mode, sampling = list(cds = list(RFP = "DMN")),
    read_lengths_per = list(RFP = 28L), libFormats = list(RFP = "ofst"),
    validate = FALSE, ...
  ))
}

# trace() only instruments the package-internal binding of load_seq_bias(), so
# a plain local `captured <- list()` is unreachable from inside it; .GlobalEnv
# is the one binding both sides can resolve unambiguously.
capture_shifts <- function(expr) {
  assign(".captured_seq_bias_shifts", list(), envir = .GlobalEnv)
  on.exit(rm(".captured_seq_bias_shifts", envir = .GlobalEnv), add = TRUE)
  trace("load_seq_bias", where = asNamespace("coverageSim"), tracer = quote(
    assign(".captured_seq_bias_shifts",
           c(get(".captured_seq_bias_shifts", envir = .GlobalEnv), shift),
           envir = .GlobalEnv)
  ), print = FALSE)
  on.exit(untrace("load_seq_bias", where = asNamespace("coverageSim")), add = TRUE)
  force(expr)
  get(".captured_seq_bias_shifts", envir = .GlobalEnv)
}

test_that("resolve_seq_bias resolves 'AUTO' by fragment_mode and site_reference, and passes everything else through", {
  expect_equal(resolve_seq_bias("AUTO", "simulated_rpf", "a_site"), load_seq_bias(shift = "a-site"))
  expect_equal(resolve_seq_bias("AUTO", "simulated_rpf", "p_site"), load_seq_bias(shift = "p-site"))
  # Non-simulated_rpf modes keep the historical P-site table regardless of site_reference.
  expect_equal(resolve_seq_bias("AUTO", "legacy_point", "a_site"), load_seq_bias(shift = "p-site"))
  expect_equal(resolve_seq_bias("AUTO", "legacy_point", "p_site"), load_seq_bias(shift = "p-site"))

  # Any other value (including NULL, which disables codon bias) passes through unchanged.
  custom <- data.table::data.table(seqs = "M", alpha = 1)
  expect_identical(resolve_seq_bias(custom, "simulated_rpf", "a_site"), custom)
  expect_null(resolve_seq_bias(NULL, "simulated_rpf", "a_site"))
})

test_that("simNGScoverage() defaults seq_bias to match fragment_geometry$site_reference in simulated_rpf mode", {
  fixture <- seq_bias_matching_fixture()
  captured <- capture_shifts({
    run_seq_bias_matching_experiment(fixture)
    run_seq_bias_matching_experiment(fixture, fragment_geometry = list(site_reference = "p_site"))
  })
  expect_equal(captured, list("a-site", "p-site"))

  # trace() only instruments the package-internal binding of load_seq_bias(),
  # so the call below that builds the explicit argument is itself invisible
  # to `captured`; what we're actually checking is that resolve_seq_bias()
  # does not call load_seq_bias() a second time internally for it.
  captured <- capture_shifts(
    run_seq_bias_matching_experiment(
      fixture, seq_bias = load_seq_bias(shift = "p-site"),
      fragment_geometry = list(site_reference = "a_site")
    )
  )
  expect_equal(captured, list())
})

test_that("the deprecated 'physical' fragment_mode alias resolves seq_bias like simulated_rpf", {
  fixture <- seq_bias_matching_fixture()
  captured <- capture_shifts(
    suppressWarnings(run_seq_bias_matching_experiment(
      fixture, fragment_mode = "physical",
      fragment_geometry = list(site_reference = "p_site")
    ))
  )
  expect_equal(captured, list("p-site"))
})

test_that("legacy_point mode keeps the historical P-site default regardless of site_reference", {
  # fragment_geometry$site_reference only affects simulated_rpf's physical
  # fragment placement; legacy_point never builds fragments, so a missing
  # seq_bias there must not change just because fragment_geometry's own
  # default (site_reference = "a_site") happens to differ.
  fixture <- seq_bias_matching_fixture()
  captured <- capture_shifts({
    run_seq_bias_matching_experiment(fixture, fragment_mode = "legacy_point")
    run_seq_bias_matching_experiment(fixture, fragment_mode = "legacy_point",
                                     fragment_geometry = list(site_reference = "p_site"))
  })
  expect_equal(captured, list("p-site", "p-site"))
})

test_that("a wrapper redeclaring its own seq_bias = 'AUTO' default still auto-selects correctly", {
  # Unlike a missing()-based check, identical(seq_bias, "AUTO") composes
  # through wrapper functions: forwarding an explicit "AUTO" still triggers
  # auto-selection, because the check is by value, not by call-site presence.
  wrapper <- function(fixture, ..., seq_bias = "AUTO") {
    run_seq_bias_matching_experiment(fixture, seq_bias = seq_bias, ...)
  }
  fixture <- seq_bias_matching_fixture()
  captured <- capture_shifts({
    wrapper(fixture)
    wrapper(fixture, fragment_geometry = list(site_reference = "p_site"))
  })
  expect_equal(captured, list("a-site", "p-site"))
})

# Run one gene's codon weights through sim_sequence_bias() the way
# simNGScoverage() does (shapes(9) smoothing, default RNase kernel) and read back
# the weight on each codon's first nucleotide, relative to an ordinary codon.
start_region_ratios <- function(alpha, start_codon_marked) {
  n <- length(alpha)
  alpha_matrix <- list(alpha)
  attr(alpha_matrix, "start_codon_marked") <- start_codon_marked
  set.seed(1)
  out <- coverageSim:::sim_sequence_bias(
    NULL, lengths = n * 3L, alpha_matrix = alpha_matrix, seq_acf = shapes(9)
  )[[1L]]
  margin <- (length(out) - 3L * n) / 2L
  per_codon <- out[margin + 1L + seq(0L, by = 3L, length.out = n)]
  per_codon / stats::median(per_codon[-c(1L, 2L, n - 1L, n)])
}

test_that("the start codon keeps its own weight through the smoothing", {
  # In an A-site table the start codon is rare, because initiation places it
  # in the P site, while the next codon carries the initiation peak. Without
  # this, the smoothing filled the start codon from that peak and it came out
  # well above an ordinary codon.
  ordinary <- 0.609
  alpha <- rep(ordinary, 60)
  alpha[1] <- 0.12 * ordinary
  alpha[2] <- 2.24 * ordinary
  restored <- start_region_ratios(alpha, start_codon_marked = TRUE)
  expect_equal(restored[1], 0.12, tolerance = 1e-8)
  # Nothing else moves: the initiation peak and every later codon are exactly
  # what the chain produces without the restoration.
  untouched <- start_region_ratios(alpha, start_codon_marked = FALSE)
  expect_gt(untouched[1], 0.5)
  expect_equal(restored[-1], untouched[-1])
})

test_that("a table without a start-codon weight leaves the first codon to the chain", {
  # Without a ### row the first weight is simply the weight of the codon ATG,
  # which says nothing about the position, so it must not be imposed.
  ordinary <- 0.609
  alpha <- rep(ordinary, 60)
  alpha[1] <- 0.12 * ordinary
  expect_equal(
    start_region_ratios(alpha, start_codon_marked = FALSE),
    start_region_ratios(alpha, start_codon_marked = NULL)
  )
  # Smoothed from its neighbours it ends up well above the 0.12 that would have
  # been imposed.
  expect_gt(start_region_ratios(alpha, start_codon_marked = FALSE)[1], 0.25)
})

test_that("simNGScoverage() tells the chain whether the table marks the start codon", {
  # The two tests above set the marker by hand. This one runs a small
  # simulation and reads what actually reaches sim_sequence_bias(), so a lost
  # marker in add_sequence_bias() fails here.
  capture_marker <- function(seq_bias) {
    assign(".captured_start_marker", list(), envir = .GlobalEnv)
    on.exit(rm(".captured_start_marker", envir = .GlobalEnv), add = TRUE)
    trace("sim_sequence_bias", where = asNamespace("coverageSim"), tracer = quote(
      assign(".captured_start_marker",
             c(get(".captured_start_marker", envir = .GlobalEnv),
               list(attr(alpha_matrix, "start_codon_marked"))),
             envir = .GlobalEnv)
    ), print = FALSE)
    on.exit(untrace("sim_sequence_bias", where = asNamespace("coverageSim")),
            add = TRUE)
    run_seq_bias_matching_experiment(seq_bias_matching_fixture(),
                                     seq_bias = seq_bias)
    get(".captured_start_marker", envir = .GlobalEnv)
  }
  with_marks <- load_seq_bias(type = "codon", shift = "a-site")
  without_marks <- with_marks[!seqs %in% c("###", "%%%", "&&&")]
  expect_true(all(vapply(capture_marker(with_marks), isTRUE, logical(1))))
  expect_false(any(vapply(capture_marker(without_marks), isTRUE, logical(1))))
  expect_gt(length(capture_marker(with_marks)), 0)
})

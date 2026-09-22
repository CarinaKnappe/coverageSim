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

test_that("default_ribosome_site_offset matches the documented P-site/A-site schedule", {
  lengths <- 20:40
  offset_p <- default_ribosome_site_offset(lengths, "p_site")
  offset_a <- default_ribosome_site_offset(lengths, "a_site")
  expect_equal(offset_a, offset_p + 3L)
  expect_equal(offset_p[lengths <= 27L], rep(11L, sum(lengths <= 27L)))
  expect_equal(offset_p[lengths > 27L & lengths <= 30L], rep(12L, sum(lengths > 27L & lengths <= 30L)))
  expect_equal(offset_p[lengths > 30L], rep(13L, sum(lengths > 30L)))
})

test_that("default_fragment_distribution() uses the shared offset schedule", {
  lengths <- c(25L, 27L, 28L, 30L, 33L)
  distribution <- default_fragment_distribution(lengths, "a_site")
  expect_equal(
    distribution$site_offset[match(lengths, distribution$fragment_length)],
    pmin(default_ribosome_site_offset(lengths, "a_site"), lengths - 1L)
  )
})

test_that("learn_fragment_geometry()'s search window is actually centered by the shared offset schedule", {
  # Confirms the call site inside learn_fragment_geometry() (not just
  # default_fragment_distribution()) still delegates to
  # default_ribosome_site_offset(), so the two cannot silently drift apart
  # again if that call site is edited later.
  # trace() runs the tracer inside default_ribosome_site_offset()'s own
  # (package-namespace) scope, so a plain local `captured <- list()` here is
  # unreachable from it via <<-; .GlobalEnv is the one binding both sides can
  # resolve unambiguously (same pattern as test_seq_bias_site_matching.R).
  assign(".captured_offset_calls", list(), envir = .GlobalEnv)
  on.exit(rm(".captured_offset_calls", envir = .GlobalEnv), add = TRUE)
  trace("default_ribosome_site_offset", where = asNamespace("coverageSim"), tracer = quote(
    assign(".captured_offset_calls",
           c(get(".captured_offset_calls", envir = .GlobalEnv),
             list(list(fragment_length, site_reference))),
           envir = .GlobalEnv)
  ), print = FALSE)
  on.exit(untrace("default_ribosome_site_offset", where = asNamespace("coverageSim")), add = TRUE)

  fixture <- make_end_learning_fixture()
  bam <- make_learning_bam(fixture, five = 1, three = 1)
  # The fixture's true site_offset is 15/16 (see helper-learning-fixtures.R);
  # a wide offset_search keeps that inside the scanned window regardless of
  # the p_site guess (12 for these lengths), so the call succeeds.
  geometry <- learn_fragment_geometry(
    bam, fixture$transcripts, fixture$cds,
    site_reference = "p_site", min_length = 28L, max_length = 29L,
    min_reads_per_length = 1L, offset_search = 6L
  )
  captured <- get(".captured_offset_calls", envir = .GlobalEnv)
  expect_true(length(captured) > 0)
  for (call in captured) {
    expect_equal(call[[2]], "p_site")
    expect_true(all(call[[1]] %in% 28:29))
  }
})

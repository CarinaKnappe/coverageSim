test_that("fragment geometry learning recovers length frequencies and A-site offsets", {
  fixture <- make_end_learning_fixture(seed = 707L)
  set.seed(708)
  bam <- make_learning_bam(fixture, five = 1, three = 1)
  geometry <- learn_fragment_geometry(
    bam, fixture$transcripts, fixture$cds,
    min_length = 28L, max_length = 29L,
    min_reads_per_length = 1L, offset_search = 2L
  )
  expect_s3_class(geometry, "covsim_fragment_geometry")
  expect_equal(geometry$source, "learned")
  expect_equal(geometry$site_reference, "a_site")
  expect_equal(geometry$distribution$fragment_length, c(28L, 29L))
  expect_equal(geometry$distribution$site_offset, c(15L, 16L))
  expect_equal(sum(geometry$distribution$probability), 1)
  diagnostics <- attr(geometry, "diagnostics")
  selected <- merge(
    geometry$distribution[, c("fragment_length", "site_offset")],
    diagnostics, by = c("fragment_length", "site_offset")
  )
  expect_true(all(selected$frame_fraction == 1))
  expect_true(all(selected$cds_reads > 0))
  expect_silent(normalize_fragment_geometry(geometry))
})

test_that("fragment geometry learning validates support and search settings", {
  fixture <- make_end_learning_fixture(seed = 711L)
  bam <- make_learning_bam(fixture, five = 1, three = 1)
  expect_error(learn_fragment_geometry(
    bam, fixture$transcripts, fixture$cds,
    min_reads_per_length = 1e9
  ), "reached")
  expect_error(learn_fragment_geometry(
    bam, fixture$transcripts, fixture$cds,
    min_reads_per_length = 1L, offset_search = -1L
  ), "offset_search")
})

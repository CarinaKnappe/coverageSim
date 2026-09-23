test_that("transcript_models() requires a fasta_file and fails fast without one", {
  # build_transcript_models()'s fasta_file defaults to NULL and only
  # conditionally attaches a sequence; transcript_models() itself must still
  # reject a missing fasta_file instead of silently returning sequence-less
  # models, since every one of its callers (simNGScoverage(), learn_end_bias())
  # needs the attached sequence.
  fixture <- make_end_learning_fixture()
  expect_error(
    transcript_models(fixture$transcripts, NULL),
    "fasta_file is required"
  )
})

test_that("transcript_models() and transcript_geometry_models() share coordinates but differ in sequence", {
  fixture <- make_end_learning_fixture()
  with_sequence <- transcript_models(fixture$transcripts, fixture$fasta)
  without_sequence <- transcript_geometry_models(fixture$transcripts)

  expect_true(all(vapply(
    with_sequence, function(m) is.character(m$sequence) && nzchar(m$sequence), logical(1)
  )))
  expect_true(all(vapply(without_sequence, function(m) is.null(m$sequence), logical(1))))

  shared_fields <- c("transcript_id", "exons", "cumulative_start", "length", "strand")
  expect_equal(
    lapply(with_sequence, `[`, shared_fields),
    lapply(without_sequence, `[`, shared_fields)
  )
})

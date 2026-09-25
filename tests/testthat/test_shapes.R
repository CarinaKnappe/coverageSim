test_that("shapes() only builds the working autocor_window quote", {
  # i = 1..4 used to build a cosine-wave formula (bquote(abs(cos(x))^...))
  # referencing variables (x, b, l) that are never bound where
  # sim_sequence_bias() evaluates it -- only alpha_vec is -- so it crashed
  # with "object 'x' not found" as soon as it was actually used, even
  # though shapes() itself never errored (building the quote doesn't
  # evaluate it). Confirmed identical to upstream/master: this was broken
  # from the original implementation, not something introduced later.
  # Those branches are now removed entirely: every i other than 0 routes to
  # the same working mechanism (a real, manuscript-validated tRNA/wobble
  # neighbor-correlation smoothing), just with i as the window size.
  expect_null(shapes(0))
  for (i in c(1, 2, 3, 4, 5, 9, 15)) {
    expr <- shapes(i)
    expect_true(is.call(expr))
    alpha_vec <- c(1, 5, 2, 8, 1, 6, 3, 9, 2, 4)
    result <- eval(expr)
    expect_true(is.numeric(result))
    expect_true(all(is.finite(result)))
    expect_equal(length(result), length(alpha_vec))
  }
})

test_that("shapes() defaults to the same window size simNGScoverage() itself uses", {
  expect_equal(shapes(), shapes(9))
})

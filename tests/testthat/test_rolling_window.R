test_that("rolling_weighted_sum applies centred weights in their given order", {
  x <- c(0, 0, 1:6, 0, 0)
  weights <- c(2, 1, 3)
  # Position 2: 2 * 0 + 1 * 0 + 3 * 1 = 3; position 3: 2 * 0 + 1 * 1 + 3 * 2 = 7.
  expect_equal(rolling_weighted_sum(x, weights),
               c(NA, 3, 7, 13, 19, 25, 31, 16, 12, NA))
})

test_that("rolling_weighted_sum matches a per-window sum, including NAs", {
  set.seed(11)
  x <- stats::rexp(200)
  x[c(5, 90, 91)] <- NA
  weights <- stats::runif(7)
  expected <- data.table::frollapply(x, 7L, function(i) sum(i * weights),
                                     align = "center", fill = NA)
  expect_equal(rolling_weighted_sum(x, weights), expected)
})

test_that("rolling_weighted_sum handles vectors shorter than the window", {
  expect_identical(rolling_weighted_sum(c(1, 2), c(1, 1, 1)), c(NA_real_, NA_real_))
  expect_identical(rolling_weighted_sum(numeric(), 1), numeric())
})

test_that("rolling_weighted_sum refuses an even number of weights", {
  expect_error(rolling_weighted_sum(1:5, c(1, 1)), "odd number of weights")
})

test_that("autocor_window equals the mean of x / distance penalty per window", {
  per_window_mean <- function(vec, max.lag, penalty = 1.5) {
    window <- max.lag * 2 + 1
    half <- ceiling(window / 2)
    penalty_weights <- c(seq(half, 2), 1, seq(2, half))^penalty
    padding <- rep(0, floor(window / 2) * 2)
    data.table::frollapply(c(padding, vec, padding), window,
                           function(x) mean(x / penalty_weights, na.rm = TRUE),
                           align = "center", fill = NA)
  }
  set.seed(12)
  vec <- stats::rexp(60)
  expect_equal(autocor_window(vec, 3), per_window_mean(vec, 3))
  expect_equal(autocor_window(vec, 9), per_window_mean(vec, 9))
  # Missing values are left out of the mean, as with mean(na.rm = TRUE).
  vec[c(4, 30, 31)] <- NA
  expect_equal(autocor_window(vec, 3), per_window_mean(vec, 3))
  # padding.rm keeps only the computed inner positions.
  expect_equal(autocor_window(vec, 3, padding.rm = TRUE),
               per_window_mean(vec, 3)[-c(1:6, 67:72)])
})

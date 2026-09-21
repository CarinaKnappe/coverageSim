test_that("rgdirmnom defaults to extraDistr::rdirmnom exactly", {
  alpha <- matrix(c(5, 3, 1, 1e-24, 2, 2, 2, 2), nrow = 2, byrow = TRUE)
  size <- c(10L, 20L)

  set.seed(1)
  observed <- rgdirmnom(2, size = size, alpha = alpha)
  set.seed(1)
  expected <- extraDistr::rdirmnom(n = 2, size = size, alpha = alpha)

  expect_equal(observed, expected)
})

test_that("row sums are preserved -- dirichlet family, zero-inflation on", {
  set.seed(2)
  alpha <- matrix(c(5, 3, 1, 4, 1e-24, 1e-24,
                    2, 2, 2, 2, 2, 2), nrow = 2, byrow = TRUE)
  size <- c(37L, 51L)
  observed <- rgdirmnom(2, size = size, alpha = alpha, zero_inflation = 0.7)

  expect_equal(unname(rowSums(observed)), size)
})

test_that("row sums are preserved -- generalized family, various gdm_scale", {
  set.seed(3)
  alpha <- matrix(c(5, 3, 1, 4, 6, 2,
                    2, 2, 2, 2, 2, 2), nrow = 2, byrow = TRUE)
  size <- c(40L, 100L)

  for (scale in c(0.3, 1, 4)) {
    observed <- rgdirmnom(2, size = size, alpha = alpha,
                          family = "generalized", gdm_scale = scale)
    expect_equal(unname(rowSums(observed)), size)
  }
})

test_that("row sums are preserved -- generalized family with zero-inflation and ragged padding", {
  set.seed(4)
  alpha <- matrix(c(5, 3, 1, 1e-24, 1e-24, 1e-24,
                    2, 2, 2, 2, 2, 2), nrow = 2, byrow = TRUE)
  size <- c(15L, 60L)
  observed <- rgdirmnom(2, size = size, alpha = alpha, family = "generalized",
                        gdm_scale = 1.7, zero_inflation = 0.5)

  expect_equal(unname(rowSums(observed)), size)
  # The padded tail of the shorter gene must never receive reads.
  expect_equal(observed[1, 4:6], c(0, 0, 0))
})

test_that("zero_inflation increases the observed proportion of zero counts", {
  set.seed(5)
  alpha <- matrix(rep(2, 20), nrow = 1)
  size <- 200L

  reps <- 300
  zero_share <- function(zi) {
    draws <- do.call(rbind, lapply(seq_len(reps), function(i) {
      rgdirmnom(1, size = size, alpha = alpha, zero_inflation = zi)
    }))
    mean(draws == 0)
  }

  off <- zero_share(0)
  on <- zero_share(0.85)
  expect_gt(on, off)
})

test_that("generalized family at gdm_scale = 1 reproduces standard DMN moments", {
  set.seed(6)
  alpha <- c(8, 4, 2, 1, 1)
  alpha_matrix <- matrix(alpha, nrow = 1)
  size <- 500L
  reps <- 2000

  draws_std <- extraDistr::rdirmnom(n = reps, size = rep(size, reps),
                                    alpha = matrix(rep(alpha, each = reps), nrow = reps))
  draws_gen <- do.call(rbind, lapply(seq_len(reps), function(i) {
    rgdirmnom(1, size = size, alpha = alpha_matrix, family = "generalized", gdm_scale = 1)
  }))

  expected_mean <- size * alpha / sum(alpha)
  expect_equal(colMeans(draws_gen), expected_mean, tolerance = 0.1)
  # Loose distributional sanity check against the standard DMN's own draws
  # (same theoretical distribution at gdm_scale = 1, compared via SDs since
  # RNG streams differ between the two algorithms).
  expect_equal(apply(draws_gen, 2, sd), apply(draws_std, 2, sd), tolerance = 0.35)
})

test_that("all-eligible-positions-excluded guard keeps at least one live position", {
  set.seed(7)
  alpha <- matrix(c(1, 1, 1e-24), nrow = 1)
  size <- 30L

  for (i in 1:20) {
    observed <- rgdirmnom(1, size = size, alpha = alpha, zero_inflation = 0.999)
    expect_equal(sum(observed), size)
    expect_equal(observed[1, 3], 0)
  }
})

test_that("rgdirmnom validates its inputs", {
  alpha <- matrix(c(1, 1, 1, 1), nrow = 1)
  expect_error(rgdirmnom(1, size = c(1, 2), alpha = alpha))
  expect_error(rgdirmnom(1, size = 10, alpha = alpha, zero_inflation = 1))
  expect_error(rgdirmnom(1, size = 10, alpha = alpha, zero_inflation = -0.1))
  expect_error(rgdirmnom(1, size = 10, alpha = alpha, gdm_scale = 0))
  expect_error(rgdirmnom(1, size = 10, alpha = alpha, gdm_scale = -1))
  expect_error(rgdirmnom(1, size = 10, alpha = alpha, family = "nope"))
})

test_that("apply_zero_inflation_mask_vector matches the matrix version row-for-row", {
  alpha <- c(5, 3, 1, 4, 1e-24, 1e-24)
  set.seed(10)
  vector_result <- apply_zero_inflation_mask_vector(alpha, 0.6)
  set.seed(10)
  matrix_result <- apply_zero_inflation_mask(matrix(alpha, nrow = 1), 0.6, 1e-24)

  expect_equal(vector_result, as.vector(matrix_result))
  # Genuine tail padding is untouched regardless of the random draw.
  expect_equal(vector_result[5:6], c(1e-24, 1e-24))
})

test_that("rgdirichlet_generalized returns a valid probability vector and reduces to Dirichlet at gdm_scale = 1", {
  set.seed(11)
  alpha <- c(8, 4, 2, 1, 1)
  reps <- 2000
  draws <- do.call(rbind, lapply(seq_len(reps), function(i) {
    rgdirichlet_generalized(alpha, gdm_scale = 1)
  }))

  expect_true(all(draws >= 0))
  expect_equal(rowSums(draws), rep(1, reps))
  expect_equal(colMeans(draws), alpha / sum(alpha), tolerance = 0.08)
})

test_that("rgdirichlet_generalized treats padded/masked entries as structurally absent", {
  alpha <- c(3, 2, 1e-24, 1e-24)
  set.seed(12)
  for (i in 1:20) {
    p <- rgdirichlet_generalized(alpha, gdm_scale = 1.4)
    expect_equal(sum(p), 1, tolerance = 1e-8)
    expect_equal(p[3:4], c(0, 0))
  }
})

test_that("draw_site_probabilities defaults are unchanged by the new ZIGDM arguments", {
  weights <- c(5, 3, 1, 4, 2)
  set.seed(13)
  observed <- draw_site_probabilities(weights, dirichlet = TRUE)
  set.seed(13)
  explicit_defaults <- draw_site_probabilities(
    weights, dirichlet = TRUE,
    family = "dirichlet", gdm_scale = 1, zero_inflation = 0
  )

  expect_equal(observed, explicit_defaults)
  expect_equal(sum(observed), 1, tolerance = 1e-8)
})

test_that("draw_site_probabilities: zero_inflation and family = 'generalized' work end-to-end", {
  weights <- c(5, 3, 1, 4, 2, 6, 1, 3)
  set.seed(14)
  zi_result <- draw_site_probabilities(weights, dirichlet = TRUE, zero_inflation = 0.8)
  expect_equal(sum(zi_result), 1, tolerance = 1e-8)
  expect_true(sum(zi_result < 1e-6) >= 1)

  set.seed(15)
  gen_result <- draw_site_probabilities(weights, dirichlet = TRUE,
                                        family = "generalized", gdm_scale = 2)
  expect_equal(sum(gen_result), 1, tolerance = 1e-8)
  expect_true(all(gen_result >= 0))

  # dirichlet = FALSE (the MN sampling path) is untouched by the new
  # arguments -- it never even looks at them.
  expect_equal(draw_site_probabilities(weights, zero_inflation = 0.9),
              weights / sum(weights))
})

test_that("last live column absorbs all remaining reads when trailing columns are masked", {
  set.seed(16)
  alpha <- matrix(c(4, 1e-24, 3, 1e-24, 1e-24,
                    1e-24, 2, 1e-24, 1e-24, 1e-24,
                    1, 2, 3, 4, 5), nrow = 3, byrow = TRUE)
  size <- c(500000L, 123456L, 1000L)
  for (scale in c(0.2, 1, 5)) {
    observed <- rgdirmnom(3, size = size, alpha = alpha,
                          family = "generalized", gdm_scale = scale)
    expect_equal(unname(rowSums(observed)), size)
    expect_equal(observed[1, c(2, 4, 5)], c(0, 0, 0))
    expect_equal(observed[2, c(1, 3, 4, 5)], c(0, 0, 0, 0))
    expect_equal(observed[2, 2], 123456)
  }

  p <- rgdirichlet_generalized(c(4, 1e-24, 3, 1e-24, 1e-24), gdm_scale = 2)
  expect_equal(sum(p), 1)
  expect_equal(p[c(2, 4, 5)], c(0, 0, 0))
})

test_that("last_live_column finds the last column above pad_value per row", {
  alpha <- matrix(c(4, 1e-24, 3, 1e-24,
                    1e-24, 1e-24, 1e-24, 1e-24,
                    1, 1, 1, 1), nrow = 3, byrow = TRUE)
  expect_equal(last_live_column(alpha, 1e-24), c(3L, 0L, 4L))
})

test_that("rgdirmnom rejects invalid sizes, alphas and rows without any live position", {
  alpha <- matrix(c(2, 3, 1, 1), nrow = 1)
  expect_error(rgdirmnom(1, size = -3, alpha = alpha, family = "generalized"))
  expect_error(rgdirmnom(1, size = 2.5, alpha = alpha))
  expect_error(rgdirmnom(1, size = NA_real_, alpha = alpha))
  expect_error(rgdirmnom(1, size = 10, alpha = -alpha))
  expect_error(rgdirmnom(1, size = 10, alpha = matrix(c(2, 0, 3), nrow = 1)))
  expect_error(rgdirmnom(1, size = 10, alpha = alpha * Inf))
  expect_error(rgdirmnom(1, size = 10, alpha = alpha, pad_value = 0))

  all_padded <- matrix(c(1e-24, 1e-24), nrow = 1)
  expect_error(rgdirmnom(1, size = 10, alpha = all_padded, family = "generalized"),
               "at least one position")
  # Nothing to distribute: a zero-size row without live positions is allowed.
  expect_equal(unname(rgdirmnom(1, size = 0, alpha = all_padded,
                                family = "generalized")), matrix(c(0, 0), nrow = 1))
})

test_that("vectorised rgdirichlet_generalized handles long vectors and single positions", {
  set.seed(17)
  long <- rgdirichlet_generalized(rep(0.5, 50000), gdm_scale = 1.3)
  expect_equal(sum(long), 1, tolerance = 1e-8)
  expect_true(all(long >= 0))
  expect_equal(rgdirichlet_generalized(2, gdm_scale = 1), 1)
})

test_that("gdm_scale != 1 shifts the expected proportions as documented", {
  set.seed(18)
  alpha <- c(1, 1, 1)
  reps <- 4000
  # For alpha = (1, 1, 1): E[p1] = a1 / (a1 + gdm_scale * (a2 + a3)).
  for (scale in c(0.5, 3)) {
    expected_p1 <- 1 / (1 + scale * 2)
    counts <- rgdirmnom(reps, size = rep(100L, reps),
                        alpha = matrix(alpha, nrow = reps, ncol = 3, byrow = TRUE),
                        family = "generalized", gdm_scale = scale)
    expect_equal(mean(counts[, 1]) / 100, expected_p1, tolerance = 0.05)

    p1 <- vapply(seq_len(reps), function(i) {
      rgdirichlet_generalized(alpha, gdm_scale = scale)[1]
    }, numeric(1))
    expect_equal(mean(p1), expected_p1, tolerance = 0.05)
  }
})

test_that("rgdirichlet_generalized rejects vectors without any live position", {
  expect_error(rgdirichlet_generalized(c(1e-24, 1e-24), gdm_scale = 1), "at least one position")
  expect_error(draw_site_probabilities(c(1e-24, 1e-24), dirichlet = TRUE, family = "generalized"),
               "at least one position")
})

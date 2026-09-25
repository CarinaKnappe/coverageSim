test_that("frollapply_compat applies a centered rolling function", {
  x <- 1:7
  observed <- frollapply_compat(
    x,
    window = 3L,
    FUN = sum,
    align = "center",
    fill = NA
  )

  expect_equal(observed, c(NA, 6, 9, 12, 15, 18, NA))
})

test_that("frollapply_compat surfaces FUN's own error instead of masking it", {
  # The previous implementation retried 4 argument-name combinations on any
  # error, including one genuinely thrown by FUN itself; since all 4
  # retries failed identically, only the last (a confusing, unrelated
  # "unused argument"-style message) ever reached the caller, never FUN's
  # actual error.
  broken_fun <- function(window_values) stop("deliberately broken FUN")
  expect_error(
    frollapply_compat(1:7, window = 3L, FUN = broken_fun, align = "center", fill = NA),
    "deliberately broken FUN"
  )
})

test_that("frollapply_compat does not warn about a deprecated argument name", {
  expect_no_warning(
    frollapply_compat(1:7, window = 3L, FUN = sum, align = "center", fill = NA)
  )
})

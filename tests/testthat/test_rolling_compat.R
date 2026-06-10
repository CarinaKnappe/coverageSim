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

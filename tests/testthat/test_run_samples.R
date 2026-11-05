test_that("resolve_cores accepts whole numbers of at least 1 only", {
  expect_identical(resolve_cores(1), 1L)
  for (bad in list(0, -1, 1.5, NA_real_, c(1, 2), "2")) {
    expect_error(resolve_cores(bad), "single whole number")
  }
})

test_that("run_samples gives the same draws and random state for any core count", {
  simulate_sample <- function(s) c(s, stats::runif(3))
  seeds <- c(11L, 22L, 33L)
  set.seed(1)
  serial <- run_samples(simulate_sample, seeds, cores = 1L)
  state_serial <- .Random.seed
  set.seed(1)
  skip_on_os("windows")
  parallel <- run_samples(simulate_sample, seeds, cores = 2L)
  expect_identical(parallel, serial)
  expect_identical(.Random.seed, state_serial)
  # Each sample draws from its own seed, not from the run's random state.
  set.seed(22L)
  expect_identical(serial[[2]], c(2, stats::runif(3)))
})

test_that("run_samples passes on warnings and errors of every sample", {
  warn_sample <- function(s) {
    warning("sample ", s, " warns")
    s
  }
  skip_on_os("windows")
  for (cores in c(1L, 2L)) {
    expect_warning(
      expect_warning(result <- run_samples(warn_sample, 1:2, cores), "sample 1 warns"),
      "sample 2 warns"
    )
    expect_identical(result, list(1L, 2L))
    expect_error(
      run_samples(function(s) if (s == 2) stop("sample 2 broke") else s, 1:2, cores),
      "sample 2 broke"
    )
  }
})

test_that("run_samples passes on warnings also when a sample fails", {
  skip_on_os("windows")
  failing <- function(s) {
    warning("diagnostic of sample ", s)
    if (s == 2) stop("sample 2 broke")
    s
  }
  for (cores in c(1L, 2L)) {
    expect_warning(
      expect_warning(
        expect_error(run_samples(failing, 1:2, cores), "sample 2 broke"),
        "diagnostic of sample 1"
      ),
      "diagnostic of sample 2"
    )
  }
})

test_that("autocor_window needs a lag of at least 1", {
  expect_error(autocor_window(1:3, 0), "at least 1")
})

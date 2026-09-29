test_that("current_code_version() reports the saved state of a clean repository", {
  skip_if(nchar(Sys.which("git")) == 0, "git is not available")
  repo <- tempfile("code-version-clean-")
  dir.create(repo)
  writeLines("x <- 1", file.path(repo, "a.R"))
  git <- function(...) system2("git", c("-C", shQuote(repo), ...),
                              stdout = FALSE, stderr = FALSE)
  git("init", "--quiet")
  git("config", "user.email", "test@example.com")
  git("config", "user.name", "Test")
  git("add", "a.R")
  git("commit", "--quiet", "-m", "first")

  version <- current_code_version(repo)

  expect_equal(nrow(version), 1L)
  expect_match(version$commit, "^[0-9a-f]{40}$")
  expect_equal(version$commit_short, substr(version$commit, 1L, 7L))
  expect_match(version$commit_date, "^\\d{4}-\\d{2}-\\d{2}$")
  expect_false(version$unsaved_changes)
  expect_equal(version$unsaved_files, "")
  expect_true(version$fully_reproducible)
})

test_that("current_code_version() names the files that were never saved", {
  skip_if(nchar(Sys.which("git")) == 0, "git is not available")
  repo <- tempfile("code-version-dirty-")
  dir.create(repo)
  writeLines("x <- 1", file.path(repo, "a.R"))
  git <- function(...) system2("git", c("-C", shQuote(repo), ...),
                              stdout = FALSE, stderr = FALSE)
  git("init", "--quiet")
  git("config", "user.email", "test@example.com")
  git("config", "user.name", "Test")
  git("add", "a.R")
  git("commit", "--quiet", "-m", "first")
  # An edit that was never saved, and a file that was never added at all: both
  # mean the run used code that the project history does not contain.
  writeLines("x <- 2", file.path(repo, "a.R"))
  writeLines("y <- 3", file.path(repo, "b.R"))

  version <- current_code_version(repo)

  expect_true(version$unsaved_changes)
  expect_false(version$fully_reproducible)
  expect_true(grepl("a.R", version$unsaved_files, fixed = TRUE))
  expect_true(grepl("b.R", version$unsaved_files, fixed = TRUE))
})

test_that("current_code_version() says so plainly outside version control", {
  plain <- tempfile("code-version-none-")
  dir.create(plain)

  version <- current_code_version(plain)

  expect_true(is.na(version$commit))
  expect_false(version$fully_reproducible)
  expect_equal(nrow(version), 1L)
})

test_that("write_code_version() writes a readable record next to the output", {
  skip_if(nchar(Sys.which("git")) == 0, "git is not available")
  out <- file.path(tempfile("code-version-out-"), "reads")

  version <- write_code_version(out, repo_dir = tempfile("no-repo-"))

  written <- data.table::fread(file.path(out, "code_version.tsv"))
  expect_equal(nrow(written), 1L)
  expect_true("fully_reproducible" %in% names(written))
  expect_equal(names(written), names(version))
})

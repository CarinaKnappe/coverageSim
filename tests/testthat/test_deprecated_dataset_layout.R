test_that("deprecated human-genome-only dataset is clearly separated", {
  data_root <- testthat::test_path(
    "..", "tests_and_demos", "real_human_riboseq"
  )

  testthat::skip_if_not(dir.exists(data_root), "Large demo dataset is unavailable")

  deprecated_root <- file.path(data_root, "deprecated")
  canonical_root <- file.path(
    data_root,
    "coverageSim_from_human_genome_only",
    "human_genome_only_v2"
  )

  expect_false(dir.exists(file.path(data_root, "human_genome_only_covsim")))
  expect_true(dir.exists(file.path(deprecated_root, "human_genome_only_covsim")))
  expect_true(file.exists(file.path(deprecated_root, "README.md")))
  expect_true(dir.exists(canonical_root))
  expect_true(file.exists(file.path(data_root, "INFO.txt")))

  marker <- readLines(file.path(deprecated_root, "README.md"), warn = FALSE)
  expect_true(any(grepl("Do not use", marker, fixed = TRUE)))
  expect_true(any(grepl("human_genome_only_v2", marker, fixed = TRUE)))

  overview <- readLines(file.path(data_root, "INFO.txt"), warn = FALSE)
  expect_true(any(grepl("coverageSim_from_real_input", overview, fixed = TRUE)))
  expect_true(any(grepl("coverageSim_from_human_genome_only", overview, fixed = TRUE)))
  expect_true(any(grepl("SRR32491292.sorted.bam", overview, fixed = TRUE)))
})

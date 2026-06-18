repo_dir <- normalizePath(
  Sys.getenv("COVSIM_REPO", unset = testthat::test_path("..", "..")),
  mustWork = TRUE
)

source(file.path(repo_dir, "rust_helpers", "R", "G_end_dropout_helpers.R"))

test_that("G-end dropout labels are stable", {
  expect_equal(g_end_dropout_label(0.30), "drop_G_30")
  expect_equal(g_end_dropout_label(0.60), "drop_G_60")
  expect_equal(g_end_dropout_label(1.00), "drop_G_100")
})

test_that("experiment paths are rewritten to the dropout dataset", {
  lines <- c(
    '"fasta","old/base/genome/genome.fa"',
    '"RFP","","1","WT","","old/base/reads/RFP_WT_1.bam"'
  )

  rewritten <- rewrite_experiment_base(lines, "old/base", "new/base")

  expect_equal(
    rewritten,
    c(
      '"fasta","new/base/genome/genome.fa"',
      '"RFP","","1","WT","","new/base/reads/RFP_WT_1.bam"'
    )
  )
})

test_that("G-end AWK program removes reads ending in G", {
  skip_if(Sys.which("awk") == "", "awk is not available")

  sam_file <- tempfile(fileext = ".sam")
  stats_file <- tempfile(fileext = ".txt")
  writeLines(
    c(
      "@HD\tVN:1.6\tSO:coordinate",
      "r1\t0\tchr1\t1\t255\t3M\t*\t0\t0\tAAG\tIII",
      "r2\t0\tchr1\t2\t255\t3M\t*\t0\t0\tAAA\tIII",
      "r3\t0\tchr1\t3\t255\t3M\t*\t0\t0\tCCG\tIII",
      "r4\t0\tchr1\t4\t255\t3M\t*\t0\t0\tTTT\tIII"
    ),
    sam_file
  )

  awk_file <- tempfile(fileext = ".awk")
  writeLines(g_end_dropout_awk_program(), awk_file)

  filtered <- system2(
    "awk",
    c(
      "-v", "drop=1",
      "-v", "seed=1",
      "-v", paste0("stats=", stats_file),
      "-f", awk_file,
      sam_file
    ),
    stdout = TRUE
  )

  expect_true(any(grepl("^@HD", filtered)))
  expect_false(any(grepl("^r1\\t", filtered)))
  expect_true(any(grepl("^r2\\t", filtered)))
  expect_false(any(grepl("^r3\\t", filtered)))
  expect_true(any(grepl("^r4\\t", filtered)))

  stats <- readLines(stats_file)
  expect_true("total\t4" %in% stats)
  expect_true("g_end\t2" %in% stats)
  expect_true("removed\t2" %in% stats)
  expect_true("kept\t2" %in% stats)
})

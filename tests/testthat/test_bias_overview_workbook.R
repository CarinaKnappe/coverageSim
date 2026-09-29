test_that("bias overview workbook documents configurable and fixed biases", {
  repo_dir <- normalizePath(testthat::test_path("..", ".."))
  data_root <- Sys.getenv(
    "COVSIM_DATA_ROOT",
    unset = file.path(dirname(repo_dir), "coverageSim_data")
  )
  workbook <- file.path(data_root, "references", "bias", "coverageSim_bias_overview.xlsx")
  testthat::skip_if_not(
    file.exists(workbook),
    "Local bias workbook is stored outside the repository"
  )
  expect_true(file.exists(workbook))

  archive <- utils::unzip(workbook, list = TRUE)
  expect_true(all(c(
    "xl/workbook.xml",
    "xl/worksheets/sheet1.xml",
    "xl/worksheets/sheet2.xml",
    "xl/worksheets/sheet3.xml"
  ) %in% archive$Name))

  extraction_dir <- tempfile("bias-overview-")
  dir.create(extraction_dir)
  shared_strings <- intersect("xl/sharedStrings.xml", archive$Name)
  utils::unzip(
    workbook,
    files = c(
      "xl/workbook.xml",
      "xl/worksheets/sheet1.xml",
      "xl/worksheets/sheet2.xml",
      "xl/worksheets/sheet3.xml",
      shared_strings
    ),
    exdir = extraction_dir
  )
  workbook_xml <- paste(readLines(
    file.path(extraction_dir, "xl", "workbook.xml"), warn = FALSE
  ), collapse = " ")
  sheet_xml <- paste(vapply(seq_len(3L), function(i) {
    paste(readLines(file.path(
      extraction_dir, "xl", "worksheets", paste0("sheet", i, ".xml")
    ), warn = FALSE), collapse = " ")
  }, character(1)), collapse = " ")
  if (length(shared_strings)) {
    sheet_xml <- paste(sheet_xml, paste(readLines(
      file.path(extraction_dir, shared_strings), warn = FALSE
    ), collapse = " "))
  }

  expect_match(workbook_xml, "Manuscript Steps", fixed = TRUE)
  expect_match(workbook_xml, "Bias Overview", fixed = TRUE)
  expect_match(workbook_xml, "Development Plan", fixed = TRUE)
  expect_match(sheet_xml, "User-tunable?", fixed = TRUE)
  expect_match(sheet_xml, "Hard-coded or implicit details", fixed = TRUE)
  expect_match(sheet_xml, "Current concern / interpretation", fixed = TRUE)
  expect_match(sheet_xml, "stop-codon peak", fixed = TRUE)
  expect_match(sheet_xml, "five_prime_bias", fixed = TRUE)
  expect_match(sheet_xml, "A-site fragment geometry", fixed = TRUE)
  expect_match(sheet_xml, "motif-wise median", fixed = TRUE)
  expect_match(sheet_xml, "R1-R10", fixed = TRUE)
  expect_false(grepl("effectively selected R1", sheet_xml, fixed = TRUE))
})

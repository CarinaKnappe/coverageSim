test_that("load_seq_bias returns the default AA p-site bias table", {
  dt <- load_seq_bias()
  expect_s3_class(dt, "data.table")
  expect_equal(unique(dt$variable), "R2")
  expect_true(all(c("seqs", "alpha") %in% colnames(dt)))
})

test_that("load_seq_bias supports alternate built-in bias tables", {
  aa_stop <- load_seq_bias(bias = "stop_codon")
  codon_all <- load_seq_bias(type = "codon", shift = "a-site", bias = "all")

  expect_equal(unique(aa_stop$variable), "R1")
  expect_true(length(unique(codon_all$variable)) > 1)
})

test_that("load_seq_bias rejects unsupported arguments", {
  expect_error(load_seq_bias(type = "peptide"))
  expect_error(load_seq_bias(shift = "e-site"))
  expect_error(load_seq_bias(bias = "unknown"))
})

context("Full simulations")
library(coverageSim)


test_that("Default genome simulation works", {
  simGenome6 <- simGenome(n = 6, genome_name = "artificial_uorf_1",
                          max_uorfs = 0, cds_length = c(rep.int(330, 3), rep.int(300, 3)),
                          export_txdb = TRUE)
})

test_that("uORF (n=1) genome simulation works", {
  simGenome_uORFs_1 <- simGenome(n = 6, genome_name = "artificial_uorf_1",
                                 max_uorfs = 1, cds_length = c(rep.int(330, 3), rep.int(300, 3)),
                                 export_txdb = TRUE)
})

test_that("uORF (n=2) genome simulation works", {
  simGenome_uORFs_2 <- simGenome(n = 6, genome_name = "artificial_uorf_2",
                                 max_uorfs = 2, cds_length = c(rep.int(330, 3), rep.int(300, 3)),
                                 export_txdb = TRUE)
})

test_that("Default gene count table simulation (RFP) works", {
  gene_count_table <- simCountTables(loadRegion(simGenome6["txdb"], "cds"),
                                    libtypes = "RFP", print_statistics = FALSE,
                                    interceptMean = 10)
})

test_that("Default gene count table simulation (RFP & uORFs) works", {
  gene_count_table_uorfs <- simCountTables(loadRegion(simGenome_uORFs_1["txdb"], "cds"),
                                    libtypes = "RFP", print_statistics = FALSE,
                                    interceptMean = 10)
})

test_that("Default gene count table simulation (RFP & uORFs #2) works", {
  gene_count_table_uorfs_2 <- simCountTables(loadRegion(simGenome_uORFs_2["txdb"], "cds"),
                                           libtypes = "RFP", print_statistics = FALSE,
                                           interceptMean = 10)
})

test_that("Default gene region count table simulation (RFP) works", {
  region_count_table <- simCountTablesRegions(gene_count_table,
                                              regionsToSample = c("cds"))
})

test_that("Default gene region count table simulation (RFP & uORFs) works", {
  region_count_table_uorfs <- simCountTablesRegions(gene_count_table_uorfs,
                                              regionsToSample = c("cds", "uorf"))
})

test_that("Default gene region count table simulation (RFP & uORFs #2) works", {
  region_count_table_uorfs_2 <- simCountTablesRegions(gene_count_table_uorfs_2,
                                                    regionsToSample = c("cds", "uorf"))
})

test_that("Default read coverage simulation (RFP) works", {
  df_default <- simNGScoverage(simGenome6, region_count_table[,1],
                 exp_name = "aritifical_default",
                 exp_save_dir = tempdir())
})

test_that("Default read coverage simulation (RFP & uORFs) works", {
  df_uorfs <- simNGScoverage(simGenome_uORFs_1, region_count_table_uorfs[,1],
                 exp_name = "aritifical_uorfs",
                 exp_save_dir = tempdir())
})

test_that("Default read coverage simulation (RFP & uORFs #2) works", {
  df_uorfs2 <- simNGScoverage(simGenome_uORFs_2, region_count_table_uorfs_2[,1],
                             exp_name = "aritifical_uorfs_2",
                             exp_save_dir = tempdir(), debug_coverage = T)
})

test_that("Loading simulated library (RFP & uORFs) works", {
  gr <- fimport(filepath(df_uorfs, "default"))
  expect_is(gr, "GRanges")
})

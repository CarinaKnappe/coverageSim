quiet_count_tables <- function(...) {
  suppressMessages(simCountTables(..., print_statistics = FALSE, plot_PCA = FALSE))
}

test_that("simCountTables works for every library type and condition layout", {
  set.seed(1)
  layouts <- list(
    list(libtypes = c("RFP", "RNA"), conditions = c("WT", "Mutant"), replicates = 2),
    list(libtypes = "RFP", conditions = c("A", "B", "C"), replicates = 2),
    list(libtypes = c("RFP", "RNA"), conditions = c("A", "B", "C"), replicates = 2),
    list(libtypes = "RFP", conditions = "WT", replicates = 3),
    list(libtypes = c("RFP", "RNA"), conditions = "WT", replicates = 2),
    list(libtypes = "RFP", conditions = "WT", replicates = 1)
  )
  for (layout in layouts) {
    counts <- do.call(quiet_count_tables, c(list(n = 15), layout))
    expected_samples <- length(layout$libtypes) * length(layout$conditions) * layout$replicates
    expect_equal(dim(counts), c(15L, expected_samples))
    truth <- as.data.frame(S4Vectors::mcols(counts))
    expect_true(all(paste0("trueDisp_", layout$libtypes) %in% names(truth)))
    expect_equal(sum(grepl("^trueBeta", names(truth))), max(1L, length(layout$conditions) - 1L))
  }
})

test_that("simCountTablesRegions() works for a single-gene count table", {
  # mat[, colnames(mat) == region] (Sim_count_tables.R) used to drop to a
  # plain vector for exactly one gene (R's default `[` drop = TRUE), which
  # as.matrix() then rebuilt as a column matrix (samples x 1) instead of the
  # intended (1 gene x samples) -- colnames<- then failed with "length of
  # 'dimnames' [2] not equal to array extent".
  set.seed(1)
  counts <- quiet_count_tables(
    n = 1, libtypes = c("RFP", "RNA"), conditions = "WT", replicates = 1
  )
  region_counts <- suppressMessages(simCountTablesRegions(
    counts, regionsToSample = c("leader", "cds")
  ))
  expect_equal(dim(SummarizedExperiment::assay(region_counts, "cds")), c(1L, 2L))
  expect_equal(dim(SummarizedExperiment::assay(region_counts, "leader")), c(1L, 2L))
  expect_equal(colnames(SummarizedExperiment::assay(region_counts, "cds")),
               colnames(SummarizedExperiment::assay(counts)))
})

test_that("every sample uses the dispersion of its own library type", {
  set.seed(2)
  # RFP gets an almost Poisson dispersion, RNA a very large one.
  by_library <- function(x) {
    x[, 1] <- 1e-6
    x[, 2] <- 100
    x
  }
  counts <- quiet_count_tables(
    n = 3000, libtypes = c("RFP", "RNA"), conditions = c("WT", "Mutant"),
    replicates = 2, interceptMean = 8, interceptSD = 0.5,
    betaSD = 0, betaLibSD = c(RFP = 0, RNA = 0), dispMeanRel = by_library
  )
  matrix_counts <- SummarizedExperiment::assay(counts)
  spread <- function(a, b) stats::var(log((matrix_counts[, a] + 1) / (matrix_counts[, b] + 1)))
  # Replicates of the same library type and condition share the same mean.
  expect_lt(spread("RFP_WT_1", "RFP_WT_2"), 0.1)
  expect_lt(spread("RFP_Mutant_1", "RFP_Mutant_2"), 0.1)
  expect_gt(spread("RNA_WT_1", "RNA_WT_2"), 1)
  expect_gt(spread("RNA_Mutant_1", "RNA_Mutant_2"), 1)
})

test_that("the stored true condition effect is the simulated condition effect", {
  set.seed(3)
  constant_dispersion <- function(x) x * 0 + 1e-8
  counts <- quiet_count_tables(
    n = 400, libtypes = "RFP", conditions = c("WT", "Mutant"), replicates = 1,
    interceptMean = 16, interceptSD = 0.2, betaSD = 1,
    betaLibSD = c(RFP = 0), dispMeanRel = constant_dispersion
  )
  truth <- as.data.frame(S4Vectors::mcols(counts))
  observed <- log2(SummarizedExperiment::assay(counts)[, "RFP_Mutant_1"] /
                     SummarizedExperiment::assay(counts)[, "RFP_WT_1"])
  expect_equal(unname(observed), truth$trueBeta, tolerance = 0.01)
  expect_gt(stats::sd(truth$trueBeta), 0.5)

  # With library-specific noise the stored intercept is the first library's coefficient.
  with_noise <- quiet_count_tables(
    n = 400, libtypes = c("RFP", "RNA"), conditions = c("WT", "Mutant"), replicates = 1,
    interceptMean = 16, interceptSD = 0.2, betaSD = 1,
    betaLibSD = c(RFP = 1, RNA = 1), dispMeanRel = constant_dispersion
  )
  noise_truth <- as.data.frame(S4Vectors::mcols(with_noise))
  expect_equal(unname(log2(SummarizedExperiment::assay(with_noise)[, "RFP_WT_1"])),
               noise_truth$trueIntercept, tolerance = 0.01)

  single <- quiet_count_tables(n = 20, libtypes = "RFP", conditions = "WT", replicates = 2)
  expect_true(all(S4Vectors::mcols(single)$trueBeta == 0))
})

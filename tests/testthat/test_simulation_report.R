# A small run through all four steps, in one temporary folder.
run_four_steps <- function(seed = 7, cds_exons = 2) {
  out <- tempfile("coverageSim-report-")
  dirs <- file.path(out, c("genome", "reads", "experiment"))
  for (dir in dirs) dir.create(dir, recursive = TRUE)
  set.seed(seed)
  genome <- suppressWarnings(suppressMessages(simGenome(
    # 12 genes: with fewer, draws taken in a different order can still give
    # the same strands by chance, and the comparison below would miss it.
    n = 12, out_dir = dirs[1], max_uorfs = 1, cds_exons = cds_exons, debug_on = FALSE
  )))
  counts <- suppressMessages(simCountTables(
    ORFik::loadRegion(genome["txdb"], "cds"), libtypes = "RFP",
    plot_PCA = FALSE, print_statistics = FALSE, interceptMean = 6
  ))
  regions <- suppressMessages(simCountTablesRegions(counts))
  suppressWarnings(suppressMessages(simNGScoverage(
    genome, regions[, 1], out_dir = dirs[2], exp_name = basename(out),
    exp_save_dir = dirs[3], validate = FALSE, libFormats = list(RFP = "bam"),
    ground_truth = TRUE
  )))
  list(out = out, genome = genome, counts = regions, rng_after = .Random.seed,
       library = readLines(file.path(dirs[2], "RFP_WT_1_ground_truth.tsv")),
       report = file.path(dirs[2], "simulation_report.tsv"))
}

test_that("the simulation report lists every setting of every step", {
  run <- run_four_steps()
  report <- data.table::fread(run$report, sep = "\t")
  steps <- list(simGenome = simGenome, simCountTables = simCountTables,
                simCountTablesRegions = simCountTablesRegions,
                simNGScoverage = simNGScoverage)
  for (step_name in names(steps)) {
    recorded <- unique(sub("[$[].*$", "", report$setting[report$step == step_name]))
    recorded <- setdiff(recorded, "(random state at start)")
    expect_setequal(recorded, setdiff(names(formals(steps[[step_name]])), "..."))
  }
  # Written by the caller versus filled in by the function.
  expect_false(report[step == "simGenome" & setting == "n", default_used])
  expect_true(report[step == "simGenome" & setting == "genome_name", default_used])
  expect_equal(report[step == "simCountTables" & setting == "plot_PCA", value],
               "not used in this run")
  expect_true(all(c("R", "coverageSim", "coverageSim commit") %in%
                    report[step == "versions", setting]))
  expect_equal(report[setting == "(random state at start)", value],
               rep("stored in simulation_report.rds", 4))
  expect_true(file.exists(file.path(dirname(run$report), "simulation_report.rds")))
  expect_length(readRDS(file.path(dirname(run$report), "simulation_report.rds"))$records, 4L)
})

test_that("recording the settings does not change the simulation", {
  # Some defaults draw random numbers when first used, so a record that read a
  # setting too early would change every result after it. Compare against runs
  # where recording reads nothing at all.
  for (cds_exons in c(2, 1)) {
    with_record <- run_four_steps(cds_exons = cds_exons)
    local({
      local_mocked_bindings(make_step_record = function(step, ...) {
        list(step = step, settings = list(), not_used = character(0),
             default_used = character(0), rng_state_at_start = NULL)
      })
      without_record <<- run_four_steps(cds_exons = cds_exons)
    })
    expect_identical(with_record$rng_after, without_record$rng_after)
    expect_identical(
      unname(tools::md5sum(with_record$genome["gtf"])),
      unname(tools::md5sum(without_record$genome["gtf"]))
    )
    expect_identical(unname(SummarizedExperiment::assays(with_record$counts)),
                     unname(SummarizedExperiment::assays(without_record$counts)))
    expect_identical(with_record$library, without_record$library)
    expect_identical(
      unname(tools::md5sum(with_record$genome["genome"])),
      unname(tools::md5sum(without_record$genome["genome"]))
    )
  }
})

test_that("a genome made outside the run is reported as missing, not invented", {
  run <- run_four_steps()
  record_file <- genome_record_path(run$genome["genome"])
  expect_true(file.exists(record_file))
  unlink(record_file)
  reads <- file.path(run$out, "reads2")
  dir.create(reads)
  suppressWarnings(suppressMessages(simNGScoverage(
    run$genome, run$counts[, 1], out_dir = reads, exp_name = basename(reads),
    exp_save_dir = file.path(run$out, "experiment"), validate = FALSE,
    libFormats = list(RFP = "bam")
  )))
  report <- data.table::fread(file.path(reads, "simulation_report.tsv"), sep = "\t")
  expect_equal(report[step == "simGenome", setting], "(no record)")
})

test_that("settings a run does not use are not evaluated for the record", {
  # If the record evaluated them, each of these calls would stop.
  set.seed(1)
  counts <- suppressMessages(simCountTables(
    n = 5, libtypes = "RFP", conditions = "WT", replicates = 2,
    betaSD = stop("betaSD was evaluated"), plot_PCA = stop("plot_PCA was evaluated"),
    print_statistics = FALSE
  ))
  record <- S4Vectors::metadata(counts)$simulation_record[[1]]
  expect_setequal(record$not_used, c("betaSD", "plot_PCA"))
  genome <- suppressWarnings(suppressMessages(simGenome(
    n = 2, out_dir = tempfile("genome-"), max_uorfs = 0, debug_on = FALSE,
    seqnames = c("chrA", "chrB"), chromosome_count = stop("chromosome_count was evaluated")
  )))
  record <- readRDS(genome_record_path(genome["genome"]))
  expect_equal(record$not_used, "chromosome_count")
  expect_equal(record$settings$seqnames, c("chrA", "chrB"))
})

test_that("the readable report copes with every kind of setting value", {
  values <- list(
    nothing = NULL, empty = character(0), missing = NA, factor = factor(c("a", "b")),
    many = seq_len(50), same = rep(120, 3000), words = paste0("gene", 1:20),
    fun = function(x) 4 / x + 0.1, expression = quote(rep(c(1, 0, 0), length.out = x)),
    table = data.table::data.table(variable = "median", seqs = "AAA", alpha = 1),
    list_column = data.table::data.table(id = 1:2, values = list(1, 2:3)),
    ranges = GenomicRanges::GRanges("chr1", IRanges::IRanges(1:3, width = 5)),
    nested = list(cds = list(RFP = 0.75, RNA = 0.4), unnamed = list(1, 2)),
    odd_names = stats::setNames(list(1, 2), c(NA, ""))
  )
  rows <- data.table::rbindlist(Map(setting_rows, names(values), values))
  expect_true(all(nzchar(rows$value)))
  expect_equal(rows[setting == "same", value], "3000 values, all 120")
  expect_equal(rows[setting == "nested$cds$RFP", value], "0.75")
  expect_true(all(c("odd_names[[1]]", "odd_names[[2]]") %in% rows$setting))
})

test_that("the genome record never takes the place of the genome itself", {
  for (fasta in c("g.fasta", "g.fa", "g.fna", "g.fasta.gz", "genome")) {
    expect_false(genome_record_path(fasta) == fasta)
    expect_match(genome_record_path(fasta), "_simulation_record\\.rds$")
  }
})

# Running the samples of simNGScoverage(), one after another or side by side.

# `cores` as a whole number. Windows cannot fork processes, which
# parallel::mclapply() needs, so there only one core is used.
resolve_cores <- function(cores) {
  if (!is.numeric(cores) || length(cores) != 1L || is.na(cores) ||
      cores < 1 || cores != round(cores)) {
    stop("cores must be a single whole number of at least 1")
  }
  cores <- as.integer(cores)
  if (cores > 1L && .Platform$OS.type == "windows") {
    warning("cores > 1 is not supported on Windows; using 1 core")
    cores <- 1L
  }
  cores
}

# Evaluate `code` with the random generator seeded by `seed`, then put the
# random state back to what it was before. A process started by
# parallel::mclapply() begins without a random state, so there may be none to
# put back.
with_seed <- function(seed, code) {
  old_state <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(if (is.null(old_state)) {
    rm(".Random.seed", envir = globalenv())
  } else {
    assign(".Random.seed", old_state, envir = globalenv())
  })
  set.seed(seed)
  code
}

# Call simulate_sample(s) for every sample s, with seeds[s] as its seed, on up to
# `cores` processes, and return the results in sample order.
#
# Each sample puts the random state back afterwards, so the state after this
# call is the same for any number of cores. Warnings raised inside a separate
# process would be lost, so they are collected for every sample and raised
# again here, after all samples have finished and before an error of a sample
# is passed on.
run_samples <- function(simulate_sample, seeds, cores = 1L) {
  run_one <- function(s) {
    warnings <- character()
    value <- tryCatch(
      withCallingHandlers(
        with_seed(seeds[s], simulate_sample(s)),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) e
    )
    list(value = value, warnings = warnings)
  }
  samples <- seq_along(seeds)
  results <- if (cores > 1L) {
    parallel::mclapply(samples, run_one, mc.cores = cores,
                       mc.preschedule = FALSE)
  } else {
    lapply(samples, run_one)
  }
  # A process that was killed returns NULL or a "try-error" text instead.
  finished <- vapply(results, is.list, logical(1))
  for (warning_text in unlist(lapply(results[finished], `[[`, "warnings"))) {
    warning(warning_text, call. = FALSE)
  }
  for (s in samples) {
    if (!finished[s]) {
      stop("Sample ", s, " did not finish: its process stopped without a ",
           "result, for example because memory ran out. Try fewer cores.")
    }
    if (inherits(results[[s]]$value, "error")) stop(results[[s]]$value)
  }
  lapply(results, `[[`, "value")
}

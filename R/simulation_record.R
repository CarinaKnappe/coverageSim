# Simulation record: what each step ran with
#
# Every step function (simGenome(), simCountTables(), simCountTablesRegions(),
# simNGScoverage()) records the settings it used -- the ones the caller wrote
# down and the defaults it filled in itself -- together with the random number
# state it started from. simNGScoverage() collects the records of all four
# steps and writes them next to the simulated libraries, so a stored run says
# by itself how it was made.
#
# Recording must not change the simulation. Some defaults draw random numbers
# themselves (simGenome()'s strand, for instance) and R only evaluates a
# setting when it is first used, so reading one early would change the order
# of the draws and with it every result. Each step therefore takes its record
# at its very end, when it has used every setting in its own order, and leaves
# alone the settings it did not use in this run. Settings a step rewrites on
# the way are kept as they were given.
#
# When adding to the code:
# - A new argument of a step function, or a new entry in a list setting such
#   as fragment_geometry, appears in the report by itself (a test checks that
#   every argument is listed).
# - An argument that a step only uses in some runs, and whose default draws
#   random numbers or stops with an error, must be added to that step's
#   `not_used` list for those runs; otherwise reading it at the end would draw
#   or stop. Fixed defaults (5, TRUE, "MN") need nothing.
# - A new step function needs the same two lines at its start (rng_state,
#   supplied) and make_step_record() at its end.
# - Values computed outside the step functions, in a workflow script for
#   instance, only appear if they are passed to one of them.

# Packages whose versions can change a simulated result for the same seed.
simulation_record_packages <- c(
  "coverageSim", "ORFik", "Biostrings", "GenomicRanges", "IRanges",
  "data.table", "extraDistr", "DESeq2", "txdbmaker", "rtracklayer"
)

# The random number state right now, read without drawing from it; NULL when
# nothing has set or used the generator in this R session yet.
current_rng_state <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
}

# One step's record. `fun` is the step function, `env` its own environment, `supplied`
# the names of the arguments the caller wrote down (from match.call(), taken
# before any of them is used), `rng_state` the state from the step's first
# line. Settings listed in `not_used` are not read, because reading them would
# evaluate them. Settings in `summarize_only` are large objects (a count table,
# for instance) that are described, not stored.
make_step_record <- function(step, fun, env, supplied, rng_state,
                             not_used = character(0),
                             summarize_only = character(0)) {
  setting_names <- setdiff(names(formals(fun)), "...")
  used <- setdiff(setting_names, not_used)
  values <- mget(used, envir = env)
  values[intersect(summarize_only, used)] <- lapply(
    values[intersect(summarize_only, used)], describe_object
  )
  values <- functions_as_text(values)
  list(
    step = step,
    settings = values,
    not_used = intersect(not_used, setting_names),
    default_used = setdiff(setting_names, supplied),
    rng_state_at_start = rng_state,
    rng_kind = RNGkind()
  )
}

# A function carries the environment it was made in -- for a default argument,
# the whole running step with all its intermediate tables -- so keeping it
# would store all of that. Its source text says what it does.
functions_as_text <- function(x) {
  if (is.function(x)) return(describe_setting_value(x))
  if (is.list(x) && !is.data.frame(x)) x[] <- lapply(x, functions_as_text)
  x
}

# A short text for an object too large or too structured to store as it is.
describe_object <- function(x) {
  if (methods::is(x, "SummarizedExperiment")) {
    assay_names <- assayNames(x)
    return(paste0(class(x)[1], ": ", nrow(x), " genes x ", ncol(x),
                  " libraries, ", length(assays(x)), " assay(s)",
                  if (length(assay_names)) paste0(": ", paste(assay_names, collapse = ", "))))
  }
  if (methods::is(x, "GRanges") || methods::is(x, "GRangesList")) {
    return(paste0(class(x)[1], " with ", length(x), " elements"))
  }
  paste0(class(x)[1], " of length ", length(x))
}

# Versions of R and of the packages that compute, plus the code state of
# coverageSim itself (commit and whether anything was unsaved).
simulation_versions <- function() {
  package_versions <- vapply(simulation_record_packages, function(package) {
    version <- tryCatch(as.character(utils::packageVersion(package)),
                        error = function(e) NA_character_)
    version
  }, character(1))
  package_dir <- tryCatch(find.package("coverageSim"), error = function(e) ".")
  list(
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    package_versions = package_versions,
    code_version = current_code_version(package_dir)
  )
}

# Attach a step record to a count table, after the records it already carries.
add_step_record <- function(count_table, record) {
  records <- S4Vectors::metadata(count_table)$simulation_record
  S4Vectors::metadata(count_table)$simulation_record <- c(records, list(record))
  count_table
}

# Where simGenome() keeps its record: next to the genome, named after it, so
# two genomes in one folder keep separate records.
genome_record_path <- function(genome_fasta) {
  paste0(tools::file_path_sans_ext(genome_fasta, compression = TRUE),
         "_simulation_record.rds")
}

# Flatten one setting into rows of setting name and readable value. Nested
# lists become one row per entry (region_proportion$cds$RFP), tables and
# functions a short description.
setting_rows <- function(name, value) {
  if (is.list(value) && !is.data.frame(value) && length(value) > 0L) {
    entry_names <- names(value)
    if (is.null(entry_names)) entry_names <- rep("", length(value))
    unnamed <- is.na(entry_names) | entry_names == ""
    entry_names[unnamed] <- paste0("[[", which(unnamed), "]]")
    rows <- Map(function(entry, entry_name) {
      separator <- if (startsWith(entry_name, "[[")) "" else "$"
      setting_rows(paste0(name, separator, entry_name), entry)
    }, value, entry_names)
    return(data.table::rbindlist(rows))
  }
  data.table::data.table(setting = name, value = describe_setting_value(value))
}

describe_setting_value <- function(x) {
  if (is.null(x)) return("NULL")
  if (is.function(x) || is.language(x)) {
    text <- paste(trimws(deparse(x)), collapse = " ")
    return(if (nchar(text) > 200L) paste0(substr(text, 1L, 197L), "...") else text)
  }
  if (is.data.frame(x)) {
    # Codon-bias tables name their source in "variable" (R1 to R10, median, ...).
    source <- if ("variable" %in% names(x)) {
      paste0("; variable: ", paste(unique(as.character(x$variable)), collapse = ", "))
    }
    return(paste0("table with ", nrow(x), " rows; columns: ",
                  paste(names(x), collapse = ", "), source))
  }
  if (is.list(x)) return("empty list")
  if (!is.atomic(x)) return(describe_object(x))
  if (length(x) == 0L) return("empty")
  if (length(x) <= 10L) {
    text <- as.character(x)
    if (!is.null(names(x))) text <- paste0(names(x), " = ", text)
    return(paste(text, collapse = ", "))
  }
  if (length(unique(x)) == 1L) {
    return(paste0(length(x), " values, all ", as.character(x[1])))
  }
  if (is.numeric(x)) {
    return(paste0(length(x), " values, min ", min(x), ", median ",
                  stats::median(x), ", max ", max(x)))
  }
  paste0(length(x), " values: ", paste(utils::head(as.character(x), 3L),
                                       collapse = ", "), ", ...")
}

# The readable table of one record: one row per setting.
record_table <- function(record) {
  rows <- data.table::rbindlist(c(
    list(data.table::data.table(setting = character(0), value = character(0))),
    Map(setting_rows, names(record$settings), record$settings)
  ))
  top_level <- sub("[$[].*$", "", rows$setting)
  rows[, default_used := top_level %in% record$default_used]
  if (length(record$not_used)) {
    rows <- rbind(rows, data.table::data.table(
      setting = record$not_used, value = "not used in this run",
      default_used = record$not_used %in% record$default_used
    ))
  }
  rows <- rbind(rows, data.table::data.table(
    setting = "(random state at start)",
    value = if (is.null(record$rng_state_at_start)) {
      "none: no seed was set and nothing had drawn random numbers yet"
    } else {
      "stored in simulation_report.rds"
    },
    default_used = NA
  ))
  rows[, step := record$step]
  data.table::setcolorder(rows, c("step", "setting", "value", "default_used"))
  rows[]
}

# Written by simNGScoverage(): the records of all steps it can find, plus the
# versions, as simulation_report.rds (complete) and simulation_report.tsv
# (readable).
write_simulation_report <- function(out_dir, genome_fasta, count_table, step4_record) {
  genome_record_file <- genome_record_path(genome_fasta)
  genome_record <- if (file.exists(genome_record_file)) {
    readRDS(genome_record_file)
  } else {
    NULL
  }
  records <- c(list(genome_record),
               S4Vectors::metadata(count_table)$simulation_record,
               list(step4_record))
  records <- Filter(Negate(is.null), records)
  versions <- simulation_versions()
  report <- list(records = records, versions = versions,
                 recorded_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))

  recorded_steps <- vapply(records, `[[`, character(1), "step")
  missing_steps <- setdiff(c("simGenome", "simCountTables", "simCountTablesRegions"),
                           recorded_steps)
  version_rows <- data.table::data.table(
    step = "versions",
    setting = c("R", names(versions$package_versions), "coverageSim commit",
                "coverageSim unsaved changes"),
    value = c(versions$r_version, versions$package_versions,
              versions$code_version$commit,
              if (isTRUE(versions$code_version$unsaved_changes)) {
                paste("yes:", versions$code_version$unsaved_files)
              } else if (isFALSE(versions$code_version$unsaved_changes)) {
                "no"
              } else {
                "unknown (not run from a git checkout)"
              }),
    default_used = NA
  )
  missing_rows <- if (length(missing_steps)) {
    data.table::data.table(
      step = missing_steps, setting = "(no record)",
      value = "this step was not run through coverageSim in this run, or its record was not kept",
      default_used = NA
    )
  }
  table <- data.table::rbindlist(c(lapply(records, record_table),
                                   list(missing_rows, version_rows)),
                                 use.names = TRUE)
  # The libraries are already written at this point; a report that cannot be
  # saved should not throw them away, but it must not go unnoticed either.
  written <- tryCatch({
    data.table::fwrite(table, file.path(out_dir, "simulation_report.tsv"), sep = "\t")
    saveRDS(report, file.path(out_dir, "simulation_report.rds"))
    TRUE
  }, error = function(e) {
    warning("The simulation finished, but its report could not be written to ",
            out_dir, ": ", conditionMessage(e), call. = FALSE)
    FALSE
  })
  invisible(if (written) report)
}

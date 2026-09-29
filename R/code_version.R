#' Record which version of the code produced a simulation
#'
#' Captures enough about the currently loaded coverageSim source to answer, months
#' later, the question "which code produced this figure?". Write the result next to
#' every simulation output that anyone might come back to.
#'
#' Recording the commit alone is not sufficient here, and that is the reason this
#' function exists rather than a one-line call to `git rev-parse`. The analysis
#' workflows load the package with `devtools::load_all()`, which loads the files as
#' they currently sit in the working directory -- including edits that were never
#' saved into the project history. A commit identifier taken at that moment
#' describes the last saved state, not necessarily the state that ran. So the
#' record also states whether anything was unsaved, and which files, because that
#' is exactly the difference between a result somebody can reproduce and one
#' nobody can.
#'
#' @param repo_dir path to the coverageSim source directory whose version should be
#'   recorded. Defaults to the working directory, which is where
#'   `devtools::load_all()` loads from in the analysis workflows.
#' @return A one-row `data.table` with these columns:
#'   \describe{
#'     \item{commit}{full identifier of the last saved state, or `NA` when
#'       `repo_dir` is not under version control or git is unavailable.}
#'     \item{commit_short}{the same identifier abbreviated, for printing.}
#'     \item{commit_date}{when that state was saved (ISO date).}
#'     \item{branch}{the line of development it was saved on.}
#'     \item{unsaved_changes}{`TRUE` if source files differed from the last saved
#'       state when this was called -- meaning the run used code that is not in the
#'       project history.}
#'     \item{unsaved_files}{those files, comma-separated; empty string if none.}
#'     \item{fully_reproducible}{`TRUE` only when the commit is known *and* nothing
#'       was unsaved. This is the single column to look at when judging whether a
#'       stored result can be recreated.}
#'     \item{recorded_at}{when this record was written.}
#'     \item{package_version, r_version}{versions of coverageSim and R.}
#'   }
#' @examples
#' \dontrun{
#' version <- current_code_version()
#' data.table::fwrite(version, file.path(out_dir, "code_version.tsv"), sep = "\t")
#' }
#' @export
current_code_version <- function(repo_dir = ".") {
  git <- function(...) {
    output <- suppressWarnings(try(
      system2("git", c("-C", shQuote(normalizePath(repo_dir, mustWork = FALSE)), ...),
              stdout = TRUE, stderr = FALSE),
      silent = TRUE
    ))
    status <- attr(output, "status")
    if (inherits(output, "try-error") || (!is.null(status) && status != 0)) {
      return(character(0))
    }
    output
  }

  commit <- git("rev-parse", "HEAD")
  if (!length(commit)) {
    # Not under version control, or git is not installed. Say so plainly instead
    # of inventing an identifier: a missing value is honest, a guessed one is not.
    return(data.table::data.table(
      commit = NA_character_, commit_short = NA_character_,
      commit_date = NA_character_, branch = NA_character_,
      unsaved_changes = NA, unsaved_files = "",
      fully_reproducible = FALSE,
      recorded_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      package_version = as.character(utils::packageVersion("coverageSim")),
      r_version = paste(R.version$major, R.version$minor, sep = ".")
    ))
  }

  # --porcelain lists one line per file that differs from the last saved state.
  # Untracked files are deliberately included: a workflow can just as easily have
  # run against a source file that was never added to the project at all.
  changed <- git("status", "--porcelain")
  changed_files <- if (length(changed)) trimws(substring(changed, 4L)) else character(0)

  data.table::data.table(
    commit = commit[[1]],
    commit_short = substr(commit[[1]], 1L, 7L),
    commit_date = { d <- git("log", "-1", "--format=%ad", "--date=short")
                    if (length(d)) d[[1]] else NA_character_ },
    branch = { b <- git("rev-parse", "--abbrev-ref", "HEAD")
               if (length(b)) b[[1]] else NA_character_ },
    unsaved_changes = length(changed_files) > 0L,
    unsaved_files = paste(changed_files, collapse = ", "),
    fully_reproducible = length(changed_files) == 0L,
    recorded_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    package_version = as.character(utils::packageVersion("coverageSim")),
    r_version = paste(R.version$major, R.version$minor, sep = ".")
  )
}

#' Write the code version next to a simulation output
#'
#' Convenience wrapper around [current_code_version()] for the analysis workflows:
#' it writes `code_version.tsv` into an output directory and returns the record
#' invisibly, so a workflow can record its provenance in one line.
#'
#' @param out_dir directory to write `code_version.tsv` into. Created if missing.
#' @param repo_dir passed to [current_code_version()].
#' @return The record, invisibly.
#' @examples
#' \dontrun{
#' write_code_version(out_reads_dir)
#' }
#' @export
write_code_version <- function(out_dir, repo_dir = ".") {
  version <- current_code_version(repo_dir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(version, file.path(out_dir, "code_version.tsv"), sep = "\t")
  if (isTRUE(version$unsaved_changes)) {
    # Worth a message rather than silence: the person starting a long simulation
    # can still stop and save their work, which is the difference between a
    # reproducible result and one that only exists as numbers.
    message(
      "Note: ", length(strsplit(version$unsaved_files, ", ")[[1]]),
      " file(s) differ from the last saved state, so this run cannot be exactly ",
      "reproduced from the project history alone. Recorded in code_version.tsv."
    )
  }
  invisible(version)
}

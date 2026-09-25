# data.table::frollapply()'s first two parameters have been named
# differently across data.table versions (X/N vs. the now-deprecated x/n).
# Rather than blindly retrying all four X/x x N/n combinations and
# swallowing every error along the way (a real error thrown by FUN's own
# body was previously masked -- all four attempts fail identically, and
# only the last, argument-matching-sounding error message ever surfaced,
# not FUN's actual one), inspect the installed frollapply()'s own formals
# once to find its current parameter names, then make exactly one call
# with them. X/N is preferred when both exist (current data.table warns
# that x/n is deprecated); x/n is only used as a fallback for older
# data.table versions that never had X/N.
frollapply_compat <- function(x, window, FUN, ..., align = "center", fill = NA) {
  params <- names(formals(data.table::frollapply))
  x_name <- if ("X" %in% params) "X" else "x"
  n_name <- if ("N" %in% params) "N" else "n"
  args <- c(
    stats::setNames(list(x, window), c(x_name, n_name)),
    list(FUN = FUN, ..., align = align, fill = fill)
  )
  do.call(data.table::frollapply, args)
}

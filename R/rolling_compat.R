frollapply_compat <- function(x, window, FUN, ..., align = "center", fill = NA) {
  common_args <- list(FUN = FUN, ..., align = align, fill = fill)
  attempts <- list(
    c(list(X = x, N = window), common_args),
    c(list(x = x, N = window), common_args),
    c(list(X = x, n = window), common_args),
    c(list(x = x, n = window), common_args)
  )

  last_error <- NULL
  for (args in attempts) {
    result <- try(
      do.call(data.table::frollapply, args),
      silent = TRUE
    )
    if (!inherits(result, "try-error")) {
      return(result)
    }
    last_error <- attr(result, "condition")
  }

  stop(last_error)
}

# Weighted sum over a centred window, for every position of `x`: position i
# gets weights[1] * x[i - r] + ... + weights[length(weights)] * x[i + r], with
# r = (length(weights) - 1) / 2. Positions whose window runs past either end of
# `x`, or whose window contains an NA, get NA.
#
# stats::filter() computes this in compiled code for the whole vector at once.
# The previous version used data.table::frollapply() with an R function that
# was called once per position, which took about a third of the time of placing
# reads (profile of 2026-11-04).
rolling_weighted_sum <- function(x, weights) {
  if (length(weights) %% 2L != 1L) {
    stop("rolling_weighted_sum() needs an odd number of weights")
  }
  if (length(x) < length(weights)) return(rep(NA_real_, length(x)))
  # filter() multiplies the window from its right end backwards, hence rev().
  as.numeric(stats::filter(as.numeric(x), rev(weights),
                           method = "convolution", sides = 2L))
}

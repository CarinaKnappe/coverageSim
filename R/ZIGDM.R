#' Draw from a (generalized) Dirichlet-multinomial, optionally zero-inflated
#'
#' Extends \code{extraDistr::rdirmnom()} in two independent, composable ways,
#' matching the Zero-Inflated Generalized Dirichlet-Multinomial ("ZIGDM")
#' idea named as future work in the coverageSim manuscript (the standard
#' Dirichlet-multinomial used in Step 4 slightly underestimates skew/the 3rd
#' statistical moment relative to real Ribo-seq data):
#' \itemize{
#'   \item \strong{Zero-inflation}: for every draw, each non-padding position
#'   is independently excluded (forced to a structural zero, i.e. its alpha
#'   is set to \code{pad_value}) with probability \code{zero_inflation},
#'   before sampling. This adds a genuine two-component mixture (structural
#'   zero vs. Dirichlet-multinomial-drawn value) on top of the smooth
#'   compound-multinomial distribution, which is what increases higher
#'   moments (skew) relative to a plain Dirichlet-multinomial draw.
#'   \item \strong{Generalized Dirichlet} (Connor-Mosimann stick-breaking):
#'   instead of a single Dirichlet(alpha) draw per row, sequentially draws
#'   break fractions \code{V_k ~ Beta(alpha_k, beta_k)}, where
#'   \code{beta_k = gdm_scale * sum(alpha[(k+1):K])}. \code{gdm_scale = 1}
#'   reproduces the standard Dirichlet exactly (a known identity); values
#'   away from 1 give a more flexible position-to-position covariance
#'   structure than a single Dirichlet allows. The same stick-breaking chain
#'   is combined with sequential binomial thinning to draw multinomial
#'   counts directly (the standard construction of a multinomial draw from a
#'   sequence of conditional binomials), which is at the same time the
#'   sampling algorithm for the generalized Dirichlet-multinomial
#'   distribution.
#' }
#' With both extensions at their defaults (\code{family = "dirichlet"},
#' \code{zero_inflation = 0}), this calls \code{extraDistr::rdirmnom()}
#' directly and is therefore exactly backward compatible with all existing
#' callers.
#'
#' Note on \code{gdm_scale}: unlike \code{dmn_alpha_scale} (which only
#' changes variance, never the expected coverage), \code{gdm_scale != 1}
#' also shifts the expected per-position proportions slightly, because it
#' changes how the "remaining probability mass" is split at each step of the
#' stick-breaking chain. At \code{gdm_scale = 1} this does not apply (exact
#' reduction to the standard Dirichlet). Treat \code{gdm_scale} as a reshape
#' control, not a pure variance knob, when moving it away from 1.
#' @param n integer, number of independent draws (rows of \code{alpha}).
#' @param size integer vector of length \code{n}, total count per draw.
#' @param alpha numeric matrix (\code{n} x \code{K}), Dirichlet/generalized-
#' Dirichlet shape parameters per position, all finite and strictly positive
#' (use \code{pad_value}, not 0, for absent positions). Positions with
#' \code{alpha <= pad_value} are treated as structurally absent (matches the
#' ragged-matrix padding convention used by \code{pack_alpha_rows()}) and are
#' never eligible for zero-inflation exclusion, since they are excluded
#' already.
#' @param family character, \code{"dirichlet"} (default) or
#' \code{"generalized"}.
#' @param gdm_scale positive numeric, default 1. Only affects the draw when
#' \code{family == "generalized"}, but is validated in every case. See Details
#' above.
#' @param zero_inflation numeric in [0, 1), default 0 (off). Probability
#' that an eligible position (\code{alpha > pad_value}) is excluded for a
#' given draw. If every eligible position in a row would be excluded, the
#' single highest-alpha position is kept instead, so a region with
#' \code{size > 0} always has somewhere to sample its reads (mirrors the
#' \code{boundary_action = "renormalize"} convention used elsewhere in this
#' package, e.g. in \code{fragment_geometry.R}).
#' @param pad_value numeric, default 1e-24, matches \code{pack_alpha_rows()}.
#' Masked and padded positions keep \code{pad_value} as their alpha, so it
#' must be negligibly small compared to real alpha values: the
#' \code{"dirichlet"} family relies on it to give them (practically) no
#' probability mass, whereas the \code{"generalized"} family forces exact
#' zeros there.
#' @return integer (numeric) matrix (\code{n} x \code{K}) of sampled counts;
#' every row sums to the corresponding \code{size}.
#' @export
#' @examples
#' alpha <- matrix(c(5, 3, 1, 1e-24, 2, 2, 2, 2), nrow = 2, byrow = TRUE)
#' # Exactly extraDistr::rdirmnom() -- default, fully backward compatible
#' rgdirmnom(2, size = c(10, 20), alpha = alpha)
#' # With zero-inflation and a generalized Dirichlet
#' rgdirmnom(2, size = c(10, 20), alpha = alpha,
#'           family = "generalized", gdm_scale = 1.5, zero_inflation = 0.3)
rgdirmnom <- function(n, size, alpha,
                      family = c("dirichlet", "generalized"),
                      gdm_scale = 1, zero_inflation = 0,
                      pad_value = 1e-24) {
  family <- match.arg(family)
  stopifnot(is.matrix(alpha), is.numeric(alpha), nrow(alpha) == n,
           length(size) == n, is.numeric(size), !anyNA(size),
           all(is.finite(size)), all(size >= 0), all(size == round(size)),
           all(is.finite(alpha)), all(alpha > 0),
           is.numeric(gdm_scale), length(gdm_scale) == 1L,
           is.finite(gdm_scale), gdm_scale > 0,
           is.numeric(zero_inflation), length(zero_inflation) == 1L,
           is.finite(zero_inflation), zero_inflation >= 0, zero_inflation < 1,
           is.numeric(pad_value), length(pad_value) == 1L,
           is.finite(pad_value), pad_value > 0)
  if (any(size > 0 & rowSums(alpha > pad_value) == 0L)) {
    stop("every row with size > 0 needs at least one position with alpha > pad_value",
         call. = FALSE)
  }

  if (family == "dirichlet" && zero_inflation == 0) {
    return(extraDistr::rdirmnom(n = n, size = size, alpha = alpha))
  }

  if (zero_inflation > 0) {
    alpha <- apply_zero_inflation_mask(alpha, zero_inflation, pad_value)
  }

  if (family == "dirichlet") {
    return(extraDistr::rdirmnom(n = n, size = size, alpha = alpha))
  }
  rgdirmnom_generalized(n, size, alpha, gdm_scale, pad_value)
}

#' Randomly mask eligible positions to structural zeros (zero-inflation)
#' @param alpha numeric matrix (n x K) of Dirichlet/GD alpha values.
#' @param zero_inflation numeric in [0, 1), exclusion probability.
#' @param pad_value numeric, positions at or below this are already
#' structurally absent and are left untouched.
#' @return alpha, with a random subset of eligible entries per row set to
#' pad_value; at least one eligible entry per row is always kept.
#' @keywords internal
apply_zero_inflation_mask <- function(alpha, zero_inflation, pad_value) {
  eligible <- alpha > pad_value
  exclude <- matrix(stats::runif(length(alpha)) < zero_inflation,
                    nrow = nrow(alpha)) & eligible

  survives <- eligible & !exclude
  all_excluded <- rowSums(survives) == 0L & rowSums(eligible) > 0L
  if (any(all_excluded)) {
    masked_alpha <- alpha
    masked_alpha[!eligible] <- -Inf
    keep <- max.col(masked_alpha[all_excluded, , drop = FALSE], ties.method = "first")
    exclude[cbind(which(all_excluded), keep)] <- FALSE
  }

  alpha[exclude] <- pad_value
  alpha
}

#' Draw a generalized Dirichlet-multinomial sample via stick-breaking
#'
#' Sequential Beta-Binomial chain: draws the K-1 Connor-Mosimann break
#' fractions and thins the remaining count at each step. At
#' \code{gdm_scale = 1} this is mathematically the standard
#' Dirichlet-multinomial (see \code{rgdirmnom()} Details).
#' @param n,size,alpha,pad_value see \code{rgdirmnom()}.
#' @param gdm_scale see \code{rgdirmnom()}.
#' @return integer (numeric) matrix (n x K), rows sum to \code{size}.
#' @keywords internal
rgdirmnom_generalized <- function(n, size, alpha, gdm_scale, pad_value) {
  K <- ncol(alpha)
  counts <- matrix(0, nrow = n, ncol = K)
  remaining_size <- as.numeric(size)
  remaining_alpha_sum <- rowSums(alpha)
  # The last live (non-padded, non-masked) column of each row takes all
  # remaining reads by construction, so they never depend on the numerical
  # behaviour of rbeta() with a vanishing second shape parameter and can
  # never spill into a structurally absent trailing column.
  last_live <- last_live_column(alpha, pad_value)

  for (k in seq_len(K - 1L)) {
    a_k <- alpha[, k]
    b_k <- gdm_scale * (remaining_alpha_sum - a_k)
    # A row is "degenerate" at this column when the column itself is
    # structurally absent for that row -- either genuine tail padding (a
    # shorter gene, pack_alpha_rows()) or a position zero-inflation just
    # masked out (apply_zero_inflation_mask() also writes pad_value). Either
    # way no reads may land there, so skip the Beta draw and leave
    # remaining_size untouched.
    degenerate <- a_k <= pad_value | remaining_size <= 0

    v_k <- rep.int(0, n)
    live <- !degenerate
    if (any(live)) {
      # b_k is mathematically always > 0 here (there is always at least one
      # more column's worth of positive alpha or pad_value left), but floor
      # it defensively against floating-point cancellation.
      v_k[live] <- stats::rbeta(sum(live), a_k[live], pmax(b_k[live], .Machine$double.eps))
    }
    v_k[k == last_live] <- 1

    draw_k <- rep.int(0, n)
    can_draw <- remaining_size > 0 & !degenerate
    if (any(can_draw)) {
      draw_k[can_draw] <- stats::rbinom(sum(can_draw), remaining_size[can_draw], v_k[can_draw])
    }

    counts[, k] <- draw_k
    remaining_size <- remaining_size - draw_k
    remaining_alpha_sum <- remaining_alpha_sum - a_k
  }
  counts[, K] <- remaining_size
  counts
}

# Index of the last column per row whose alpha is above pad_value (i.e. not
# structurally absent); 0 if the row has none.
last_live_column <- function(alpha, pad_value) {
  live <- alpha > pad_value
  ifelse(rowSums(live) > 0L, ncol(alpha) - max.col(live[, ncol(alpha):1, drop = FALSE], ties.method = "first") + 1L, 0L)
}

# Validators for the simNGScoverage() dmn_gdm_scale/dmn_zero_inflation
# arguments, mirroring validate_dmn_alpha_scale() in
# Sim_reads_from_counts_helpers.R. rgdirmnom() itself re-checks these values
# too (defense in depth); these run earlier, before the simulation pipeline
# starts, for a fast, clear failure.
validate_dmn_gdm_scale <- function(scale) {
  if (!is.numeric(scale) || length(scale) != 1L || is.na(scale) ||
      !is.finite(scale) || scale <= 0) {
    stop("dmn_gdm_scale must be one finite number greater than zero", call. = FALSE)
  }
  invisible(scale)
}

validate_dmn_zero_inflation <- function(zero_inflation) {
  if (!is.numeric(zero_inflation) || length(zero_inflation) != 1L || is.na(zero_inflation) ||
      !is.finite(zero_inflation) || zero_inflation < 0 || zero_inflation >= 1) {
    stop("dmn_zero_inflation must be one finite number in [0, 1)", call. = FALSE)
  }
  invisible(zero_inflation)
}

#' Zero-inflation masking for a single alpha/weight vector
#'
#' Vector counterpart of \code{apply_zero_inflation_mask()}, used by
#' \code{draw_site_probabilities()} (the deferred, simulated-RPF
#' fragment-geometry sampling path), which works with one per-gene weight
#' vector at a time rather than a padded n x K matrix. Delegates to the
#' matrix version via a one-row matrix so both call sites share the same
#' exclusion/guard logic.
#' @param alpha numeric vector of Dirichlet/GD weights for a single gene.
#' @param zero_inflation numeric in [0, 1), exclusion probability.
#' @param pad_value numeric, positions at or below this are already
#' structurally absent (e.g. sim_sequence_bias() maps zero to 1e-24) and are
#' left untouched, matching apply_zero_inflation_mask()'s convention.
#' @return alpha, with a random subset of eligible entries set to pad_value;
#' at least one eligible entry is always kept.
#' @keywords internal
apply_zero_inflation_mask_vector <- function(alpha, zero_inflation, pad_value = 1e-24) {
  as.vector(apply_zero_inflation_mask(matrix(alpha, nrow = 1), zero_inflation, pad_value))
}

#' Draw one generalized-Dirichlet probability vector
#'
#' Same Connor-Mosimann stick-breaking chain as \code{rgdirmnom_generalized()}
#' (see \code{rgdirmnom()} Details), but returns the drawn probability vector
#' directly instead of thinning a read-count budget through it. Used by
#' \code{draw_site_probabilities()}, which needs latent site probabilities
#' (later multiplied by a region's read count) rather than integer counts.
#' At \code{gdm_scale = 1} this is exactly the standard Dirichlet(alpha)
#' distribution (same identity as \code{rgdirmnom_generalized()}); values
#' away from 1 carry the same expected-value caveat described there.
#' @param alpha numeric vector of positive shape values (already
#' zero-inflation-masked if desired; masked/padded entries at or below
#' \code{pad_value} are treated as structurally absent, same convention as
#' \code{rgdirmnom_generalized()}).
#' @param gdm_scale positive numeric.
#' @param pad_value numeric, default 1e-24.
#' @return a numeric probability vector, same length as \code{alpha},
#' summing to 1.
#' @keywords internal
rgdirichlet_generalized <- function(alpha, gdm_scale, pad_value = 1e-24) {
  K <- length(alpha)
  # The K - 1 break fractions are independent Beta draws, so the whole chain
  # is vectorised: no sequential dependence except the final cumulative product.
  # tail_sum[k] = sum(alpha[(k + 1):K]), computed without subtractive cancellation.
  tail_sum <- rev(cumsum(rev(alpha)))[-1L]
  head_alpha <- alpha[-K]
  live <- head_alpha > pad_value
  last_live <- last_live_column(matrix(alpha, nrow = 1), pad_value)
  if (last_live == 0L) {
    stop("alpha needs at least one position above pad_value", call. = FALSE)
  }

  v <- numeric(K - 1L)
  if (any(live)) {
    v[live] <- stats::rbeta(sum(live), head_alpha[live],
                            pmax(gdm_scale * tail_sum[live], .Machine$double.eps))
  }
  if (last_live < K) v[last_live] <- 1

  survival <- cumprod(c(1, 1 - v))
  c(v * survival[-K], survival[K])
}

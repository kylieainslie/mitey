#' Calculate Sample Weighted Variance
#'
#' Computes the sample weighted variance of a numeric vector using precision weights.
#' This function implements the standard unbiased weighted variance estimator commonly
#' used in statistical applications where observations have different precisions
#' or levels of confidence.
#'
#' The weighted variance is calculated using the formula:
#' \deqn{s^2_w = \frac{\sum w_i}{\sum w_i^2 - (\sum w_i)^2} \sum w_i (x_i - \bar{x}_w)^2}
#'
#' where \eqn{\bar{x}_w} is the weighted mean and \eqn{w_i} are the weights.
#'
#' @param x numeric vector; the data values for which to calculate weighted variance.
#'          Missing values are allowed if \code{na.rm = TRUE}
#' @param w numeric vector; the precision weights corresponding to each observation in \code{x}.
#'          These represent how much confidence to place in each measurement (higher = more trusted).
#'          Must be the same length as \code{x}. Should be non-negative
#' @param na.rm logical; if \code{TRUE}, missing values (both \code{NA} and \code{NaN})
#'              are removed before computation. If \code{FALSE} (default), missing
#'              values will cause the function to return \code{NA}
#'
#' @return numeric; the weighted sample variance. Returns \code{NA} if insufficient
#'         data or if \code{na.rm = FALSE} and missing values are present
#'
#' @details
#' This function uses precision weights (also called reliability weights), which represent how much confidence or trust to place in each observation, rather than frequency weights that represent how many times to count each observation.
#'
#' The denominator correction (\eqn{\sum w_i^2 - (\sum w_i)^2}) provides an
#' unbiased estimator for precision weights. Examples of precision weights include
#' probabilities (0-1), measurement confidence scores, or inverse error variances.
#'
#' @seealso \code{\link{weighted.mean}} for weighted mean calculation,
#'          \code{\link{var}} for unweighted sample variance
#'
#' @keywords internal
#' @importFrom stats weighted.mean
#'
#' @examples
#' \dontrun{
#' values <- c(2.1, 3.5, 1.8, 4.2, 2.9)
#' weights <- c(0.8, 0.3, 0.9, 0.5, 0.7)
#' weighted_var(values, weights)
#' }

weighted_var <- function(x, w, na.rm = FALSE) {
  if (!is.numeric(x) || !is.numeric(w)) {
    stop("x and/or w must be numeric")
  }

  if (length(x) != length(w)) {
    stop("x and w must have the same length")
  }

  if (na.rm) {
    i <- !is.na(x) & !is.na(w)
    x <- x[i]
    w <- w[i]
  }
  if (length(x) < 2L) {
    return(NA_real_)
  }
  sum.w <- sum(w)
  denom <- sum.w^2 - sum(w^2)
  if (sum.w == 0 || denom <= 0) {
    return(NA_real_)
  }
  rtn <- (sum.w * (sum(w * (x - weighted.mean(x, w))^2))) / denom
  return(rtn)
}

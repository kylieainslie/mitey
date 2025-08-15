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
#' @export
#' @importFrom stats weighted.mean
#'
#' @examples
#' # Example 1: Basic weighted variance calculation
#' values <- c(2.1, 3.5, 1.8, 4.2, 2.9)
#' weights <- c(0.8, 0.3, 0.9, 0.5, 0.7)
#' weighted_var(values, weights)
#'
#' # Example 2: Compare with unweighted variance
#' x <- 1:10
#' equal_weights <- rep(1, 10)
#' unweighted_var <- var(x)
#' weighted_var_equal <- weighted_var(x, equal_weights)
#'
#' # Example 3: Using precision weights
#' # Measurements with different levels of confidence
#' measurements <- c(10.2, 9.8, 10.5, 9.9, 10.1)
#' confidence_levels <- c(0.9, 0.6, 0.8, 0.95, 0.7)  # How much we trust each measurement
#'
#' precision_weighted_var <- weighted_var(measurements, confidence_levels)
#'
#' # Example 4: Handling missing values
#' x_with_na <- c(1, 2, NA, 4, 5)
#' weights_with_na <- c(0.2, 0.3, 0.1, 0.8, 0.4)
#'
#' # This will return NA
#' result_na <- weighted_var(x_with_na, weights_with_na, na.rm = FALSE)
#'
#' # This will calculate after removing NA
#' result_removed <- weighted_var(x_with_na, weights_with_na, na.rm = TRUE)
#'
#' # Example 5: Different weight patterns and their effects
#' data_points <- c(10, 15, 20, 25, 30)
#'
#' # Equal weights (should approximate unweighted variance)
#' equal_wts <- rep(1, 5)
#' var_equal <- weighted_var(data_points, equal_wts)
#'
#' # Emphasizing central values
#' central_wts <- c(0.1, 0.3, 1.0, 0.3, 0.1)
#' var_central <- weighted_var(data_points, central_wts)
#'
#' # Emphasizing extreme values
#' extreme_wts <- c(1.0, 0.1, 0.1, 0.1, 1.0)
#' var_extreme <- weighted_var(data_points, extreme_wts)

weighted_var <- function(x, w, na.rm = FALSE) {
  if (!is.numeric(x) | !is.numeric(w)) {
    stop("x and/or w must be numeric")
  }

  if (length(x) != length(w)) {
    stop("x and w must have the same length")
  }

  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  sum.w <- sum(w)
  rtn <- (sum.w * (sum(w * (x - weighted.mean(x, w))^2))) / (sum.w^2 - sum(w^2))
  return(rtn)
}

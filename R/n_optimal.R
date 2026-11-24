#' Heuristics for the Number of Transmission Routes
#'
#' Simple data-driven starting points for choosing \code{n_routes_max} in
#' \code{\link{compare_n_routes}}. They are rough guides, not estimators: use
#' them to set the range of models to compare, then choose between models
#' using AIC/BIC.
#'
#' @param dat numeric vector; index case-to-case (ICC) intervals in days.
#'   Must contain at least two non-missing values.
#'
#' @details
#' \itemize{
#'   \item \code{n_optimal1()} returns the 95th percentile of \code{dat}
#'     divided by its mean, truncated to an integer. \code{dat} must have a
#'     positive mean.
#'   \item \code{n_optimal2()} returns the maximum of \code{dat} divided by
#'     its standard deviation, truncated to an integer. \code{dat} must not be
#'     constant.
#' }
#' Either heuristic can return a value below 2, the smallest valid
#' \code{n_routes_max}; in that case start from 2.
#'
#' @return A single integer.
#'
#' @seealso \code{\link{compare_n_routes}}
#'
#' @name n_optimal
#' @examples
#' set.seed(42)
#' icc_data <- round(pmax(c(
#'   abs(rnorm(30, mean = 0,  sd = 2)),
#'   rnorm(80, mean = 12, sd = 3),
#'   rnorm(30, mean = 24, sd = 4)
#' ), 0))
#'
#' n_optimal1(icc_data)
#' n_optimal2(icc_data)
NULL

#' @rdname n_optimal
#' @export
n_optimal1 <- function(dat) {
  check_n_optimal_dat(dat)
  if (mean(dat) <= 0) {
    stop("dat must have a positive mean.")
  }
  as.integer(quantile(dat, 0.95) / mean(dat))
}

#' @rdname n_optimal
#' @export
n_optimal2 <- function(dat) {
  check_n_optimal_dat(dat)
  if (sd(dat) == 0) {
    stop("dat must not be constant (standard deviation is 0).")
  }
  as.integer(max(dat) / sd(dat))
}

check_n_optimal_dat <- function(dat) {
  if (!is.numeric(dat) || length(dat) < 2) {
    stop("dat must be a numeric vector with at least two values.")
  }
  if (anyNA(dat) || any(!is.finite(dat))) {
    stop("dat must not contain missing or infinite values.")
  }
}

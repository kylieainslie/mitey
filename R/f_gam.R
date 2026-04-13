#' Calculate serial interval mixture density assuming underlying gamma distribution
#'
#' This function computes the weighted mixture density for serial intervals based on
#' different transmission routes in an outbreak. It implements part of the Vink et al.
#' (2014) method for serial interval estimation, assuming an underlying gamma
#' distribution for the serial interval.
#'
#' The function models n_routes distinct transmission routes:
#' \itemize{
#'   \item Co-primary (CP): Cases infected simultaneously from the same source
#'   \item Primary-secondary (PS): Direct transmission from index case
#'   \item Primary-tertiary (PT): Transmission through one intermediate case
#'   \item And so on up to n_routes
#' }
#'
#' @param x quantile or vector of quantiles (time in days since index case symptom onset)
#' @param weights numeric vector of length n_routes - 1; probability weights for each
#'   transmission route starting from Co-primary. The weight for the last route is
#'   derived as 1 - sum(weights).
#' @param mu mean serial interval in days (must be positive)
#' @param sigma standard deviation of serial interval in days (must be positive)
#'
#' @details
#' The weights vector must have length n_routes - 1, with the remaining probability
#' (1 - sum(weights)) assigned to the last route. The function converts mean and
#' standard deviation to gamma distribution shape (k) and scale (theta) parameters:
#' \deqn{k = \mu^2 / \sigma^2}
#' \deqn{\theta = \sigma^2 / \mu}
#' Route i uses shape i * k and scale theta.
#'
#' @returns Vector of weighted density values corresponding to input quantiles x.
#'
#' @references
#' Vink, M. A., Bootsma, M. C. J., & Wallinga, J. (2014). Serial intervals of
#' respiratory infectious diseases: A systematic review and analysis.
#' American Journal of Epidemiology, 180(9), 865-875.
#'
#' @seealso \code{\link{si_estim}}, \code{\link{plot_si_fit}}, \code{\link{f_norm}}
#' @keywords internal
#' @importFrom stats dgamma
#' @examples
#' \dontrun{
#' # 4 routes (default behaviour)
#' x <- seq(0.1, 30, by = 0.1)
#' density_values <- f_gam(x, weights = c(0.1, 0.6, 0.2), mu = 6.5, sigma = 2.8)
#' plot(x, density_values, type = "l")
#'
#' # 5 routes
#' density_values5 <- f_gam(x, weights = c(0.1, 0.6, 0.15, 0.10), mu = 6.5, sigma = 2.8)
#' plot(x, density_values5, type = "l")
#' }
#'
f_gam <- function(
    x,
    weights,
    mu,
    sigma
) {
  n_routes <- length(weights) + 1L

  if (any(weights < 0)) {
    stop("All weights must be non-negative.")
  }
  if (sum(weights) > 1) {
    stop("Sum of weights must not exceed 1.")
  }

  k     <- (mu^2) / (sigma^2)
  theta <- (sigma^2) / mu

  # Last route weight derived from the others
  w_last <- 1 - sum(weights)

  # Component 1: Co-primary (special Bessel case)
  term1 <- weights[1] /
    sqrt(pi) *
    2^(3/2 - k) *
    theta^(-0.5 - k) *
    x^(-0.5 + k) *
    besselK(x / theta, 0.5 - k) *
    1 / gamma(k)
  result <- term1

  # Routes 2 to n_routes - 1
  if (n_routes >= 3) {
    for (i in 2:(n_routes - 1)) {
      result <- result + weights[i] * dgamma(x, shape = i * k, scale = theta)
    }
  }

  # Last route
  result <- result + w_last * dgamma(x, shape = (n_routes - 1) * k, scale = theta)

  return(result)
}

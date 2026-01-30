#' Calculate serial interval mixture density assuming underlying gamma distribution
#'
#' This function computes the weighted mixture density for serial intervals based on
#' different transmission routes in an outbreak. It implements part of the Vink et al.
#' (2014) method for serial interval estimation, assuming an underlying gamma
#' distribution for the serial interval.
#'
#' The function models four distinct transmission routes:
#' \itemize{
#'   \item Co-primary (CP): Cases infected simultaneously from the same source
#'   \item Primary-secondary (PS): Direct transmission from index case
#'   \item Primary-tertiary (PT): Transmission through one intermediate case
#'   \item Primary-quaternary (PQ): Transmission through two intermediate cases
#' }
#'
#' Each route contributes to the overall serial interval distribution with different
#' means and variances. The co-primary component uses a modified gamma distribution
#' to account for simultaneous infections, while subsequent generations follow
#' gamma distributions with progressively longer means and larger variances.
#'
#' This function is primarily used internally by \code{\link{si_estim}} when
#' \code{dist = "gamma"} is specified, and by \code{\link{plot_si_fit}} for
#' visualizing fitted distributions.
#'
#' @param x quantile or vector of quantiles (time in days since index case symptom onset)
#' @param w1 probability weight of being a co-primary case
#' @param w2 probability weight of being a primary-secondary case
#' @param w3 probability weight of being a primary-tertiary case
#' @param mu mean serial interval in days (must be positive)
#' @param sigma standard deviation of serial interval in days (must be positive)
#'
#' @details
#' The weights w1, w2, and w3 must sum to <= 1, with the remaining probability
#' (1 - w1 - w2 - w3) assigned to primary-quaternary cases. The function converts
#' the mean and standard deviation to gamma distribution shape (k) and scale (\\theta)
#' parameters using the method of moments:
#' \deqn{k = \mu^2 / \sigma^2}
#' \deqn{\theta = \sigma^2 / \mu}
#'
#' @returns Vector of weighted density values corresponding to input quantiles x.
#'   Returns the sum of densities from all four transmission routes.
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
#' x <- seq(0.1, 30, by = 0.1)
#' density_values <- f_gam(x, w1 = 0.1, w2 = 0.6, w3 = 0.2, mu = 6.5, sigma = 2.8)
#' plot(x, density_values, type = "l")
#' }
#'
f_gam <- function(
  x,
  w1,
  w2,
  w3,
  mu,
  sigma
) {
  k <- (mu^2) / (sigma^2)
  theta <- (sigma^2) / mu

  term1 <- w1 /
    sqrt(pi) *
    2^(3 / 2 - k) *
    theta^(-0.5 - k) *
    x^(-0.5 + k) *
    besselK(x / theta, 0.5 - k) *
    1 /
    gamma(k)
  term2 <- w2 * dgamma(x, k, scale = theta)
  term3 <- w3 * dgamma(x, 2 * k, scale = theta)
  term4 <- (1 - w1 - w2 - w3) * dgamma(x, 3 * k, scale = theta)

  return(term1 + term2 + term3 + term4)
}

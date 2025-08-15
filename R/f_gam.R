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
#' @param w1 probability weight of being a co-primary case (0 ≤ w1 ≤ 1)
#' @param w2 probability weight of being a primary-secondary case (0 ≤ w2 ≤ 1)
#' @param w3 probability weight of being a primary-tertiary case (0 ≤ w3 ≤ 1)
#' @param mu mean serial interval in days (must be positive)
#' @param sigma standard deviation of serial interval in days (must be positive)
#'
#' @details
#' The weights w1, w2, and w3 must sum to ≤ 1, with the remaining probability
#' (1 - w1 - w2 - w3) assigned to primary-quaternary cases. The function converts
#' the mean and standard deviation to gamma distribution shape (k) and scale (θ)
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
#' @export
#' @importFrom stats dgamma
#' @examples
#' # Example: Plot serial interval mixture density for influenza-like outbreak
#'
#' # Set parameters for a typical respiratory infection
#' mu <- 6.5      # Mean serial interval of 6.5 days
#' sigma <- 2.8   # Standard deviation of 2.8 days
#'
#' # Set transmission route weights
#' w1 <- 0.1      # 10% co-primary cases
#' w2 <- 0.6      # 60% primary-secondary cases
#' w3 <- 0.2      # 20% primary-tertiary cases
#' # Remaining 10% are primary-quaternary cases (1 - w1 - w2 - w3 = 0.1)
#'
#' # Create sequence of time points
#' x <- seq(0.1, 30, by = 0.1)
#'
#' # Calculate mixture density
#' density_values <- f_gam(x, w1, w2, w3, mu, sigma)
#'
#' # Plot the result
#' plot(x, density_values, type = "l", lwd = 2, col = "red",
#'      xlab = "Days", ylab = "Density",
#'      main = "Serial Interval Mixture Density (Gamma Distribution)")
#' grid()
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

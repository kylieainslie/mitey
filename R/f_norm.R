#' Calculate serial interval mixture density assuming underlying normal distribution
#'
#' This function computes the weighted mixture density for serial intervals based on
#' different transmission routes in an outbreak. It implements part of the Vink et al.
#' (2014) method for serial interval estimation, assuming an underlying normal
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
#' means and variances. The co-primary component uses a half-normal distribution
#' to model simultaneous infections (preventing negative serial intervals), while
#' subsequent generations follow normal distributions with means that are multiples
#' of the base serial interval.
#'
#' This function is primarily used internally by \code{\link{si_estim}} when
#' \code{dist = "normal"} is specified (the default), and by \code{\link{plot_si_fit}}
#' for visualizing fitted distributions. The normal distribution assumption allows
#' for negative serial intervals, which may be more realistic for some pathogens.
#'
#' @param x quantile or vector of quantiles (time in days since index case symptom onset)
#' @param w1 probability weight of being a co-primary case (0 ≤ w1 ≤ 1)
#' @param w2 probability weight of being a primary-secondary case (0 ≤ w2 ≤ 1)
#' @param w3 probability weight of being a primary-tertiary case (0 ≤ w3 ≤ 1)
#' @param mu mean serial interval in days (can be any real number)
#' @param sigma standard deviation of serial interval in days (must be positive)
#'
#' @details
#' The weights w1, w2, and w3 must sum to ≤ 1, with the remaining probability
#' (1 - w1 - w2 - w3) assigned to primary-quaternary cases. The transmission
#' route distributions are parameterized as:
#' \itemize{
#'   \item Co-primary: Half-normal with scale parameter derived from sigma
#'   \item Primary-secondary: Normal(μ, σ)
#'   \item Primary-tertiary: Normal(2μ, √2σ)
#'   \item Primary-quaternary: Normal(3μ, √3σ)
#' }
#'
#' @returns Vector of weighted density values corresponding to input quantiles x.
#'   Returns the sum of densities from all four transmission routes.
#'
#' @references
#' Vink, M. A., Bootsma, M. C. J., & Wallinga, J. (2014). Serial intervals of
#' respiratory infectious diseases: A systematic review and analysis.
#' American Journal of Epidemiology, 180(9), 865-875.
#'
#' @seealso \code{\link{si_estim}}, \code{\link{plot_si_fit}}, \code{\link{f_gam}}
#' @export
#' @importFrom stats dnorm
#' @importFrom fdrtool dhalfnorm
#' @examples
#' # Example: Plot serial interval mixture density for scabies outbreak
#'
#' # Set parameters based on scabies epidemiology (longer serial interval)
#' mu <- 123     # Mean serial interval of 123 days (from Ainslie et al.)
#' sigma <- 32   # Standard deviation of 32 days
#'
#' # Set transmission route weights typical for scabies
#' w1 <- 0.15    # 15% co-primary cases
#' w2 <- 0.50    # 50% primary-secondary cases
#' w3 <- 0.25    # 25% primary-tertiary cases
#' # Remaining 10% are primary-quaternary cases (1 - w1 - w2 - w3 = 0.1)
#'
#' # Create sequence of time points
#' x <- seq(0, 400, by = 1)
#'
#' # Calculate mixture density
#' density_values <- f_norm(x, w1, w2, w3, mu, sigma)
#'
#' # Plot the result
#' plot(x, density_values, type = "l", lwd = 2, col = "red",
#'      xlab = "Days", ylab = "Density",
#'      main = "Serial Interval Mixture Density (Normal Distribution)")
#' grid()
#'
f_norm <- function(
  x,
  w1,
  w2,
  w3,
  mu,
  sigma
) {
  term1 <- w1 * dhalfnorm(x, sqrt(pi / 2) / (sqrt(2) * sigma))
  term2 <- w2 * dnorm(x, mean = mu, sd = sigma)
  term3 <- w3 * dnorm(x, mean = 2 * mu, sd = sqrt(2) * sigma)
  term4 <- (1 - w1 - w2 - w3) * dnorm(x, mean = 3 * mu, sd = sqrt(3) * sigma)

  return(term1 + term2 + term3 + term4)
}

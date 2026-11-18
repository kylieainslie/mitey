
#' Calculate serial interval mixture density assuming underlying lognormal distribution
#'
#' This function computes the weighted mixture density for serial intervals based on
#' different transmission routes in an outbreak. It extends the Vink et al.
#' (2014) method for serial interval estimation by assuming an underlying lognormal
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
#' Each route contributes to the overall mixture density with different
#' distributional characteristics. The co-primary component uses a half-normal
#' approximation, the primary-secondary component follows the lognormal
#' distribution, and the primary-tertiary and primary-quaternary components
#' use Fenton-Wilkinson approximations for sums of independent lognormal
#' serial intervals.
#'
#' This function is used by \code{\link{plot_si_fit}} for visualizing the
#' fitted continuous mixture density when \code{dist = "lognormal"}.
#' It is not called directly by the EM algorithm in \code{\link{si_estim}},
#' which evaluates discretised component probabilities separately.
#'
#' @param x quantile or vector of quantiles (time in days since index case symptom onset)
#' @param w1 probability weight of being a co-primary case
#' @param w2 probability weight of being a primary-secondary case
#' @param w3 probability weight of being a primary-tertiary case
#' @param mu mean serial interval in days (must be positive)
#' @param sigma standard deviation of serial interval in days (must be positive)
#'
#' @details
#' The weights w1, w2, and w3 must sum to <= 1, with the
#' remaining probability (1 - w1 - w2 - w3) assigned to primary-quaternary cases.
#' The mean and standard deviation are specified on the natural scale.
#' For CP,  a half-normal approximation is used with scale determined by sigma.
#' For PS, the parameters are converted to the log scale:
#' \deqn{sdlog = \sqrt{\log(1 + \sigma^2 / \mu^2)}}
#' \deqn{meanlog = \log(\mu) - sdlog^2 / 2}
#' For PT and PQ, the Fenton-Wilkinson approximation matches the
#' moments of sums of two and three independent lognormal intervals,
#' respectively.
#'
#' @returns Vector of weighted density values corresponding to input quantiles x.
#'   Returns the sum of densities from all four transmission routes.
#'
#' @references
#' Vink, M. A., Bootsma, M. C. J., & Wallinga, J. (2014). Serial intervals of
#' respiratory infectious diseases: A systematic review and analysis.
#' American Journal of Epidemiology, 180(9), 865-875.
#'
#' @seealso \code{\link{si_estim}}, \code{\link{plot_si_fit}},
#'          \code{\link{f_norm}}, \code{\link{f_gam}}
#' @keywords internal
#' @importFrom stats dlnorm
#' @importFrom fdrtool dhalfnorm
#' @examples
#' \dontrun{
#' x <- seq(0, 60, by = 0.1)
#' density_values <- f_lnorm(x, w1 = 0.2, w2 = 0.5, w3 = 0.2, mu = 15, sigma = 3)
#' plot(x, density_values, type = "l")
#' }
#'

f_lnorm <- function(
    x,
    w1,
    w2,
    w3,
    mu,
    sigma
) {
  if (mu <= 0 || sigma <= 0 || !is.finite(mu) || !is.finite(sigma)) {
    return(rep(0, length(x)))
  }
  
  sdlog <- sqrt(log(1 + sigma^2 / mu^2))
  meanlog <- log(mu) - 0.5 * sdlog^2
  
  sdlog_pt <- sqrt(log(1 + sigma^2 / (2 * mu^2)))
  meanlog_pt <- log(2 * mu) - 0.5 * sdlog_pt^2
  
  sdlog_pq <- sqrt(log(1 + sigma^2 / (3 * mu^2)))
  meanlog_pq <- log(3 * mu) - 0.5 * sdlog_pq^2
  
  term1 <- w1 * ifelse(
    x >= 0,
    dhalfnorm(x, theta = sqrt(pi / 2) / (sqrt(2) * sigma)),
    0
  )
  term2 <- w2 * ifelse(x > 0, dlnorm(x, meanlog = meanlog, sdlog = sdlog), 0)
  term3 <- w3 * ifelse(x > 0, dlnorm(x, meanlog = meanlog_pt, sdlog = sdlog_pt), 0)
  term4 <- (1 - w1 - w2 - w3) *
    ifelse(x > 0, dlnorm(x, meanlog = meanlog_pq, sdlog = sdlog_pq), 0)
  
  term1 + term2 + term3 + term4
}
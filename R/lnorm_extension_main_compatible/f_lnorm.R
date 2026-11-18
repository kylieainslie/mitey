#' Calculate serial interval mixture density assuming underlying lognormal distribution
#'
#' This function computes the weighted mixture density for serial intervals based on
#' different transmission routes in an outbreak. It extends the Vink et al.
#' (2014) method by assuming an underlying lognormal distribution for the serial interval.
#'
#' The function models n_routes distinct transmission routes:
#' \itemize{
#'   \item Co-primary (CP): Cases infected simultaneously from the same source
#'   \item Primary-secondary (PS): Direct transmission from index case
#'   \item Primary-tertiary (PT): Transmission through one intermediate case
#'   \item And so on up to n_routes
#' }
#'
#' The co-primary component uses a half-normal approximation, the primary-secondary
#' component follows the exact lognormal distribution, and higher transmission routes
#' use Fenton-Wilkinson approximations for sums of independent lognormal intervals.
#'
#' This function is used by \code{\link{plot_si_fit}} to visualise the fitted
#' continuous mixture density. The EM algorithm in \code{\link{si_estim}}
#' evaluates discretised component probabilities separately.
#'
#' @param x quantile or vector of quantiles (time in days since index case symptom onset)
#' @param weights numeric vector of length n_routes - 1; probability weights for each
#'   transmission route starting from Co-primary. The weight for the last route is
#'   derived as 1 - sum(weights).
#' @param mu mean serial interval in days (must be positive)
#' @param sigma standard deviation of serial interval in days (must be positive)
#'
#' @details
#' The mean and standard deviation are specified on the natural scale. Co-primary
#' uses a half-normal approximation with scale determined by sigma. For a positive
#' route consisting of i serial intervals, the lognormal approximation matches
#' the natural-scale mean i * mu and variance i * sigma^2. The case i = 1 (PS)
#' is exact; i >= 2 uses the Fenton-Wilkinson approximation.
#'
#' @returns Vector of weighted density values corresponding to input quantiles x.
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
#' # 4 routes (default behaviour)
#' x <- seq(0, 60, by = 0.1)
#' density_values <- f_lnorm(x, weights = c(0.2, 0.5, 0.2), mu = 15, sigma = 3)
#' plot(x, density_values, type = "l")
#'
#' # 5 routes
#' density_values5 <- f_lnorm(x, weights = c(0.2, 0.5, 0.15, 0.1), mu = 15, sigma = 3)
#' plot(x, density_values5, type = "l")
#' }
#'
f_lnorm <- function(
    x,
    weights,
    mu,
    sigma
) {
  n_routes <- length(weights) + 1L

  if (length(weights) < 1L || any(!is.finite(weights))) {
    stop("Weights must be a non-empty vector of finite values.")
  }
  if (any(weights < 0)) {
    stop("All weights must be non-negative.")
  }
  if (sum(weights) > 1) {
    stop("Sum of weights must not exceed 1.")
  }
  if (!is.finite(mu) || !is.finite(sigma) || mu <= 0 || sigma <= 0) {
    return(rep(0, length(x)))
  }

  # Last route weight derived from the others
  w_last <- 1 - sum(weights)

  # Component 1: Co-primary (half-normal approximation)
  result <- weights[1] * ifelse(
    x >= 0,
    dhalfnorm(x, theta = sqrt(pi / 2) / (sqrt(2) * sigma)),
    0
  )

  # Routes 2 to n_routes - 1
  if (n_routes >= 3) {
    for (i in 2:(n_routes - 1)) {
      k <- i - 1
      sdlog <- sqrt(log1p(sigma^2 / (k * mu^2)))
      meanlog <- log(k * mu) - 0.5 * sdlog^2
      result <- result + weights[i] * dlnorm(x, meanlog = meanlog, sdlog = sdlog)
    }
  }

  # Last route
  k <- n_routes - 1
  sdlog <- sqrt(log1p(sigma^2 / (k * mu^2)))
  meanlog <- log(k * mu) - 0.5 * sdlog^2
  result <- result + w_last * dlnorm(x, meanlog = meanlog, sdlog = sdlog)

  return(result)
}

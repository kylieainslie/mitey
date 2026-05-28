#' Calculate serial interval mixture density assuming underlying normal distribution
#'
#' This function computes the weighted mixture density for serial intervals based on
#' different transmission routes in an outbreak. It implements part of the Vink et al.
#' (2014) method for serial interval estimation, assuming an underlying normal
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
#' @param mu mean serial interval in days (can be any real number)
#' @param sigma standard deviation of serial interval in days (must be positive)
#'
#' @details
#' The weights vector must have length n_routes - 1, with the remaining probability
#' (1 - sum(weights)) assigned to the last route. The transmission route distributions
#' are parameterized as:
#' \itemize{
#'   \item Co-primary: Half-normal with scale parameter derived from sigma
#'   \item Route i (i >= 2): Normal(i * mu, sqrt(i) * sigma) and Normal(-i * mu, sqrt(i) * sigma)
#' }
#'
#' @returns Vector of weighted density values corresponding to input quantiles x.
#'
#' @references
#' Vink, M. A., Bootsma, M. C. J., & Wallinga, J. (2014). Serial intervals of
#' respiratory infectious diseases: A systematic review and analysis.
#' American Journal of Epidemiology, 180(9), 865-875.
#'
#' @seealso \code{\link{si_estim}}, \code{\link{plot_si_fit}}, \code{\link{f_gam}}
#' @keywords internal
#' @importFrom stats dnorm
#' @importFrom fdrtool dhalfnorm
#' @examples
#' \dontrun{
#' # 4 routes (default behaviour)
#' x <- seq(0, 400, by = 1)
#' density_values <- f_norm(x, weights = c(0.15, 0.50, 0.25), mu = 123, sigma = 32)
#' plot(x, density_values, type = "l")
#'
#' # 5 routes
#' density_values5 <- f_norm(x, weights = c(0.15, 0.50, 0.20, 0.10), mu = 123, sigma = 32)
#' plot(x, density_values5, type = "l")
#' }
#'
f_norm <- function(
    x,
    weights,
    mu,
    sigma,
    n_routes
) {

  if (any(weights < 0)) {
    stop("All weights must be non-negative.")
  }
  if (sum(weights) > 1  + .Machine$double.eps * 100) {
    stop("Sum of weights must not exceed 1.")
  }

  if (length(weights)!=n_routes -1) {
    stop("f norm takes n_routes - 1 weights, more have been given ")
  }


  # Component 1: Co-primary (half-normal, special case)
  result <- weights[1] * dhalfnorm(x, sqrt(pi / 2) / (sqrt(2) * sigma))

  if (n_routes >2){
    for (i in 2:(n_routes - 1)) {
      result <- result +
        weights[i] * dnorm(x, mean =  (i-1) * mu, sd = sqrt(i-1) * sigma)
    }

  }




  last_weight<-1-sum(weights)
  result <- result + last_weight*dnorm(x, mean =  (n_routes -1) * mu, sd = sqrt(n_routes -1) * sigma)

  return(result)
}

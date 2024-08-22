#' Helper function to plot the fit of the serial interval distribution assuming an underlying gamma distribution
#'
#' @param x quantile
#' @param w1 probability weight of being a co-primary case
#' @param w2 probability weight of being a primary-secondary case
#' @param w3 probability weight of being a primary-tertiary case
#' @param mu mean serial interval
#' @param sigma standard deviation of serial interval
#' @returns the weighted density for each value of x
#' @export
#' @importFrom stats dgamma
#'
f_gam <- function(x, w1, w2, w3, mu, sigma) {

  k <- (mu^2) / (sigma^2)
  theta <- (sigma^2) / mu

  term1 <- w1 / sqrt(pi) * 2^(3/2 - k) * theta^(-0.5 - k) * x^(-0.5 + k) * besselK(x / theta, 0.5 - k) * 1 / gamma(k)
  term2 <- w2 * dgamma(x, k, scale = theta)
  term3 <- w3 * dgamma(x, 2 * k, scale = theta)
  term4 <- (1 - w1 - w2 - w3) * dgamma(x, 3 * k, scale = theta)

  return(term1 + term2 + term3 + term4)

}

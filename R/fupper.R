#' Calculate fupper for Different Components
#'
#' This function calculates the value of fupper based on the component.
#'
#' @param x The value at which to evaluate the function.
#' @param r The value of r.
#' @param mu The mean value.
#' @param sigma The standard deviation.
#' @param comp The component number (1 to 7).
#'
#' @return The calculated value of fupper.
#' @export
fupper <- function(x, r, mu, sigma, comp) {
  if (comp == 1) return((r + 1 - x) * dhalfnorm(x, theta = sqrt(pi / 2) / (sqrt(2) * sigma)))
  if (comp == 2) return((r + 1 - x) * dnorm(x, mean = mu, sd = sigma))
  if (comp == 3) return((r + 1 - x) * dnorm(x, mean = -mu, sd = sigma))
  if (comp == 4) return((r + 1 - x) * dnorm(x, mean = 2 * mu, sd = sqrt(2) * sigma))
  if (comp == 5) return((r + 1 - x) * dnorm(x, mean = -2 * mu, sd = sqrt(2) * sigma))
  if (comp == 6) return((r + 1 - x) * dnorm(x, mean = 3 * mu, sd = sqrt(3) * sigma))
  if (comp == 7) return((r + 1 - x) * dnorm(x, mean = -3 * mu, sd = sqrt(3) * sigma))
}

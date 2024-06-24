#' Calculate flower for Different Components
#'
#' This function calculates the value of flower based on the component.
#'
#' @param x The value at which to evaluate the function.
#' @param r The value of r.
#' @param mu The mean value.
#' @param sigma The standard deviation.
#' @param comp The component number (1 to 7).
#'
#' @return The calculated value of flower.
#' @export
flower <- function(x, r, mu, sigma, comp) {
  if (comp == 1) return((x - r + 1) * dhalfnorm(x, theta = sqrt(pi / 2) / (sqrt(2) * sigma)))
  if (comp == 2) return((x - r + 1) * dnorm(x, mean = mu, sd = sigma))
  if (comp == 3) return((x - r + 1) * dnorm(x, mean = -mu, sd = sigma))
  if (comp == 4) return((x - r + 1) * dnorm(x, mean = 2 * mu, sd = sqrt(2) * sigma))
  if (comp == 5) return((x - r + 1) * dnorm(x, mean = -2 * mu, sd = sqrt(2) * sigma))
  if (comp == 6) return((x - r + 1) * dnorm(x, mean = 3 * mu, sd = sqrt(3) * sigma))
  if (comp == 7) return((x - r + 1) * dnorm(x, mean = -3 * mu, sd = sqrt(3) * sigma))
}

#' Calculate f0 for Different Components
#'
#' This function calculates the value of f0 based on the component.
#'
#' @param x The value at which to evaluate the function.
#' @param mu The mean value.
#' @param sigma The standard deviation.
#' @param comp The component number (1 to 7).
#'
#' @return The calculated value of f0.
#' @export
f0 <- function(x, mu, sigma, comp) {
  if (comp == 1) return((2 - 2 * x) * dhalfnorm(x, theta = sqrt(pi / 2) / (sqrt(2) * sigma)))
  if (comp == 2) return((2 - 2 * x) * dnorm(x, mean = mu, sd = sigma))
  if (comp == 3) return((2 - 2 * x) * dnorm(x, mean = -mu, sd = sigma))
  if (comp == 4) return((2 - 2 * x) * dnorm(x, mean = 2 * mu, sd = sqrt(2) * sigma))
  if (comp == 5) return((2 - 2 * x) * dnorm(x, mean = -2 * mu, sd = sqrt(2) * sigma))
  if (comp == 6) return((2 - 2 * x) * dnorm(x, mean = 3 * mu, sd = sqrt(3) * sigma))
  if (comp == 7) return((2 - 2 * x) * dnorm(x, mean = -3 * mu, sd = sqrt(3) * sigma))
}

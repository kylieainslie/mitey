#' Calculate flower for Different Components
#'
#' This function calculates the value of flower based on the component.
#'
#' @param x The value at which to evaluate the function.
#' @param r The value of r.
#' @param mu The mean value.
#' @param sigma The standard deviation.
#' @param comp The component number (1 to 7).
#' @param dist string; assumed distribution of the serial interval; accepts "normal" or "gamma"; defaults to "normal"
#'
#' @return The calculated value of flower.
#' @export
flower <- function(x, r, mu, sigma, comp, dist = "normal") {

  # error messages
  if(dist != "normal" && dist != "gamma"){
    stop("Incorrect distribution specified. Acceptable arguments are c('normal', 'gamma')")
  }

  if (dist == "normal"){
    if (comp == 1) return((x - r + 1) * dhalfnorm(x, theta = sqrt(pi / 2) / (sqrt(2) * sigma)))
    if (comp == 2) return((x - r + 1) * dnorm(x, mean = mu, sd = sigma))
    if (comp == 3) return((x - r + 1) * dnorm(x, mean = -mu, sd = sigma))
    if (comp == 4) return((x - r + 1) * dnorm(x, mean = 2 * mu, sd = sqrt(2) * sigma))
    if (comp == 5) return((x - r + 1) * dnorm(x, mean = -2 * mu, sd = sqrt(2) * sigma))
    if (comp == 6) return((x - r + 1) * dnorm(x, mean = 3 * mu, sd = sqrt(3) * sigma))
    if (comp == 7) return((x - r + 1) * dnorm(x, mean = -3 * mu, sd = sqrt(3) * sigma))

  } else if (dist == "gamma"){

    # convert mean and sd for normal distn into shape and scale parameters for gamma distn.
    k <- (mu^2) / (sigma^2)
    theta <- (sigma^2) / mu

    if (comp == 1) return((x - r + 1) * 1/sqrt(pi) * 2^(3/2-k) * (theta)^(-0.5-k) * x^(-0.5+k) * besselK(x/(theta),0.5-k) * 1/gamma(k))
    if (comp %in% c(2, 3)) return((x - r + 1) * dgamma(x, shape = k, scale = theta))
    if (comp %in% c(4, 5)) return((x - r + 1) * dgamma(x, shape = 2*k, scale = theta))
    if (comp %in% c(6, 7)) return((x - r + 1) * dgamma(x, shape = 3*k, scale = theta))
  }
}

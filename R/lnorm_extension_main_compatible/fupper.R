#' Calculate fupper for Different Components
#'
#' This function calculates the value of fupper based on the component.
#' For lognormal, CP uses a Half-Normal approximation and positive routes use
#' Lognormal densities with Fenton-Wilkinson moment matching.
#'
#' @param x The value at which to evaluate the function.
#' @param r The value of r.
#' @param mu The mean value.
#' @param sigma The standard deviation.
#' @param comp integer; the component number. Component 1 is Co-Primary. Even components
#'   2i are positive routes, odd components 2i+1 are negative routes (normal only).
#' @param dist string; assumed distribution of the serial interval; accepts "normal", "gamma" or "lognormal"; defaults to "normal".
#' @param wind The window censure interval
#'
#' @return The calculated value of fupper.
#' @keywords internal
#' @examples
#' \dontrun{
#' fupper(x = 15, r = 20, mu = 12, sigma = 3, comp = 2, dist = "normal")
#' fupper(x = 15, r = 20, mu = 12, sigma = 3, comp = 2, dist = "gamma")
#' }
fupper <- function(
    x,
    r,
    mu,
    sigma,
    comp,
    dist = "normal",
    wind = 1
) {

  # error messages
  if (!dist %in% c("normal", "gamma", "lognormal")) {
    stop(
      "Incorrect distribution specified. Acceptable arguments are c('normal', 'gamma', 'lognormal')"
    )
  }

  if (dist == "lognormal") {
    if (!is.finite(mu) || !is.finite(sigma) || mu <= 0 || sigma <= 0) {
      return(rep(0, length(x)))
    }
    if (comp != 1 && (comp < 2 || comp %% 2 != 0)) {
      stop("Lognormal only supports component 1 and positive even components.")
    }
  }

  # Component 1: Co-Primary route (special case, unchanged)
  if (comp == 1) {
    if (dist %in% c("normal", "lognormal")) {
      result <- (r/wind + 1 - x/wind) * dhalfnorm(x, theta = sqrt(pi / 2) / (sqrt(2) * sigma))
      result[!is.finite(result)] <- 0
      return(result)

    } else {
      k <- (mu^2) / (sigma^2)
      theta <- (sigma^2) / mu
      if (k <= 0 || theta <= 0) return(0)
      bessel_val <- besselK(x / theta, 0.5 - k)
      bessel_val[!is.finite(bessel_val)] <- 0
      return_val <- (r/wind + 1 - x/wind) *
        1 / sqrt(pi) *
        2^(3/2 - k) *
        theta^(-0.5 - k) *
        x^(-0.5 + k) *
        bessel_val *
        1 / gamma(k)
      return_val[is.nan(return_val)] <- 0
      return(return_val)
    }
  }

  # Components 2 and above: generic route logic
  if (comp %% 2 == 0) {
    i <- comp / 2
    route_mean <- i * mu
  } else {
    i <- (comp - 1) / 2
    route_mean <- -i * mu
  }
  route_sd <- sqrt(i) * sigma

  if (dist == "normal") {
    return((r/wind + 1 - x/wind) * dnorm(x, mean = route_mean, sd = route_sd))
  } else if (dist == "lognormal") {
    # Moment-matched Lognormal density for PS and higher routes
    sdlog <- sqrt(log1p((route_sd / route_mean)^2))
    meanlog <- log(route_mean) - 0.5 * sdlog^2
    return((r/wind + 1 - x/wind) *
      dlnorm(x, meanlog = meanlog, sdlog = sdlog))
  } else {
    k <- (mu^2) / (sigma^2)
    theta <- (sigma^2) / mu
    if (k <= 0 || theta <= 0) return(0)
    return((r/wind + 1 - x/wind) * dgamma(x, shape = i * k, scale = theta))
  }
}

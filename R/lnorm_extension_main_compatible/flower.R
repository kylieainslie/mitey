#' Calculate flower for Different Components
#'
#' This function calculates the value of flower based on the component.
#' For lognormal, CP uses a Half-Normal approximation and positive routes use
#' Lognormal densities with Fenton-Wilkinson moment matching.
#'
#' @param x The value at which to evaluate the function.
#' @param r The value of r.
#' @param mu The mean value.
#' @param sigma The standard deviation.
#' @param comp integer; the component number. Component 1 is Co-Primary. Even components
#'   2i are positive routes, odd components 2i+1 are negative routes (normal only).
#' @param dist string; assumed distribution of the serial interval; accepts "normal", "gamma" or "lognormal"; defaults to "normal"
#' @param wind The window censure interval
#'.
#' @return The calculated value of flower.
#' @keywords internal
#' @examples
#' \dontrun{
#' flower(x = 15, r = 10, mu = 12, sigma = 3, comp = 2, dist = "normal")
#' flower(x = 15, r = 10, mu = 12, sigma = 3, comp = 2, dist = "gamma")
#' }
flower <- function(
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
    result <- numeric(length(x))

    if (dist %in% c("normal", "lognormal")) {
      pos <- x > 0   # strictly positive — excludes x = 0
      if (any(pos)) {
        result[pos] <- (x[pos]/wind - r/wind + 1) *
          dhalfnorm(x[pos], theta = sqrt(pi / 2) / (sqrt(2) * sigma))
      }

      return(result)
    }


    else {
      k <- (mu^2) / (sigma^2)
      theta <- (sigma^2) / mu
      if (k <= 0 || theta <= 0) return(rep(0, length(x)))
      pos <- x > 0
      if (any(pos)) {
        bessel_val <- besselK(x[pos] / theta, 0.5 - k)
        bessel_val[!is.finite(bessel_val)] <- 0
        result[pos] <- (x[pos]/wind - r/wind + 1) *
          1 / sqrt(pi) *
          2^(3/2 - k) *
          theta^(-0.5 - k) *
          x[pos]^(-0.5 + k) *
          bessel_val *
          1 / gamma(k)
        result[is.nan(result)] <- 0
      }
      result[!is.finite(result)] <- 0
      return(result)
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
    result <- numeric(length(x))

    pos <- x > 0
    if (any(pos)) {
      result[pos] <- (x[pos]/wind - r/wind + 1) *
        dnorm(x[pos], mean = route_mean, sd = route_sd)
    }
    result[!is.finite(result)] <- 0
    return(result)


  } else if (dist == "lognormal") {
    # Moment-matched Lognormal density for PS and higher routes
    sdlog <- sqrt(log1p((route_sd / route_mean)^2))
    meanlog <- log(route_mean) - 0.5 * sdlog^2
    return((x/wind - r/wind + 1) *
      dlnorm(x, meanlog = meanlog, sdlog = sdlog))
  } else {
    k <- (mu^2) / (sigma^2)
    theta <- (sigma^2) / mu
    if (k <= 0 || theta <= 0) return(0)
    return((x/wind - r/wind + 1) * dgamma(x, shape = i * k, scale = theta))
  }
}

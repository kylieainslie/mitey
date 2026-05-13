#' Calculate flower for Different Components
#'
#' This function calculates the value of flower based on the component.
#'
#' @param x The value at which to evaluate the function.
#' @param r The value of r.
#' @param mu The mean value.
#' @param sigma The standard deviation.
#' @param comp integer; the component number. Component 1 is Co-Primary. Even components
#'   2i are positive routes, odd components 2i+1 are negative routes (normal only).
#' @param dist string; assumed distribution of the serial interval; accepts "normal" or "gamma"; defaults to "normal"
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
  if (dist != "normal" && dist != "gamma") {
    stop(
      "Incorrect distribution specified. Acceptable arguments are c('normal', 'gamma')"
    )
  }

  # Component 1: Co-Primary route (special case, unchanged)
  if (comp == 1) {
    result <- numeric(length(x))

    if(dist == "normal") {
      pos <- x > 0
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
    return(result)


  } else {
    k <- (mu^2) / (sigma^2)
    theta <- (sigma^2) / mu
    if (k <= 0 || theta <= 0) return(0)
    return((x/wind - r/wind + 1) * dgamma(x, shape = i * k, scale = theta))
  }
}

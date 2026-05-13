#' Calculate f0 for Different Components
#'
#' This function calculates the value of f0 based on the component, where the components
#' represent the transmission routes: Co-Primary (CP), Primary-Secondary (PS), Primary-Tertiary (PT), and Primary-Quaternary (PQ) and beyond. We split the PS, PT, PQ (and higher) routes into two parts, such that
#'  - component 1: CP route
#'  - components 2+3: PS route
#'  - components 4+5: PT route
#'  - components 2i + (2i+1): route i+1
#'
#'  If the dist = gamma, then the mean (\eqn{\mu}) and standard deviation (sigma) are converted into the
#'  shape (k) and scale (theta) parameters for the gamma distribution, such that the mean (\eqn{\mu}
#'  ) and variance (\eqn{\sigma^2}) are given by:
#' \deqn{\mu = k \times \theta}
#' \deqn{\sigma^2 = k \times \theta^2}.
#'
#' @param x numeric; the value at which to evaluate the function.
#' @param mu numeric; the mean value.
#' @param sigma numeric; the standard deviation.
#' @param comp integer; the component number. Component 1 is Co-Primary. Even components
#'   2i are positive routes, odd components 2i+1 are negative routes (normal only).
#' @param dist string; assumed distribution of the serial interval; takes "normal" or "gamma"; defaults to "normal"
#' @param wind The window censure interval
#'
#' @return The calculated value of f0.
#' @keywords internal
#' @examples
#' \dontrun{
#' # Basic example with normal distribution
#' f0(x = 0.5, mu = 12, sigma = 3, comp = 2, dist = "normal")
#'
#' # Same parameters with gamma distribution
#' f0(x = 0.5, mu = 12, sigma = 3, comp = 2, dist = "gamma")
#' }
f0 <- function(
    x,
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
    if (dist == "normal") {
      return((2/wind - 2*x/wind) * dhalfnorm(x, theta = sqrt(pi / 2) / (sqrt(2) * sigma)))
    } else {
      k <- (mu^2) / (sigma^2)
      theta <- (sigma^2) / mu
      if (k <= 0 || theta <= 0) return(0)
      bessel_val <- besselK(x / theta, 0.5 - k)
      bessel_val[!is.finite(bessel_val)] <- 0
      return_val <- (2/wind - 2*x/wind) *
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
  # Even comp = 2i  -> positive route i, mean = i*mu,    sd = sqrt(i)*sigma
  # Odd  comp = 2i+1 -> negative route i, mean = -i*mu,  sd = sqrt(i)*sigma
  if (comp %% 2 == 0) {
    i <- comp / 2
    route_mean <- i * mu
  } else {
    i <- (comp - 1) / 2
    route_mean <- -i * mu
  }
  route_sd <- sqrt(i) * sigma

  if (dist == "normal") {
    return((2/wind - 2*x/wind) * dnorm(x, mean = route_mean, sd = route_sd))
  } else {
    # Gamma: only even components (positive routes) are used
    k <- (mu^2) / (sigma^2)
    theta <- (sigma^2) / mu
    if (k <= 0 || theta <= 0) return(0)
    return((2/wind - 2*x/wind) * dgamma(x, shape = i * k, scale = theta))
  }
}

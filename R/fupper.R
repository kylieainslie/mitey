#' Calculate fupper for Different Components
#'
#' This function calculates the value of fupper based on the component.
#'
#' @param x The value at which to evaluate the function.
#' @param r The value of r.
#' @param mu The mean value.
#' @param sigma The standard deviation.
#' @param comp The component number (1 to 7).
#' @param dist string; assumed distribution of the serial interval; accepts "normal" or "gamma"; defaults to "normal".
#'
#' @return The calculated value of fupper.
#' @export
fupper <- function(
  x,
  r,
  mu,
  sigma,
  comp,
  dist = "normal"
) {
  # error messages
  if (dist != "normal" && dist != "gamma") {
    stop(
      "Incorrect distribution specified. Acceptable arguments are c('normal', 'gamma')"
    )
  }

  if (dist == "normal") {
    return(
      switch(
        comp,
        `1` = (r + 1 - x) *
          dhalfnorm(x, theta = sqrt(pi / 2) / (sqrt(2) * sigma)),
        `2` = (r + 1 - x) * dnorm(x, mean = mu, sd = sigma),
        `3` = (r + 1 - x) * dnorm(x, mean = -mu, sd = sigma),
        `4` = (r + 1 - x) * dnorm(x, mean = 2 * mu, sd = sqrt(2) * sigma),
        `5` = (r + 1 - x) * dnorm(x, mean = -2 * mu, sd = sqrt(2) * sigma),
        `6` = (r + 1 - x) * dnorm(x, mean = 3 * mu, sd = sqrt(3) * sigma),
        `7` = (r + 1 - x) * dnorm(x, mean = -3 * mu, sd = sqrt(3) * sigma)
      )
    )
  } else if (dist == "gamma") {
    # convert mean and sd for normal distn into shape and scale parameters for gamma distn.
    k <- (mu^2) / (sigma^2)
    theta <- (sigma^2) / mu

    if (k <= 0 || theta <= 0) {
      return(0)
    }

    return(
      switch(
        comp,
        `1` = {
          # Handle potential numerical issues
          bessel_val <- besselK(x / (theta), 0.5 - k)
          # Replace Inf with 0 while keeping finite values
          bessel_val[!is.finite(bessel_val)] <- 0
          return_val <- (r + 1 - x) *
            1 /
            sqrt(pi) *
            2^(3 / 2 - k) *
            theta^(-0.5 - k) *
            x^(-0.5 + k) *
            bessel_val *
            1 /
            gamma(k)
          # Replace NaN with 0 while keeping finite values
          return_val[is.nan(return_val)] <- 0
          return_val
        },
        `2` = (r + 1 - x) * dgamma(x, shape = k, scale = theta),
        `3` = (r + 1 - x) * dgamma(x, shape = k, scale = theta),
        `4` = (r + 1 - x) * dgamma(x, shape = 2 * k, scale = theta),
        `5` = (r + 1 - x) * dgamma(x, shape = 2 * k, scale = theta),
        `6` = (r + 1 - x) * dgamma(x, shape = 3 * k, scale = theta),
        `7` = (r + 1 - x) * dgamma(x, shape = 3 * k, scale = theta)
      )
    )
  }
}

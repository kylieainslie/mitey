#' Calculate f0 for Different Components
#'
#' This function calculates the value of f0 based on the component, where the components
#' represent the transmission routes: Co-Primary (CP), Primary-Secondary (PS), Primary-Tertiary (PT), and Primary-Quaternary (PQ). We split the PS, PT and PQ routes into two parts, such that
#'  - component 1: CP route
#'  - component 2+3: PS route
#'  - component 4+5: PT route
#'  - component 6+7: PQ route
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
#' @param comp integer; the component number (1 to 7).
#' @param dist string; assumed distribution of the serial interval; takes "normal" or "gamma"; defaults to "normal"
#'
#' @return The calculated value of f0.
#' @export
f0 <- function(x, mu, sigma, comp, dist = "normal") {

  # error messages
  if(dist != "normal" && dist != "gamma"){
    stop("Incorrect distribution specified. Acceptable arguments are c('normal', 'gamma')")
  }

  if (dist == "normal"){
    return(switch(comp,
           `1` = (2 - 2 * x) * dhalfnorm(x, theta = sqrt(pi / 2) / (sqrt(2) * sigma)),
           `2` = (2 - 2 * x) * dnorm(x, mean = mu, sd = sigma),
           `3` = (2 - 2 * x) * dnorm(x, mean = -mu, sd = sigma),
           `4` = (2 - 2 * x) * dnorm(x, mean = 2 * mu, sd = sqrt(2) * sigma),
           `5` = (2 - 2 * x) * dnorm(x, mean = -2 * mu, sd = sqrt(2) * sigma),
           `6` = (2 - 2 * x) * dnorm(x, mean = 3 * mu, sd = sqrt(3) * sigma),
           `7` = (2 - 2 * x) * dnorm(x, mean = -3 * mu, sd = sqrt(3) * sigma)
    ))
  } else if (dist == "gamma"){

    # convert mean and sd for normal distn into shape and scale parameters for gamma distn.
    k <- (mu^2) / (sigma^2)
    theta <- (sigma^2) / mu

    print(paste("k:", k))
    print(paste("theta:", theta))

    if (k <= 0 || theta <= 0) return(0)

    return(switch(comp,
           `1` = { # Handle potential numerical issues
             bessel_val <- besselK(x/(theta), 0.5 - k)
             print(paste("Bessel values:", bessel_val))
             # Replace Inf with 0 while keeping finite values
             bessel_val[!is.finite(bessel_val)] <- 0
             return_val <- (2 - 2 * x) * 1/sqrt(pi) * 2^(3/2 - k) * theta^(-0.5 - k) * x^(-0.5 + k) * bessel_val * 1/gamma(k)
             # Replace NaN with 0 while keeping finite values
             return_val[is.nan(return_val)] <- 0
             return_val
           },
           `2` = (2 - 2 * x) * dgamma(x, shape = k, scale = theta),
           `3` = (2 - 2 * x) * dgamma(x, shape = k, scale = theta),
           `4` = (2 - 2 * x) * dgamma(x, shape = 2 * k, scale = theta),
           `5` = (2 - 2 * x) * dgamma(x, shape = 2 * k, scale = theta),
           `6` = (2 - 2 * x) * dgamma(x, shape = 3 * k, scale = theta),
           `7` = (2 - 2 * x) * dgamma(x, shape = 3 * k, scale = theta)
    ))
  }
}

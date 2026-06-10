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
  dist = "normal"
) {
  # error messages
  if (!dist %in% c("normal", "gamma", "lognormal")) {
    stop(
      "Incorrect distribution specified. Acceptable arguments are c('normal', 'gamma', 'lognormal')"
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
  }  else if (dist == "lognormal") {
    if (mu <= 0 || sigma <= 0 || !is.finite(mu) || !is.finite(sigma)) {
      return(rep(0, length(x)))
    }
    
    sdlog <- sqrt(log(1 + sigma^2 / mu^2))
    meanlog <- log(mu) - 0.5 * sdlog^2
    
    fx <- function(z) {
      ifelse(z > 0, dlnorm(z, meanlog = meanlog, sdlog = sdlog), 0)
    }
    
    safe_integrate <- function(f, lower, upper) {
      val <- tryCatch(
        integrate(f, lower = lower, upper = upper)[[1]],
        error = function(e) 0
      )
      ifelse(is.finite(val), val, 0)
    }
    
    d_pt <- function(z) {
      if (!is.finite(z) || z <= 0) return(0)
      safe_integrate(
        f = function(t) fx(t) * fx(z - t),
        lower = 0,
        upper = z
      )
    }
    
    dens_one <- function(z) {
      if (!is.finite(z) || z <= 0) return(0)
      
      switch(
        as.character(comp),
        `1` = safe_integrate(
          f = function(t) 2 * fx(t) * fx(t + z),
          lower = 0,
          upper = Inf
        ),
        `2` = fx(z),
        `4` = d_pt(z),
        `6` = safe_integrate(
          f = function(u) d_pt(u) * fx(z - u),
          lower = 0,
          upper = z
        ),
        stop("For lognormal distribution, comp must be one of 1, 2, 4, or 6.")
      )
    }
    
    dens <- vapply(x, dens_one, numeric(1))
    
    return((r + 1 - x) * dens)
  }
}

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
        `1` = (x - r + 1) *
          dhalfnorm(x, theta = sqrt(pi / 2) / (sqrt(2) * sigma)),
        `2` = (x - r + 1) * dnorm(x, mean = mu, sd = sigma),
        `3` = (x - r + 1) * dnorm(x, mean = -mu, sd = sigma),
        `4` = (x - r + 1) * dnorm(x, mean = 2 * mu, sd = sqrt(2) * sigma),
        `5` = (x - r + 1) * dnorm(x, mean = -2 * mu, sd = sqrt(2) * sigma),
        `6` = (x - r + 1) * dnorm(x, mean = 3 * mu, sd = sqrt(3) * sigma),
        `7` = (x - r + 1) * dnorm(x, mean = -3 * mu, sd = sqrt(3) * sigma)
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
          return_val <- (x - r + 1) *
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
        `2` = (x - r + 1) * dgamma(x, shape = k, scale = theta),
        `3` = (x - r + 1) * dgamma(x, shape = k, scale = theta),
        `4` = (x - r + 1) * dgamma(x, shape = 2 * k, scale = theta),
        `5` = (x - r + 1) * dgamma(x, shape = 2 * k, scale = theta),
        `6` = (x - r + 1) * dgamma(x, shape = 3 * k, scale = theta),
        `7` = (x - r + 1) * dgamma(x, shape = 3 * k, scale = theta)
      )
    )
  } else if (dist == "lognormal") {
    if (mu <= 0 || sigma <= 0 || !is.finite(mu) || !is.finite(sigma)) {
      return(rep(0, length(x)))
    }
    
    # 1. Parameters for PS (comp 2), Lognormal dist
    sdlog <- sqrt(log(1 + sigma^2 / mu^2))
    meanlog <- log(mu) - 0.5 * sdlog^2
    
    # 2. FWA parameters for PT (comp 4) - Sum of 2 i.i.d lognormals
    # Moment Matching: E[S2] = 2*mu, Var(S2) = 2*sigma^2
    sdlog_pt <- sqrt(log(1 + (sigma^2) / (2 * mu^2)))
    meanlog_pt <- log(2 * mu) - 0.5 * sdlog_pt^2
    
    # 3. FWA parameters for PQ (comp 6) - Sum of 3 i.i.d lognormals
    sdlog_pq <- sqrt(log(1 + (sigma^2) / (3 * mu^2)))
    meanlog_pq <- log(3 * mu) - 0.5 * sdlog_pq^2
    
    # Base density function (needed only for CP integration now)
    fx <- function(z) {
      ifelse(z > 0, dlnorm(z, meanlog = meanlog, sdlog = sdlog), 0)
    }
    
    # Safe integrate wrapper (for CP path)
    safe_integrate <- function(f, lower, upper) {
      val <- tryCatch(
        integrate(f, lower = lower, upper = upper)[[1]],
        error = function(e) 0
      )
      ifelse(is.finite(val), val, 0)
    }
    
    # Evaluate density based on the transmission component
    if (comp == 1) {
      # CP path: Absolute difference, still requires numerical integration
      dens_one <- function(z) {
        if (!is.finite(z) || z <= 0) return(0)
        safe_integrate(
          f = function(t) 2 * fx(t) * fx(t + z),
          lower = 0,
          upper = Inf
        )
      }
      dens <- vapply(x, dens_one, numeric(1))
      
    } else if (comp == 2) {
      # PS path: Direct density
      dens <- ifelse(x > 0, dlnorm(x, meanlog = meanlog, sdlog = sdlog), 0)
      
    } else if (comp == 4) {
      # PT path: FWA approximation
      dens <- ifelse(x > 0, dlnorm(x, meanlog = meanlog_pt, sdlog = sdlog_pt), 0)
      
    } else if (comp == 6) {
      # PQ path: FWA approximation
      dens <- ifelse(x > 0, dlnorm(x, meanlog = meanlog_pq, sdlog = sdlog_pq), 0)
      
    } else {
      stop("For lognormal distribution, comp must be one of 1, 2, 4, or 6.")
    }
    
    # Return density multiplied by the flower triangular weight
    return((x - r + 1) * dens)
  }
  
}

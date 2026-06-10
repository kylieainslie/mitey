f_lnorm <- function(
    x,
    w1,
    w2,
    w3,
    mu,
    sigma
) {
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
    
    term1 <- w1 * safe_integrate(
      f = function(t) 2 * fx(t) * fx(t + z),
      lower = 0,
      upper = Inf
    )
    
    term2 <- w2 * fx(z)
    
    term3 <- w3 * d_pt(z)
    
    term4 <- (1 - w1 - w2 - w3) * safe_integrate(
      f = function(u) d_pt(u) * fx(z - u),
      lower = 0,
      upper = z
    )
    
    term1 + term2 + term3 + term4
  }
  
  vapply(x, dens_one, numeric(1))
}
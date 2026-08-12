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
  
  sdlog_pt <- sqrt(log(1 + sigma^2 / (2 * mu^2)))
  meanlog_pt <- log(2 * mu) - 0.5 * sdlog_pt^2
  
  sdlog_pq <- sqrt(log(1 + sigma^2 / (3 * mu^2)))
  meanlog_pq <- log(3 * mu) - 0.5 * sdlog_pq^2
  
  term1 <- w1 * ifelse(
    x >= 0,
    dhalfnorm(x, theta = sqrt(pi / 2) / (sqrt(2) * sigma)),
    0
  )
  term2 <- w2 * ifelse(x > 0, dlnorm(x, meanlog = meanlog, sdlog = sdlog), 0)
  term3 <- w3 * ifelse(x > 0, dlnorm(x, meanlog = meanlog_pt, sdlog = sdlog_pt), 0)
  term4 <- (1 - w1 - w2 - w3) *
    ifelse(x > 0, dlnorm(x, meanlog = meanlog_pq, sdlog = sdlog_pq), 0)
  
  term1 + term2 + term3 + term4
}
#' Helper function to plot the fit of the serial interval distribution assuming an underlying normal distribution
#'
#' @param x quantile
#' @param w1 probability weight of being a co-primary case
#' @param w2 probability weight of being a primary-secondary case
#' @param w3 probability weight of being a primary-tertiary case
#' @param mu mean serial interval
#' @param sigma standard deviation of serial interval
#' @returns the weighted density for each value of x
#' @export
#' @importFrom stats dnorm
f_norm <- function(x, w1, w2, w3, mu, sigma){

  term1 <- w1*dhalfnorm(x, sqrt(pi/2)/(sqrt(2)*sigma))
  term2 <- w2*dnorm(x,mean= mu, sd=sigma)
  term3 <- w3*dnorm(x,mean=2*mu,sd= sqrt(2)*sigma)
  term4 <- (1-w1-w2-w3)*dnorm(x, mean=3*mu, sd =sqrt(3)*sigma)

  return(term1 + term2 + term3 + term4)
}

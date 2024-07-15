#' Compute the likelihood of a transmission pair given the serial interval distribution
#'
#' Two distributions can be specified for the serial interval distribution: Normal and Gamma. A Normal distribution with mean \mu and standard deviation \sigma is assumed by default. If a Gamma distribution is specified, the mean \mu and standard deviation \sigma are converted to the shape and rate parameters for a Gamma distribution as:
#' \[
#' \beta = \frac{\mu}{\sigma^2}
#' \alpha = (\frac{\mu}{\beta})^2
#' \]
#' @param serial_interval numeric; time between infection or symptom onset of two individuals
#' @param mu numeric; mean of serial interval distribution
#' @param sigma numeric; standard deviation of serial interval distribution
#' @param distn string; distribution to be assumed for serial interval. Defaults to "normal", but also accepts "gamma".
#' @return value of the likelihood of the input serial_interval
#' @export
get_likelihood <- function(serial_interval, mu, sigma, distn = "normal") {
  if (distn == "normal"){
    likelihood <- dnorm(serial_interval, mean = mu, sd = sigma)
  } else if (distn == "gamma"){
    beta <- mu/sigma^2
    alpha <- (mu/beta)^2
    likelihood <- dgamma(serial_interval, shape = alpha, rate = beta)
  } else{
    stop("Incorrect serial distribution specified.
         Allowable inputs are distn = c('normal', 'gamma')")
  }

  return(likelihood)
}

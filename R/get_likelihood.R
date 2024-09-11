#' Compute the likelihood of a transmission pair given the serial interval distribution
#'
#' Two distributions can be specified for the serial interval distribution: Normal and Gamma.
#' A Normal distribution with mean \eqn{\mu} and standard deviation \eqn{\sigma} is assumed by default.
#' If a Gamma distribution is specified, the mean \eqn{\mu} and standard deviation \eqn{\sigma}
#' are converted to the shape and rate parameters for a Gamma distribution as:
#' \deqn{
#' \beta = \frac{\mu}{\sigma^2}
#' \alpha = \left(\frac{\mu}{\beta}\right)^2
#' }
#'
#' @param serial_intervals numeric vector of serial intervals
#' @param mu numeric; mean of serial interval distribution
#' @param sigma numeric; standard deviation of serial interval distribution
#' @param distn string; distribution to be assumed for serial interval. Defaults to "normal", but also accepts "gamma".
#' @return Likelihood of the transmission pair.
#' @export
#' @importFrom stats dgamma
#'
get_likelihood <- function(serial_intervals, mu, sigma, distn = "normal") {

  ### Checks
  # Ensure serial_intervals is numeric
  if(!is.numeric(serial_intervals)) {
    stop("serial_intervals must be numeric")
  }

  # Ensure no non-numeric values
  if(any(is.na(serial_intervals))) {
    stop("serial_intervals contains NA values")
  }

  if(any(is.nan(serial_intervals))) {
    stop("serial_intervals contains NaN values")
  }

  if(any(is.null(serial_intervals))) {
    stop("serial_intervals contains NULL values")
  }

  ### Compute likelihood
  if (distn == "normal"){
    likelihood <- dnorm(serial_intervals, mean = mu, sd = sigma)
  } else if (distn == "gamma"){
    beta <- mu/sigma^2
    alpha <- (mu/beta)^2
    likelihood <- dgamma(serial_intervals, shape = alpha, rate = beta)
  } else{
    stop("Incorrect serial distribution specified.
         Allowable inputs are distn = c('normal', 'gamma')")
  }

  return(likelihood)
}

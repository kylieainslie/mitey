#' Compute the likelihood of a transmission pair given the serial interval distribution
#'
#' Two distributions can be specified for the serial interval distribution: Normal and Gamma.
#' A Normal distribution with mean \code{\link{mu}} and standard deviation \code{\link{sigma}} is assumed by default.
#' If a Gamma distribution is specified, the mean \code{\link{mu}} and standard deviation \code{\link{sigma}}
#' are converted to the shape and rate parameters for a Gamma distribution as:
#' \deqn{
#' \beta = \frac{\code{\link{mu}}}{\code{\link{sigma}}^2}
#' \alpha = \left(\frac{\code{\link{mu}}}{\code{\link{beta}}}\right)^2
#' }
#'
#' @param serial_interval numeric; Time difference \( t_1 - t_2 \) (in days) between two events.
#' @param mu numeric; mean of serial interval distribution
#' @param sigma numeric; standard deviation of serial interval distribution
#' @param distn string; distribution to be assumed for serial interval. Defaults to "normal", but also accepts "gamma".
#' @param tail_cut numeric; number of days beyond which the likelihood of transmission between two events is 0. Defaults to NULL.
#' @param positive_only logical; if TRUE, only positive serial intervals are considered. If the serial interval is negative, the function returns 0.
#' @return Likelihood of the transmission pair.
#' @export
#' @importFrom stats dgamma
#' @examples
#' # Compute likelihood for a time difference of 5 days
#' get_likelihood(5, mu = 4, sigma = 2)
#'
get_likelihood <- function(serial_interval, mu, sigma, distn = "normal",
                           tail_cut = NULL, positive_only = TRUE) {

  # cut tail
  if(!is.null(tail_cut) & serial_interval > tail_cut){
    return(0)
  }

  # consider only positive serial intervals
  if(distn == "gamma"){positive_only = TRUE}

  if(positive_only & serial_interval < 0){
    return(0)
  }

  # get likelihood based on distribution
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

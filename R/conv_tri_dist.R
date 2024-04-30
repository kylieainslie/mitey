#' Convolution of the triangular distribution with the mixture component density (continuous case)
#'
#' We split the folded normal distribution for Primary-Secondary, Primary-Tertiary and Primary-Quaternary routes into two parts
#' + component 1: Co-Primary route
#' + component 2+3: Primary-Secondary route
#' + component 4+5: Primary-Tertiary route
#' + component 6+7: Primary-Quaternary route
#' @param x vector of index case to case intervals
#' @param sigma standard deviation of density distribution
#' @param r description??
#' @param mu mean of density distribution
#' @param route integer; between 1 and 7 and indicates the route of transmission.
#' @param quantity character; "mean", "lower", "upper"
#'
#' @return vector of density draws for each value of x
#' @export
#'
#' @importFrom fdrtool dhalfnorm
#' @importFrom stats dnorm
#'
#' @example
#' iccs <- 1:30
#' conv_tri_dist(x = iccs, sigma = 3, r = 10, mu = 15, route = 1)

conv_tri_dist <- function(x, sigma, r, mu, route, quantity = "mean"){

  # determine the distribution to draw from based on the route of transmission
  if(route == 1){
    dist_draw <- dhalfnorm(x, theta = sqrt(pi/2)/(sqrt(2)*sigma))
  } else if(route == 2){
    dist_draw <- dnorm(x, mean = mu, sd = sigma)
  } else if(route == 3){
    dist_draw <- dnorm(x, mean = -mu, sd = sigma)
  } else if (route == 4){
    dist_drawm <- dnorm(x, mean = 2*mu, sd = sqrt(2)*sigma)
  } else if (route == 5){
    dist_draw <- dnorm(x, mean = -2*mu, sd = sqrt(2)*sigma)
  } else if (route == 6){
    dist_draw <- dnorm(x, mean = 3*mu, sd = sqrt(3)*sigma)
  } else if (route == 7){
    dist_draw <- dnorm(x, mean = -3*mu, sd = sqrt(3)*sigma)
  } else {
    stop("Invalid route argument. route must be an integer from 1 to 7 (inclusive).")
  }

  # calculate the quantity of interest
  if(quantity == "mean"){
    rtn <- (2 - 2 * x) * dist_draw
  } else if (quantity == "lower"){
    rtn <- (x - r + 1) * dist_draw
  } else if (quantity == "upper"){
    rtn <- (r + 1 - x) * dist_draw
  } else {
    stop("Invalid string in quantity argument. Valid inputs are 'mean', 'lower', and 'upper'.")
  }

  # output
  return(rtn)
}

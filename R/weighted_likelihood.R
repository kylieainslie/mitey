#' weighted likelihood
#' each point adds to likelihood given weight belonging to component 2: w2
#'
#' @param dat data to be used for maximisation
#' @param par parameters to be maximised
#' @param tau2 numeric; mixture weights for (primary, secondary) pairs
#' @return value of the negative sum of the weighted log likelihood
#' @export
wt_loglik <- function(dat, par, tau2){
  som <- 0

  # convert par values to be appropriate for gamma dist
  k <- (par[1]^2) / (par[2]^2)
  theta <- (par[2]^2) / par[1]

  if(par[1] < 0 | par[2] < 0){
    som <- -1000000
  }
  else{
    som <- som + sum(log(tau2*dgamma(dat, shape = k, scale = theta)))
  }
  return(-som)
}

#' weighted likelihood
#' each point adds to likelihood given weight belonging to component 2: w2
#'
#' @param par parameters to be maximised
#' @param dat data to be used for maximisation
#' @param tau2 numeric; mixture weights for (primary, secondary) pairs
#' @return value of the negative sum of the weighted log likelihood
#' @export
wt_loglik <- function(par, dat, tau2){
  som <- 0

  dat <- ifelse(dat == 0, 0.00001, dat)

  # convert par values to be appropriate for gamma dist
  k <- (par[1]^2) / (par[2]^2)
  theta <- (par[2]^2) / par[1]

  if(k < 0 | theta < 0){
    som <- -1000000
  } else{

    # Calculate gamma density and add a small constant to avoid log(0)
    epsilon <- 1e-10
    densities <- tau2 * dgamma(dat, shape = k, scale = theta)
    densities <- pmax(densities, epsilon)  # Replace zeros with a small constant

    # Compute the negative log-likelihood
    som <- sum(log(densities))
    #som <- som + sum(log(tau2*dgamma(dat, shape = k, scale = theta)))
  }
  return(-som)
}

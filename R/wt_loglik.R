#' weighted likelihood for optimising parameters assuming underlying gamma distribution
#' each point adds to likelihood given weight belonging to component 2: w2
#'
#' @param par parameters to be maximised
#' @param dat data to be used for maximisation
#' @param tau2 numeric; mixture weights for (primary, secondary) pairs
#' @return value of the negative sum of the weighted log likelihood
#' @export
wt_loglik <- function(
  par,
  dat,
  tau2
) {
  som <- 0

  dat <- ifelse(dat == 0, 0.00001, dat)

  # convert par values to be appropriate for gamma dist
  k <- (par[1]^2) / (par[2]^2)
  theta <- (par[2]^2) / par[1]

  # Check if k and theta are valid
  if (k <= 0 || theta <= 0 || !is.finite(k) || !is.finite(theta)) {
    # Apply a large negative penalty if k or theta are invalid
    return(1e10)
  }

  # Calculate the log-likelihood, handling non-finite values
  log_likelihood <- sum(log(tau2 * dgamma(dat, shape = k, scale = theta)))

  if (!is.finite(log_likelihood)) {
    return(1e10) # Apply a large penalty if log-likelihood is non-finite
  }

  return(-log_likelihood) # Return negative log-likelihood for minimization
  #som <- som + sum(log(tau2*dgamma(dat, shape = k, scale = theta)))
}

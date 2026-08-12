#' Calculate Weighted Negative Log-Likelihood for Lognormal PS Component
#'
#' Computes the weighted negative log-likelihood using numerical integration
#' for the primary-secondary (PS) component of the lognormal distribution.
#'
#' @param par numeric vector; c(mu, sigma) representing the mean and sd on the natural scale.
#' @param dat numeric vector; the ICC intervals.
#' @param tau2 numeric vector; posterior weights for the PS component from the E-step.
#' @param zero_idx integer vector; indices of observations that were originally zero.
#' @return numeric; negative log-likelihood value for minimization.
#' @keywords internal
lognormal_ps_loglik <- function(par, dat, tau2, zero_idx = integer(0)) {
  mu_opt <- par[1]
  sig_opt <- par[2]
  
  # Apply a severe penalty if parameters step out of valid  bounds
  if (mu_opt <= 0 || sig_opt <= 0 || !is.finite(mu_opt) || !is.finite(sig_opt)) {
    return(1e10)
  }
  # Initialize loglikelihood
  ll <- 0
  zero_idx <- union(zero_idx, which(dat == 0))
  dat <- ifelse(dat == 0, 1, dat)

  for (i in seq_along(dat)) {
    is_zero_interval <- i %in% zero_idx

    # Calculate the integral probability for the PS component (comp = 2)
    prob <- integrate_component(
      d = if (is_zero_interval) 0 else dat[i],
      mu = mu_opt,
      sigma = sig_opt,
      comp = 2,
      dist = "lognormal",
      lower = !is_zero_interval
    )
    
    # Accumulate the weighted log-likelihood, preventing log(0)
    if (prob > 0) {
      ll <- ll + tau2[i] * log(prob)
    } else {
      ll <- ll + tau2[i] * log(.Machine$double.xmin)
    }
  }
  
  return(-ll)
}

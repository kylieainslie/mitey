#' Calculate Weighted Negative Log-Likelihood for Gamma Distribution Parameters
#'
#' Computes the weighted negative log-likelihood for gamma distribution parameters
#' in the M-step of the EM algorithm for serial interval estimation. This function
#' is used as the objective function for numerical optimization when the serial
#' interval distribution is assumed to follow a gamma distribution.
#'
#' The function converts mean and standard deviation parameters to gamma distribution
#' shape and scale parameters, then calculates the weighted log-likelihood based on
#' the posterior probabilities from the E-step. Returns the negative log-likelihood
#' for minimization by \code{\link{optim}}.
#'
#' @param par numeric vector of length 2; parameters to optimize where \code{par[1]}
#'            is the mean and \code{par[2]} is the standard deviation of the serial
#'            interval distribution
#' @param dat numeric vector; index case-to-case (ICC) intervals. Zero values are
#'            replaced with 0.00001 to avoid gamma distribution issues at zero
#' @param tau2 numeric vector; posterior probabilities (weights) that each observation
#'             belongs to the primary-secondary transmission component. These are
#'             typically derived from the E-step of the EM algorithm
#'
#' @return numeric; negative log-likelihood value for minimization. Returns a large
#'         penalty value (1e10) if parameters result in invalid gamma distribution
#'         parameters (non-positive shape/scale) or non-finite likelihood values
#'
#' @details
#' This function is used internally by \code{\link{si_estim}} when \code{dist = "gamma"}.
#' The gamma distribution is parameterized using shape (k) and scale (theta) parameters
#' derived from the mean and standard deviation:
#' \itemize{
#'   \item Shape: \eqn{k = \mu^2 / \sigma^2}
#'   \item Scale: \eqn{\theta = \sigma^2 / \mu}
#' }
#'
#' The weighted log-likelihood is calculated as:
#' \deqn{\sum_{i} \tau_{2,i} \log f(x_i | k, \theta)}
#' where \eqn{f(x_i | k, \theta)} is the gamma probability density function and
#' \eqn{\tau_{2,i}} are the weights from the E-step.
#'
#' @note This function is primarily intended for internal use within the EM algorithm.
#'       Users typically won't call this function directly but rather use
#'       \code{\link{si_estim}} with \code{dist = "gamma"}.
#'
#' @seealso \code{\link{si_estim}} for the main serial interval estimation function,
#'          \code{\link{optim}} for the optimization routine that uses this function
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' set.seed(123)
#' icc_intervals <- rgamma(50, shape = 2, scale = 3)
#' weights <- runif(50, 0.1, 1)
#' wt_loglik(c(5, 4), icc_intervals, weights)
#' }
#'
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
  log_likelihood <- sum(tau2 * log(dgamma(dat, shape = k, scale = theta)))

  if (!is.finite(log_likelihood)) {
    return(1e10) # Apply a large penalty if log-likelihood is non-finite
  }

  return(-log_likelihood) # Return negative log-likelihood for minimization
  #som <- som + sum(log(tau2*dgamma(dat, shape = k, scale = theta)))
}

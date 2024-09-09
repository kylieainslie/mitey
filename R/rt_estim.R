#' Estimate Time-varying Reproduction Number
#'
#' This function estimates the time-varying reproduction number (R_t) using the method proposed by Wallinga and Teunis (AJE 2004).
#' The function takes as input a data frame with incidence data and requires the specification of the serial interval distribution.
#' The user must input the mean and standard deviation of the serial interval distribution and select the functional form: Normal or Gamma.
#'
#' @param inc_dat data frame; data frame with incidence data. The data frame should have two columns: inc (daily incidence) and date.
#' @param mean_si numeric; mean of serial interval distribution
#' @param sd_si numeric; standard deviation of serial interval distribution
#' @param dist_si string; distribution to be assumed for serial interval. Accepts "normal" or "gamma".
#' @param cut_tail numeric; number of days beyond which the likelihood of transmission between two events is 0. Defaults to NULL.
#' @param pos_only logical; if TRUE, only positive serial intervals are considered. If the serial interval is negative, the function returns 0.
#' @param n_boostrap integer; number of bootstrap samples of the serial interval distribution
#' @return A numeric vector with the expected reproduction number for each day.
#' @export
#' @importFrom stats quantile
#' @importFrom stats rnorm
#' @importFrom stats rgamma
rt_estim <- function(inc_dat, mean_si, sd_si, dist_si = "normal",
                     cut_tail = NULL, pos_only = TRUE, n_bootstrap = 100){

# Pre-compute the likelihood values for all possible serial intervals
  serial_intervals <- outer(inc_dat$date, inc_dat$date, "-")

# Initialize a matrix to store R_t estimates from each bootstrap sample
  bootstrap_rt <- matrix(0, nrow = n_bootstrap, ncol = nrow(inc_dat))

  for (b in 1:n_bootstrap) {

    # Generate a perturbed serial interval distribution
    if (dist_si == "normal") {
      si_distribution <- rnorm(1000, mean = mean_si, sd = sd_si)
    } else if (dist_si == "gamma") {
      beta <- mean_si / sd_si^2
      alpha <- (mean_si / beta)^2
      si_distribution <- rgamma(1000, shape = alpha, rate = beta)
    }

    # get mean of perturbed serial interval distribution
    new_mean_si <- mean(si_distribution)

    # Initialize a matrix to store likelihood values
    likelihood_values <- matrix(NA, nrow = nrow(serial_intervals),
                                ncol = ncol(serial_intervals))

    # Loop over all pairs of serial intervals
    likelihood_values <- apply(serial_intervals, 1:2, get_likelihood,
                               mu = new_mean_si, sigma = sd_si, distn = dist_si,
                               tail_cut = 180, positive_only = TRUE)

  # Calculation the likelihood matrix
  nt <- nrow(inc_dat)
  # Initialize the likelihood matrix
  likelihood_mat <- matrix(0, nrow = nt, ncol = nt)

  for (i in 1:nt) {
    for (j in 1:nt) {
      if (i > j) {  # Only consider pairs where i > j
        likelihood_mat[i, j] <- likelihood_values[i, j] * inc_dat$inc[j]
      }
    }
  }

  # Get the marginal likelihood by summing each row
  marginal_likelihood <- rowSums(likelihood_mat)

  # Calculate the probability matrix
  prob_mat <- sweep(likelihood_mat, 1, marginal_likelihood, FUN = "/")
  prob_mat[is.nan(prob_mat)] <- NA

  # Calculate the expected reproduction number per day (R_t)
  expected_rt <- colSums(prob_mat, na.rm = TRUE)

  # Store the bootstrap result
  bootstrap_rt[b, ] <- expected_rt
  }

  # Calculate the mean R_t and confidence intervals across bootstrap samples
  rt_mean <- apply(bootstrap_rt, 2, mean)
  rt_lower <- apply(bootstrap_rt, 2, quantile, probs = 0.025)
  rt_upper <- apply(bootstrap_rt, 2, quantile, probs = 0.975)

  # Return a list containing the mean R_t, and its lower and upper bounds
  return(list(R_t_mean = rt_mean, R_t_lower = rt_lower, R_t_upper = rt_upper))

}

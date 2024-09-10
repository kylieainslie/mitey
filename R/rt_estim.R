#' Estimate Time-varying Reproduction Number
#'
#' This function estimates the time-varying reproduction number (R_t) using the method proposed by Wallinga and Teunis (AJE 2004).
#' The function takes as input a data frame with incidence data and requires the specification of the serial interval distribution.
#' The user must input the mean and standard deviation of the serial interval distribution and select the functional form: Normal or Gamma.
#'
#' @param inc_dat data frame; data frame with incidence data. The data frame should have two columns: inc (daily incidence) and onset_date (date of onset of symptoms).
#' @param mean_si numeric; mean of serial interval distribution
#' @param sd_si numeric; standard deviation of serial interval distribution
#' @param dist_si string; distribution to be assumed for serial interval. Accepts "normal" or "gamma".
#' @param cut_tail numeric; number of days beyond which the likelihood of transmission between two events is 0. Defaults to NULL.
#' @param pos_only logical; if TRUE, only positive serial intervals are considered. If the serial interval is negative, the function returns 0.
#' @param perturb_si_dist logical; if true a perturbed serial interval is generated based on dist_si, mean_si, and sd_si, then the mean and sd of the perturbed distribution are used in the calculation of the likelihood function.
#' @return A numeric vector with the expected reproduction number for each day.
#' @export
#' @import tidyr
#' @importFrom dplyr arrange
#' @importFrom stats rnorm
#' @importFrom stats rgamma
rt_estim <- function(inc_dat, mean_si, sd_si, dist_si = "normal",
                     cut_tail = NULL, pos_only = TRUE, perturb_si_dist = FALSE){

  # fill in missing dates, set inc = 0
  all_dates <- seq(min(inc_dat$onset_date), max(inc_dat$onset_date), by = "day")

  inc_dat <- inc_dat %>%
    complete(onset_date = all_dates,
             fill = list(inc = 0)
             ) %>%
    arrange(onset_date)

  nt <- nrow(inc_dat)

  # Pre-compute the likelihood values for all possible serial intervals
  serial_intervals <- outer(inc_dat$onset_date, inc_dat$onset_date, "-")

    # Generate a perturbed serial interval distribution
    if (perturb_si_dist){
      si_distribution <- switch(dist_si,
                                "normal" = rnorm(1000, mean = mean_si, sd = sd_si),
                                "gamma" = {
                                  beta <- mean_si / sd_si^2
                                  alpha <- (mean_si / beta)^2
                                  rgamma(1000, shape = alpha, rate = beta)
                                },
                                stop("Unsupported distribution type"))
      # get mean of perturbed serial interval distribution
      mean_si <- mean(si_distribution)
      sd_si <- sd(si_distribution)
    }

  # Loop over all pairs of serial intervals
  likelihood_values <- apply(serial_intervals, 1:2, get_likelihood,
                             mu = mean_si, sigma = sd_si, distn = dist_si,
                             tail_cut = cut_tail, positive_only = pos_only)

  # Calculation the likelihood matrix
  likelihood_mat <- likelihood_values * outer(rep(1, nt), inc_dat$inc)
  # Initialize the likelihood matrix
  # likelihood_mat <- matrix(0, nrow = nt, ncol = nt)
  #
  # for (i in 1:nt) {
  #   for (j in 1:nt) {
  #     if (i > j) {  # Only consider pairs where i > j
  #       likelihood_mat[i, j] <- likelihood_values[i, j] * inc_dat$inc[j]
  #     }
  #   }
  # }

  # Zero out upper triangle and diagonal
  likelihood_mat[upper.tri(likelihood_mat, diag = TRUE)] <- 0

  # Get the marginal likelihood by summing each row
  marginal_likelihood <- rowSums(likelihood_mat)

  # Calculate the probability matrix
  prob_mat <- sweep(likelihood_mat, 1, marginal_likelihood, FUN = "/")
  prob_mat[is.nan(prob_mat)] <- 0

  # Calculate the expected reproduction number per day (R_t)
  expected_rt <- colSums(prob_mat, na.rm = TRUE)

  # correct for right-truncation using the method of Cauchemez et al. 2006
  # Apply the right truncation correction to each observation
  onset_times <- 1:nt
  correction_factors <- sapply(onset_times, function(t) {
    right_truncation_correction(t, nt, mean_si, sd_si)
  })

  # Adjusted Rt values
  adjusted_rt <- expected_rt * correction_factors

  # Return a list containing the mean R_t, and its lower and upper bounds
  return(list(rt = expected_rt, rt_adjusted = adjusted_rt))

}

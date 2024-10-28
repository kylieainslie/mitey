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
#' @param perturb_si_dist logical; if true a perturbed serial interval is generated based on dist_si, mean_si, and sd_si, then the mean and sd of the perturbed distribution are used in the calculation of the likelihood function.
#' @return A numeric vector with the expected reproduction number for each day.
#' @export
#' @import tidyr
#' @importFrom dplyr arrange
#' @importFrom stats rnorm
#' @importFrom stats rgamma
#' @importFrom lubridate is.Date
rt_estim <- function(inc_dat, mean_si, sd_si, dist_si = "normal",
                     perturb_si_dist = FALSE){

  # fill in missing dates, set inc = 0
  if(is.Date(inc_date$onset_date)){
    all_dates <- seq(min(inc_dat$onset_date), max(inc_dat$onset_date), by = "day")
  } else{
    all_dates <- seq(min(inc_dat$onset_date), max(inc_dat$onset_date), by = 1)
  }
  inc_dat <- inc_dat %>%
    complete(onset_date = all_dates,
             fill = list(inc = 0)
             ) %>%
    arrange(.data$onset_date)

  nt <- nrow(inc_dat)

  # Pre-compute the likelihood values for all possible serial intervals (t-u)
  serial_intervals <- as.vector(outer(inc_dat$onset_date, inc_dat$onset_date, "-"),
                                mode = "numeric")

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
  # Rt = sum{b(t)g(t-u)/sum(b(t-a)g(a))},
  # where b(t) is the incidence on day t and g(t) is the probability density of
  # the serial interval distribution

  # Get likelihood for each serial interval based on probability density (g(t-u))
  likelihood_values <- get_likelihood(serial_intervals, mu = mean_si,
                                      sigma = sd_si, distn = dist_si)

  # convert vector of likelihood values into a nt x nt matrix
  likelihood_values_mat <- matrix(likelihood_values, nrow = nt, ncol = nt)

  # b(u)g(t-u): multiply columns of likelihood matrix by the incidence
  likelihood_mat <- likelihood_values_mat * outer(inc_dat$inc, rep(1, nt))

    # Zero out upper triangle and diagonal
    likelihood_mat[upper.tri(likelihood_mat, diag = TRUE)] <- 0

  # sum(b(t-a)g(a)): calculate marginal likelihood by summing over rows
  marginal_likelihood <- rowSums(likelihood_mat)

  # b(u)g(t-u)/sum(b(t-a)g(a)): divide likelihood matrix columns by marginal likelihood
  prob_mat <- sweep(likelihood_mat, 2, marginal_likelihood, FUN = "/")
  prob_mat[is.nan(prob_mat)] <- 0

  # b(t) * [b(u)g(t-u)/sum(b(t-a)g(a))]: multiply by incidence
  prob_mat2 <- prob_mat * outer(rep(1, nt), inc_dat$inc)

  # sum_t{b(t) * [b(u)g(t-u)/sum(b(t-a)g(a))]} * 1/(b(u)): sum over columns and divide by incidence
  expected_rt <- colSums(prob_mat2, na.rm = TRUE) / inc_dat$inc

  # correct for right-truncation using the method of Cauchemez et al. 2006
  # Apply the right truncation correction to each observation
  onset_times <- 1:nt
  correction_factors <- sapply(onset_times, function(t) {
    right_truncation_correction(t, nt, mean_si, sd_si)
  })

  # Adjusted Rt values
  adjusted_rt <- expected_rt * correction_factors

  # Return a list containing the mean R_t, and its lower and upper bounds
  rtn <- data.frame(onset_date = all_dates, rt = expected_rt, rt_adjusted = adjusted_rt)
  return(rtn)

}

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
#' @return A numeric vector with the expected reproduction number for each day.
#' @export
rt_estim <- function(inc_dat, mean_si, sd_si, dist_si = c("normal", "gamma"),
                     cut_tail = NULL, pos_only = TRUE){

# Pre-compute the likelihood values for all possible serial intervals
serial_intervals <- outer(inc_dat$date, inc_dat$date, "-")

# Apply get_likelihood only to positive serial_intervals
likelihood_values <- apply(serial_intervals, 1:2, get_likelihood, mu = mean_si,
                           sigma = sd_si, distn = "normal", tail_cut = 180,
                           positive_only = TRUE)

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

return(expected_rt)

}

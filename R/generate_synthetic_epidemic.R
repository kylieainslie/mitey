#' Generate synthetic incidence data from a known reproduction number
#'
#' @param true_r Vector of true reproduction numbers
#' @param si_mean Mean of serial interval
#' @param si_sd SD of serial interval
#' @param si_dist Serial interval distribution ("gamma" or "normal")
#' @param initial_cases Number of initial cases
#' @return Data frame with dates, true R, and incidence
generate_synthetic_epidemic <- function(true_r, si_mean, si_sd,
                                        si_dist = "gamma", initial_cases = 10) {
  n_days <- length(true_r)
  dates <- seq(as.Date("2023-01-01"), by = "day", length.out = n_days)

  # Generate serial interval distribution
  if (si_dist == "gamma") {
    shape <- (si_mean / si_sd)^2
    rate <- si_mean / (si_sd^2)
    max_si <- qgamma(0.99, shape, rate)
    si_prob <- dgamma(1:ceiling(max_si), shape, rate)
  } else {
    max_si <- qnorm(0.99, si_mean, si_sd)
    si_prob <- dnorm(1:ceiling(max_si), si_mean, si_sd)
    si_prob <- si_prob / sum(si_prob)  # Normalize to ensure sum equals 1
  }

  # Initialize incidence
  incidence <- rep(0, n_days)
  incidence[1] <- initial_cases

  # Generate cases using the renewal equation directly
  for (t in 2:n_days) {
    # Sum over all previous days, weighted by serial interval
    lambda_t <- 0
    for (s in 1:(t-1)) {
      interval <- t - s
      if (interval <= length(si_prob)) {
        lambda_t <- lambda_t + incidence[s] * true_r[s] * si_prob[interval]
      }
    }
    # Generate new cases from Poisson distribution with expected value lambda_t
    incidence[t] <- rpois(1, lambda_t)
  }

  # Return data frame
  return(data.frame(date = dates, true_r = true_r, incidence = incidence))
}

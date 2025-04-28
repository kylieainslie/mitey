#' Estimate time-varying reproduction number using Wallinga-Lipsitch method with bootstrap confidence intervals
#'
#' @param incidence Numeric vector of daily case counts
#' @param dates Vector of dates corresponding to the incidence data
#' @param si_mean Mean of the serial interval distribution
#' @param si_sd Standard deviation of the serial interval distribution
#' @param si_dist Distribution to use for serial interval ("gamma" or "normal")
#' @param smoothing Window size for smoothing estimates (0 for no smoothing)
#' @param bootstrap Logical; whether to compute bootstrap confidence intervals
#' @param n_bootstrap Number of bootstrap samples to generate
#' @param conf_level Confidence level for intervals (0.95 = 95% CI)
#' @return Data frame with dates, R estimates, and confidence intervals if requested
#'
wallinga_lipsitch <- function(incidence, dates, si_mean, si_sd,
                              si_dist = "gamma", smoothing = 7,
                              bootstrap = FALSE, n_bootstrap = 1000,
                              conf_level = 0.95) {

  # Input validation
  stopifnot(length(incidence) == length(dates))
  stopifnot(is.numeric(incidence))
  stopifnot(all(incidence >= 0))
  stopifnot(si_mean > 0)
  stopifnot(si_sd > 0)
  stopifnot(si_dist %in% c("gamma", "normal"))

  n <- length(incidence)

  # Create matrix of all pairwise day differences
  day_diffs <- matrix(NA, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      day_diffs[i, j] <- as.numeric(difftime(dates[i], dates[j], units = "days"))
    }
  }

  # Calculate probability density of the serial interval for each day difference
  si_prob <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      # Only consider forward transmission (positive day differences)
      if (day_diffs[i, j] > 0) {
        if (si_dist == "gamma") {
          # Convert mean and SD to shape and rate parameters
          shape <- (si_mean / si_sd)^2
          rate <- si_mean / (si_sd^2)
          si_prob[i, j] <- dgamma(day_diffs[i, j], shape = shape, rate = rate)
        } else if (si_dist == "normal") {
          si_prob[i, j] <- dnorm(day_diffs[i, j], mean = si_mean, sd = si_sd)
          # Truncate negative values for normal distribution
          if (day_diffs[i, j] <= 0) {
            si_prob[i, j] <- 0
          }
        }
      }
    }
  }

  # Function to calculate R from a given incidence vector
  calculate_r <- function(inc) {
    # Calculate the relative likelihood that case i was infected by case j
    rel_likelihood <- matrix(0, nrow = n, ncol = n)
    for (i in 1:n) {
      denominator <- sum(inc * si_prob[i, ])
      if (denominator > 0) {
        rel_likelihood[i, ] <- inc * si_prob[i, ] / denominator
      }
    }

    # Calculate expected number of secondary cases
    r_estimate <- rep(NA, n)
    for (j in 1:n) {
      if (inc[j] > 0) {
        r_estimate[j] <- sum(inc * rel_likelihood[, j]) / inc[j]
      }
    }

    # Apply smoothing if requested
    if (smoothing > 1) {
      r_smooth <- rep(NA, n)
      half_window <- floor(smoothing / 2)

      for (i in 1:n) {
        start_idx <- max(1, i - half_window)
        end_idx <- min(n, i + half_window)
        valid_values <- r_estimate[start_idx:end_idx]
        valid_values <- valid_values[!is.na(valid_values) & is.finite(valid_values)]

        if (length(valid_values) > 0) {
          r_smooth[i] <- mean(valid_values)
        }
      }

      r_estimate <- r_smooth
    }

    # Add right-truncation correction
    correction <- rep(1, n)
    for (i in 1:n) {
      days_since_case <- as.numeric(max(dates) - dates[i])
      if (si_dist == "gamma") {
        shape <- (si_mean / si_sd)^2
        rate <- si_mean / (si_sd^2)
        prob_observed <- pgamma(days_since_case, shape = shape, rate = rate)
      } else {
        prob_observed <- pnorm(days_since_case, mean = si_mean, sd = si_sd)
      }

      if (prob_observed > 0.5) {
        correction[i] <- 1 / prob_observed
      } else {
        correction[i] <- NA
      }
    }

    r_corrected <- r_estimate * correction

    # Cap extremely high values
    r_corrected[r_corrected > 10] <- NA

    return(list(r = r_estimate, r_corrected = r_corrected))
  }

  # Calculate point estimates
  point_estimates <- calculate_r(incidence)

  # Initialize result data frame
  result <- data.frame(
    date = dates,
    incidence = incidence,
    R = point_estimates$r,
    R_corrected = point_estimates$r_corrected
  )

  # Add bootstrap confidence intervals if requested
  if (bootstrap) {
    # Initialize matrices to store bootstrap results
    bootstrap_r <- matrix(NA, nrow = n_bootstrap, ncol = n)
    bootstrap_r_corrected <- matrix(NA, nrow = n_bootstrap, ncol = n)

    generate_case_bootstrap <- function(incidence) {
      # Create a case-level representation (each case is an individual unit)
      case_days <- rep(1:length(incidence), times = incidence)

      # Resample individual cases
      n_cases <- sum(incidence)
      bootstrap_cases <- sample(case_days, n_cases, replace = TRUE)

      # Convert back to daily incidence
      bootstrap_incidence <- tabulate(bootstrap_cases, nbins = length(incidence))

      return(bootstrap_incidence)
    }

    # Generate bootstrap samples
    for (i in 1:n_bootstrap) {
      bootstrap_incidence <- generate_case_bootstrap(incidence)
      bootstrap_estimates <- calculate_r(bootstrap_incidence)
      bootstrap_r[i, ] <- bootstrap_estimates$r
      bootstrap_r_corrected[i, ] <- bootstrap_estimates$r_corrected
    }

    # Calculate confidence intervals
    alpha <- (1 - conf_level) / 2
    quantiles <- c(alpha, 1 - alpha)

    # For each day, calculate quantiles of bootstrap distributions
    r_lower <- apply(bootstrap_r, 2, function(x) quantile(x, probs = quantiles[1], na.rm = TRUE))
    r_upper <- apply(bootstrap_r, 2, function(x) quantile(x, probs = quantiles[2], na.rm = TRUE))

    r_corrected_lower <- apply(bootstrap_r_corrected, 2, function(x) quantile(x, probs = quantiles[1], na.rm = TRUE))
    r_corrected_upper <- apply(bootstrap_r_corrected, 2, function(x) quantile(x, probs = quantiles[2], na.rm = TRUE))

    # Add to result data frame
    result$R_lower <- r_lower
    result$R_upper <- r_upper
    result$R_corrected_lower <- r_corrected_lower
    result$R_corrected_upper <- r_corrected_upper
  }

  return(result)
}

#' Calculate Serial Interval Probability Matrix
#'
#' Creates a matrix of probabilities based on the serial interval distribution
#' for all pairwise day differences.
#'
#' @param day_diffs Matrix of day differences between each pair of cases
#' @param si_mean Mean of the serial interval distribution
#' @param si_sd Standard deviation of the serial interval distribution
#' @param si_dist Distribution type ("gamma" or "normal")
#' @return Matrix of serial interval probabilities
#'
calculate_si_probability_matrix <- function(day_diffs, si_mean, si_sd, si_dist) {
  n <- nrow(day_diffs)
  si_prob <- matrix(0, nrow = n, ncol = n)

  for (i in 1:n) {
    for (j in 1:n) {
      ## Only consider cases where j could infect i (i.e., j comes before i)
      # This means day_diffs[i, j] > 0 (date i - date j > 0)
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

  return(si_prob)
}

#' Create Day Difference Matrix
#'
#' Creates a matrix of day differences between all pairs of cases.
#'
#' @param dates Vector of dates for each case
#' @return Matrix of day differences
#'
create_day_diff_matrix <- function(dates) {
  n <- length(dates)
  day_diffs <- matrix(NA, nrow = n, ncol = n)

  for (i in 1:n) {
    for (j in 1:n) {
      day_diffs[i, j] <- as.numeric(difftime(dates[i], dates[j], units = "days"))
    }
  }

  return(day_diffs)
}

#' Apply Smoothing to Estimates
#'
#' Applies moving average smoothing to R estimates.
#'
#' @param r_estimate Vector of R estimates
#' @param window Smoothing window size
#' @return Smoothed R estimates
#'
smooth_estimates <- function(r_estimate, window) {
  n <- length(r_estimate)
  r_smooth <- rep(NA, n)
  #half_window <- floor(window / 2)

  for (i in 1:n) {
    start_idx <- max(1, i - (window - 1))
    end_idx <- i
    valid_values <- r_estimate[start_idx:end_idx]
    valid_values <- valid_values[!is.na(valid_values) & is.finite(valid_values)]

    if (length(valid_values) > 0) {
      r_smooth[i] <- mean(valid_values)
    }
  }

  return(r_smooth)
}

#' Calculate Right Truncation Correction
#'
#' Calculates correction factors for right truncation of the time series.
#'
#' @param dates Vector of dates
#' @param si_mean Mean of the serial interval
#' @param si_sd Standard deviation of the serial interval
#' @param si_dist Distribution type ("gamma" or "normal")
#' @return Vector of correction factors
#'
calculate_truncation_correction <- function(dates, si_mean, si_sd, si_dist) {
  n <- length(dates)
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

  return(correction)
}

#' Generate Case Bootstrap Sample
#'
#' Generates a bootstrap sample of the incidence data by resampling individual cases.
#'
#' @param incidence Vector of case counts
#' @return Bootstrapped incidence vector
#'
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

#' Calculate R From Incidence and Serial Interval
#'
#' Core function to calculate R estimates from incidence data and serial interval probabilities.
#'
#' @param incidence Vector of case counts
#' @param si_prob Matrix of serial interval probabilities
#' @param dates Vector of dates
#' @param si_mean Mean of the serial interval
#' @param si_sd Standard deviation of the serial interval
#' @param si_dist Distribution type
#' @param smoothing Smoothing window size
#' @return List with R and R_corrected vectors
#'
calculate_r_estimates <- function(incidence, si_prob, dates, si_mean, si_sd, si_dist, smoothing) {
  n <- length(incidence)

  # Calculate the relative likelihood that case i was infected by case j
  rel_likelihood <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    denominator <- sum(incidence * si_prob[i, ])
    if (denominator > 0) {
      rel_likelihood[i, ] <- incidence * si_prob[i, ] / denominator
    }
  }

  # Calculate expected number of secondary cases
  r_estimate <- rep(NA, n)
  for (j in 1:n) {
    if (incidence[j] > 0) {
      r_estimate[j] <- sum(incidence * rel_likelihood[, j]) / incidence[j]
    }
  }

  # Apply smoothing if requested
  if (smoothing > 1) {
    r_estimate <- smooth_estimates(r_estimate, smoothing)
  }

  # Add right-truncation correction
  correction <- calculate_truncation_correction(dates, si_mean, si_sd, si_dist)
  r_corrected <- r_estimate * correction

  # Cap extremely high values
  r_corrected[r_corrected > 10] <- NA

  return(list(r = r_estimate, r_corrected = r_corrected))
}

#' Calculate Bootstrap Confidence Intervals
#'
#' Calculates confidence intervals using bootstrap resampling.
#'
#' @param incidence Vector of case counts
#' @param si_prob Matrix of serial interval probabilities
#' @param dates Vector of dates
#' @param si_mean Mean of the serial interval
#' @param si_sd Standard deviation of the serial interval
#' @param si_dist Distribution type
#' @param smoothing Smoothing window size
#' @param n_bootstrap Number of bootstrap samples
#' @param conf_level Confidence level
#' @return List with lower and upper bounds for R and R_corrected
#'
calculate_bootstrap_ci <- function(incidence, si_prob, dates, si_mean, si_sd, si_dist,
                                   smoothing, n_bootstrap, conf_level) {
  n <- length(incidence)

  # Initialize matrices to store bootstrap results
  bootstrap_r <- matrix(NA, nrow = n_bootstrap, ncol = n)
  bootstrap_r_corrected <- matrix(NA, nrow = n_bootstrap, ncol = n)

  # Generate bootstrap samples
  for (i in 1:n_bootstrap) {
    bootstrap_incidence <- generate_case_bootstrap(incidence)
    bootstrap_estimates <- calculate_r_estimates(
      bootstrap_incidence, si_prob, dates, si_mean, si_sd, si_dist, smoothing
    )
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

  return(list(
    r_lower = r_lower,
    r_upper = r_upper,
    r_corrected_lower = r_corrected_lower,
    r_corrected_upper = r_corrected_upper
  ))
}

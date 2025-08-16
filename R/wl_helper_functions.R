#' Calculate Serial Interval Probability Matrix
#'
#' Computes a matrix of transmission probabilities between all pairs of cases based
#' on their time differences and the specified serial interval distribution. Only
#' considers epidemiologically plausible transmission pairs (earlier to later cases).
#'
#' @param day_diffs numeric matrix; matrix of day differences between each pair of cases, where element \code{[i,j]} represents days between case i and case j
#' @param si_mean numeric; mean of the serial interval distribution in days
#' @param si_sd numeric; standard deviation of the serial interval distribution in days
#' @param si_dist character; distribution type, either "gamma" or "normal"
#' @return numeric matrix; matrix of transmission probabilities where element \code{[i,j]}
#'         represents the probability that case j infected case i based on their
#'         time difference and the serial interval distribution
#' @export
#' @examples
#' # Create sample day differences matrix
#' dates <- as.Date(c("2023-01-01", "2023-01-03", "2023-01-05"))
#' day_diffs <- create_day_diff_matrix(dates)
#'
#' # Calculate probability matrix
#' prob_matrix <- calculate_si_probability_matrix(day_diffs, si_mean = 7, si_sd = 3, si_dist = "gamma")
#'
calculate_si_probability_matrix <- function(
  day_diffs,
  si_mean,
  si_sd,
  si_dist
) {
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
#' Creates a symmetric matrix containing the time differences (in days) between
#' all pairs of cases based on their symptom onset dates.
#'
#' @param dates vector; dates of symptom onset for each case. Can be Date objects
#'              or any format coercible to dates
#' @return numeric matrix; symmetric matrix where element \code{[i,j]} represents the
#'         number of days between case i and case j (positive if i occurs after j)
#' @export
#' @examples
#' # Create day difference matrix from onset dates
#' onset_dates <- as.Date(c("2023-01-01", "2023-01-04", "2023-01-07", "2023-01-10"))
#' day_differences <- create_day_diff_matrix(onset_dates)
#' print(day_differences)
#'
create_day_diff_matrix <- function(dates) {
  n <- length(dates)
  day_diffs <- matrix(NA, nrow = n, ncol = n)

  for (i in 1:n) {
    for (j in 1:n) {
      day_diffs[i, j] <- as.numeric(difftime(
        dates[i],
        dates[j],
        units = "days"
      ))
    }
  }

  return(day_diffs)
}

#' Apply Moving Average Smoothing to R Estimates
#'
#' Applies temporal smoothing to reproduction number estimates using a centered
#' moving average window. Handles missing and infinite values appropriately.
#'
#' @param r_estimate numeric vector; reproduction number estimates to smooth.
#'                   Can contain NA or infinite values
#' @param window integer; size of the smoothing window in time units. Window is
#'               centered around each point
#' @return numeric vector; smoothed reproduction number estimates of the same
#'         length as input. Returns NA for points with insufficient valid
#'         neighboring values
#'
#' @examples
#' # This function is used internally for bootstrap confidence intervals
#' # For bootstrap reproduction number estimates, see ?wallinga_lipsitch
#' \dontrun{
#' # Smooth noisy R estimates
#' noisy_r <- c(1.2, 3.5, 1.8, 2.1, 1.6, 4.2, 1.9, 1.4)
#' smoothed_r <- smooth_estimates(noisy_r, window = 3)
#' }
#'
smooth_estimates <- function(r_estimate, window) {
  # Replace NA and infinite values with NA
  r_estimate[!is.finite(r_estimate)] <- NA

  n <- length(r_estimate)
  r_smooth <- rep(NA, n)
  half_window <- floor(window / 2)

  for (i in 1:n) {
    start_idx <- max(1, i - half_window)
    end_idx <- min(n, i + half_window)

    window_values <- r_estimate[start_idx:end_idx]
    if (sum(!is.na(window_values)) > 0) {
      r_smooth[i] <- mean(window_values, na.rm = TRUE)
    }
  }

  return(r_smooth)
}

#' Calculate Right-Truncation Correction Factors
#'
#' Computes correction factors to adjust reproduction number estimates for
#' right-truncation bias. This bias occurs because cases near the end of the
#' observation period may have generated secondary cases that are not yet observed.
#'
#' @param dates vector; dates corresponding to each case
#' @param si_mean numeric; mean of the serial interval distribution in days
#' @param si_sd numeric; standard deviation of the serial interval distribution in days
#' @param si_dist character; distribution type, either "gamma" or "normal"
#' @return numeric vector; correction factors for each case. Values > 1 indicate
#'         upward adjustment needed. Returns NA when correction would be unreliable
#'         (probability of observation <= 0.5)
#' @export
#' @examples
#' # Calculate truncation correction for recent cases
#' case_dates <- seq(as.Date("2023-01-01"), as.Date("2023-01-20"), by = "day")
#' corrections <- calculate_truncation_correction(
#'   case_dates, si_mean = 7, si_sd = 3, si_dist = "gamma"
#'   )
#'
#' # Show how correction increases for more recent cases
#' tail(corrections, 5)
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

#' Generate Bootstrap Sample of Case Incidence
#'
#' Creates a bootstrap sample by resampling individual cases with replacement,
#' then reconstructing daily incidence counts. This maintains the temporal
#' distribution while introducing sampling variation for uncertainty estimation.
#'
#' @param incidence numeric vector; daily case counts (non-negative integers)
#' @return numeric vector; bootstrapped daily incidence of the same length as input.
#'         Total number of cases remains the same but their temporal distribution varies
#'
#' @examples
#' # This function is used internally for bootstrap confidence intervals
#' # For bootstrap reproduction number estimates, see ?wallinga_lipsitch
#'
#' \dontrun{
#' # Example usage through main interface:
#' dates <- seq(as.Date("2023-01-01"), by = "day", length.out = 10)
#' incidence <- c(1, 3, 5, 8, 12, 15, 10, 6, 3, 1)
#'
#' # Get bootstrap confidence intervals
#' results <- wallinga_lipsitch(
#'   incidence = incidence,
#'   dates = dates,
#'   si_mean = 7,
#'   si_sd = 3,
#'   si_dist = "gamma",
#'   bootstrap = TRUE,
#'   n_bootstrap = 100
#' )
#'
#' # The bootstrap sampling happens automatically inside wallinga_lipsitch()
#' head(results)
#' }
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

#' Calculate Reproduction Number Estimates
#'
#' Implements the Wallinga-Lipsitch algorithm to estimate case reproduction
#' numbers from incidence data and serial interval probabilities. This function
#' performs the likelihood calculations for retrospective reproduction
#' number estimation.
#'
#' The algorithm calculates the probability that each earlier case infected each
#' later case based on their time difference and the serial interval distribution.
#' These probabilities are then aggregated to estimate the expected number of
#' secondary cases generated by cases on each day.
#'
#' @param incidence numeric vector; daily case counts. Must be non-negative integers.
#'                   Days with zero cases will have R estimates of NA
#' @param si_prob numeric matrix; serial interval probability matrix from
#'                \code{\link{calculate_si_probability_matrix}}. Element \code{[i,j]}
#'                represents the probability that case j infected case i
#' @param dates vector; dates corresponding to incidence data. Used for
#'              right-truncation correction calculations
#' @param si_mean numeric; mean of the serial interval distribution in days
#' @param si_sd numeric; standard deviation of the serial interval distribution in days
#' @param si_dist character; distribution type for serial interval, either "gamma"
#'                or "normal"
#' @param smoothing integer; window size for temporal smoothing (0 = no smoothing).
#'                  When > 1, applies centered moving average to reduce noise
#'
#' @return named list with two numeric vectors of the same length as \code{incidence}:
#' \itemize{
#'   \item \code{r}: Raw case reproduction number estimates. Returns NA for days
#'         with zero cases or single-case epidemics
#'   \item \code{r_corrected}: Estimates with right-truncation correction applied.
#'         Values > 10 are capped at NA to avoid unrealistic estimates
#' }
#'
#' @details
#' The Wallinga-Lipsitch method works by:
#' \enumerate{
#'   \item Computing transmission likelihoods from earlier to later cases
#'   \item Normalizing these likelihoods to create proper probabilities
#'   \item Aggregating probabilities to estimate expected secondary cases per primary case
#'   \item Applying right-truncation correction for cases near the observation end
#' }
#'
#' The right-truncation correction accounts for the fact that cases near the end
#' of the observation period may have generated secondary cases that occur after
#' data collection ended.
#'
#' @seealso \code{\link{wallinga_lipsitch}} for the main user interface,
#'          \code{\link{calculate_si_probability_matrix}} for probability matrix creation,
#'          \code{\link{calculate_truncation_correction}} for correction details
#'
#' @examples
#' # This function is used internally by wallinga_lipsitch()
#' # For complete reproduction number estimation, see ?wallinga_lipsitch
#'
#' \dontrun{
#' # Example usage through main interface:
#' dates <- seq(as.Date("2023-01-01"), by = "day", length.out = 10)
#' incidence <- c(1, 2, 4, 6, 8, 6, 4, 2, 1, 0)
#'
#' results <- wallinga_lipsitch(
#'   incidence = incidence,
#'   dates = dates,
#'   si_mean = 7,
#'   si_sd = 3,
#'   si_dist = "gamma"
#' )
#'
#' # Access the reproduction number estimates:
#' results$R_corrected
#' }
#'
calculate_r_estimates <- function(
  incidence,
  si_prob,
  dates,
  si_mean,
  si_sd,
  si_dist,
  smoothing
) {
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
      if (sum(incidence) == incidence[j]) {
        r_estimate[j] <- NA
      } else {
        r_estimate[j] <- sum(incidence * rel_likelihood[, j]) / incidence[j]
      }
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

  return(list(
    r = r_estimate,
    r_corrected = r_corrected
  ))
}

#' Calculate Bootstrap Confidence Intervals for R Estimates
#'
#' Generates bootstrap confidence intervals for reproduction number estimates by
#' resampling the incidence data multiple times and calculating quantiles of the
#' resulting R distributions.
#'
#' @param incidence numeric vector; daily case counts
#' @param si_prob numeric matrix; serial interval probability matrix
#' @param dates vector; dates corresponding to incidence data
#' @param si_mean numeric; mean of the serial interval distribution
#' @param si_sd numeric; standard deviation of the serial interval distribution
#' @param si_dist character; distribution type, either "gamma" or "normal"
#' @param smoothing integer; window size for temporal smoothing
#' @param n_bootstrap integer; number of bootstrap samples to generate
#' @param conf_level numeric; confidence level (between 0 and 1)
#' @return named list with confidence interval bounds:
#' \itemize{
#'   \item \code{r_lower, r_upper}: Confidence intervals for raw R estimates
#'   \item \code{r_corrected_lower, r_corrected_upper}: Confidence intervals for corrected R estimates
#' }
#'
#' @examples
#' # This function is used internally by wallinga_lipsitch() when bootstrap=TRUE
#' # For complete examples with confidence intervals, see ?wallinga_lipsitch
#'
#' \dontrun{
#' # Example usage through main interface:
#' dates <- seq(as.Date("2023-01-01"), by = "day", length.out = 10)
#' incidence <- c(1, 2, 4, 6, 8, 6, 4, 2, 1, 0)
#'
#' results <- wallinga_lipsitch(
#'   incidence = incidence,
#'   dates = dates,
#'   si_mean = 7,
#'   si_sd = 3,
#'   si_dist = "gamma",
#'   bootstrap = TRUE,
#'   n_bootstrap = 100
#' )
#'
#' # Access confidence intervals:
#' results$R_corrected_lower
#' results$R_corrected_upper
#' }
#'
calculate_bootstrap_ci <- function(
  incidence,
  si_prob,
  dates,
  si_mean,
  si_sd,
  si_dist,
  smoothing,
  n_bootstrap,
  conf_level
) {
  n <- length(incidence)

  # Initialize matrices to store bootstrap results
  bootstrap_r <- matrix(NA, nrow = n_bootstrap, ncol = n)
  bootstrap_r_corrected <- matrix(NA, nrow = n_bootstrap, ncol = n)

  # Generate bootstrap samples
  for (i in 1:n_bootstrap) {
    bootstrap_incidence <- generate_case_bootstrap(incidence)
    bootstrap_estimates <- calculate_r_estimates(
      bootstrap_incidence,
      si_prob,
      dates,
      si_mean,
      si_sd,
      si_dist,
      smoothing
    )
    bootstrap_r[i, ] <- bootstrap_estimates$r
    bootstrap_r_corrected[i, ] <- bootstrap_estimates$r_corrected
  }

  # Calculate confidence intervals
  alpha <- (1 - conf_level) / 2
  quantiles <- c(alpha, 1 - alpha)

  # For each day, calculate quantiles of bootstrap distributions
  r_lower <- apply(bootstrap_r, 2, function(x) {
    quantile(x, probs = quantiles[1], na.rm = TRUE)
  })
  r_upper <- apply(bootstrap_r, 2, function(x) {
    quantile(x, probs = quantiles[2], na.rm = TRUE)
  })

  r_corrected_lower <- apply(bootstrap_r_corrected, 2, function(x) {
    quantile(x, probs = quantiles[1], na.rm = TRUE)
  })
  r_corrected_upper <- apply(bootstrap_r_corrected, 2, function(x) {
    quantile(x, probs = quantiles[2], na.rm = TRUE)
  })

  return(list(
    r_lower = r_lower,
    r_upper = r_upper,
    r_corrected_lower = r_corrected_lower,
    r_corrected_upper = r_corrected_upper
  ))
}

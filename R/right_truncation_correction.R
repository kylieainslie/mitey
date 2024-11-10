#' Correct for right truncation in estimates of time-varying reproduction number
#'
#' This function calculates the correction factor for each observed time of symptom onset, t. It uses the cumulative distribution function (`pnorm`) of the normal distribution to calculate the probability that the serial interval is less than the time lag T_now - t, where T_now is the current time.
#'
#' @param t integer; time point of infection or symptom onset
#' @param T_now interger; current time point
#' @param mean_serial_interval numeric; mean of serial interval distribution
#' @param sd_serial_interval numeric; standard deviation of serial interval distribution
#' @param distribution character string; assumed serial interval distribution, takes arguments "normal" or "gamma".
#' @return vector of correction weights
#' @export
#' @importFrom stats pnorm
#' @importFrom stats pgamma
right_truncation_correction <- function(t, T_now, mean_serial_interval,
                                        sd_serial_interval, distribution = "normal") {
  if (distribution == "normal") {
    # Probability that the serial interval is less than the time lag (T - t) using normal distribution
    prob <- pnorm((T_now - t - mean_serial_interval) / sd_serial_interval)
  } else if (distribution == "gamma") {
    # Assuming mean and sd to shape and scale parameters for gamma distribution
    shape <- (mean_serial_interval / sd_serial_interval)^2
    scale <- (sd_serial_interval^2) / mean_serial_interval
    # Probability using gamma distribution
    prob <- pgamma(T_now - t, shape = shape, scale = scale)
  } else {
    stop("Unsupported distribution. Use 'normal' or 'gamma'.")
  }

  w <- 1 / prob
  return(w)  # Correction factor
}


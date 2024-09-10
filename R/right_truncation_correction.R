#' Correct for right truncation in estimates of time-varying reproduction number
#'
#' This function calculates the correction factor for each observed time of symptom onset, t. It uses the cumulative distribution function (`pnorm`) of the normal distribution to calculate the probability that the serial interval is less than the time lag T_now - t, where T_now is the current time.
#'
#' @param t integer; time point of infection or symptom onset
#' @param T_now interger; current time point
#' @param mean_serial_interval numeric; mean of serial interval distribution
#' @param sd_serial_interval numeric; standard deviation of serial interval distribution
#' @return vector of correction weights
#' @export
#' @importFrom stats pnorm
right_truncation_correction <- function(t, T_now, mean_serial_interval, sd_serial_interval) {
  # Probability that the serial interval is less than the time lag (T - t)
  prob <- pnorm((T_now - t - mean_serial_interval) / sd_serial_interval)
  w <- 1 / prob
  # Return the correction factor
  return(w)  # Correction factor as per Cauchemez et al.
}

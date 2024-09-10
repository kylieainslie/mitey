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
#' @param n_bootstrap integer; number of bootstrap samples of the serial interval distribution
#' @return A numeric vector with the expected reproduction number for each day.
#' @export
#' @importFrom purrr map
#' @importFrom stats quantile
#' @importFrom stats rnorm
#' @importFrom stats rgamma
rt_estim_w_boot <- function(inc_dat, mean_si, sd_si, dist_si = "normal",
                            cut_tail = NULL, pos_only = TRUE, n_bootstrap = 100){

  # get bootstrap samples
  boot_samples <- map(1:n_bootstrap, ~inc_dat[sample(nrow(inc_dat), replace = TRUE), ])

  # Apply rt_estim to each bootstrap sample
  boot_rt <- map(boot_samples, ~rt_estim(inc_dat = .x, mean_si = mean_si, sd_si = sd_si, dist_si = dist_si,
                                         cut_tail = cut_tail, pos_only = pos_only))

  # Calculate the mean R_t and confidence intervals across bootstrap samples
  rt_mean <- lapply(boot_rt, 2, mean)
  rt_lower <- lapply(boot_rt, 2, quantile, probs = 0.025)
  rt_upper <- lapply(bootstrap_rt, 2, quantile, probs = 0.975)



}

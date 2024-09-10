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
#' @import dplyr
#' @import tidyr
#' @importFrom purrr map
#' @importFrom stats quantile
#' @importFrom stats rnorm
#' @importFrom stats rgamma
rt_estim_w_boot <- function(inc_dat, mean_si, sd_si, dist_si = "normal",
                            cut_tail = NULL, pos_only = TRUE, n_bootstrap = 100){

  # uncount inc_dat, so that individual cases can be samples from
  inc_dat_uncount <- inc_dat %>%
    uncount(inc) %>%
    mutate(inc = 1)

  # get bootstrap samples
  boot_samples <- map(1:n_bootstrap, ~inc_dat_uncount[sample(nrow(inc_dat_uncount), replace = TRUE), ])

  # a little wrangling of each bootstrapped data frame:
  # 1. order by onset_date
  # 2. collapse rows with the same onset_date and increase the count of inc
  # 3. complete the time series by filling in missing dates
  # get all dates
  all_dates <- seq(min(inc_dat$onset_date), max(inc_dat$onset_date), by = "day")

  boot_samples_wrangled <- map(boot_samples, ~.x %>%
                                arrange(onset_date) %>%
                                group_by(onset_date) %>%
                                summarise(inc = sum(inc), .groups = "drop") %>%
                                complete(onset_date = all_dates,
                                         fill = list(inc = 0)
                                )
                             )

  # Apply rt_estim to each bootstrap sample
  boot_rt <- map(boot_samples_wrangled, ~rt_estim(inc_dat = .x, mean_si = mean_si, sd_si = sd_si, dist_si = dist_si,
                                         cut_tail = cut_tail, pos_only = pos_only))

  # Calculate the mean R_t and confidence intervals across bootstrap samples
  rt_mean <- lapply(boot_rt, 2, mean)
  rt_lower <- lapply(boot_rt, 2, quantile, probs = 0.025)
  rt_upper <- lapply(bootstrap_rt, 2, quantile, probs = 0.975)

  return(lis)

}

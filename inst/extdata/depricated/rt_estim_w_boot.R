#' Estimate Time-varying Reproduction Number
#'
#' This function estimates the time-varying reproduction number (R_t) using the method proposed by Wallinga and Teunis (AJE 2004).
#' The function takes as input a data frame with incidence data and requires the specification of the serial interval distribution.
#' The user must input the mean and standard deviation of the serial interval distribution and select the functional form: Normal or Gamma.
#'
#' @param inc_dat data frame; data frame with incidence data. The data frame should have two columns: inc (daily incidence) and onset_date. Onset_date does not have to be in date format, it can be a column of days, such as c(1,2,3,4, ... ),
#' @param mean_si numeric; mean of serial interval distribution
#' @param sd_si numeric; standard deviation of serial interval distribution
#' @param dist_si string; distribution to be assumed for serial interval. Accepts "normal" or "gamma".
#' @param n_bootstrap integer; number of bootstrap samples of the serial interval distribution
#' @return A named list of two data frames. The first data frame (`results`) contains the mean, mediean, 2.5th percentile, and 97.5th percentile from all of the boot strapped samples. The second data frame (`boot_samples`) contains all boot strap samples for each time point.
#' @export
#' @import dplyr
#' @import tidyr
#' @importFrom purrr map
#' @importFrom stats quantile
#' @importFrom stats rnorm
#' @importFrom stats rgamma
#' @importFrom stats median
rt_estim_w_boot <- function(inc_dat, mean_si, sd_si, dist_si = "normal",
                            n_bootstrap = 100){

  # uncount inc_dat, so that individual cases can be samples from
  if(!is.integer(inc_dat$inc)){
    warning("The inputted incidence data is not an integer. Incidence will be rounded to the nearest whole number. This may affect results.")
  }

  inc_dat_uncount <- inc_dat %>%
    mutate(inc = as.integer(round(.data$inc))) %>%
    uncount(.data$inc) %>%
    mutate(inc = 1)

  # get bootstrap samples
  boot_samples <- map(1:n_bootstrap,
                      ~inc_dat_uncount[sample(nrow(inc_dat_uncount), replace = TRUE), ])

  # a little wrangling of each bootstrapped data frame:
  # 1. order by onset_date
  # 2. collapse rows with the same onset_date and increase the count of inc
  # 3. complete the time series by filling in missing dates
  # get all dates
  if(is.Date(inc_dat$onset_date)){
    all_dates <- seq(min(inc_dat$onset_date), max(inc_dat$onset_date), by = "day")
  } else {
    all_dates <- seq(min(inc_dat$onset_date), max(inc_dat$onset_date), by = 1)
  }

  boot_samples_wrangled <- map(boot_samples, ~.x %>%
                                arrange(.data$onset_date) %>%
                                group_by(.data$onset_date) %>%
                                summarise(inc = sum(inc), .groups = "drop") %>%
                                complete(onset_date = all_dates,
                                         fill = list(inc = 0)
                                )
                             )

  # Apply rt_estim to each bootstrap sample
  boot_rt <- map(boot_samples_wrangled,
                 ~rt_estim(inc_dat = .x, mean_si = mean_si, sd_si = sd_si,
                           dist_si = dist_si))

  df_boot_rt <- bind_rows(boot_rt, .id = "sample")

  df_boot_results <- df_boot_rt %>%
    group_by(.data$onset_date) %>%
    summarise(
      mean_rt = mean(.data$rt, na.rm = TRUE),
      median_rt = median(.data$rt, na.rm = TRUE),
      lower_rt = quantile(.data$rt, 0.025, na.rm = TRUE),
      upper_rt = quantile(.data$rt, 0.975, na.rm = TRUE),
      mean_rt_adjusted = mean(.data$rt_adjusted, na.rm = TRUE),
      median_rt_adjusted = median(.data$rt_adjusted, na.rm = TRUE),
      lower_rt_adjusted = quantile(.data$rt_adjusted, 0.025, na.rm = TRUE),
      upper_rt_adjusted = quantile(.data$rt_adjusted, 0.975, na.rm = TRUE)
    )

  return(list(results = df_boot_results, boot_samples = df_boot_rt))
}

# Calculate reproduction number from scabies time series

# 1. We have reported cases of scabies reported weekly. We will first randomly
#    assign each date of symptom onset to a day of the reported week. This will
#    give us a daily time series instead of weekly.
#
# 2. We will use the method proposed by Wallinga and Teunis 2004 to estimate the
#    time-varying reproduction number
#    a. we will first assume a Normal serial interval distribution
#    b. next, we will assume a Gamma serial interval distribution

# load required packages -------------------------------------------------------
library(dplyr)
library(tidyr)
library(lubridate)
library(ISOweek)
library(readxl)
library(ggplot2)
library(zoo)
library(EpiEstim)
library(devtools)
load_all()

# read in data -----------------------------------------------------------------
nivel_wkly_data <-
  read_xlsx("~/Dropbox/Kylie/Projects/RIVM/Projects/scabies/data/scabies_data_weekly.xlsx") %>%
  # fix/translate variable names
  rename(diagnosis_code = `Diagnose (ICPC)`,
         year = `ISO-jaar`,
         week_num = `ISO-weeknr (ma-zo)`,
         pop_size = `Aantal populatie`,
         cases = `Aantal prevalente cases`,
         prev_per_100000 = `Prevalentie per 100.000`) %>%
  # drop diagnosis var
  select(-diagnosis_code) %>%
  # create new var that combines year and week
  mutate(yr_wk = paste(year, week_num, sep = "_"),
         year = as.factor(year))

# 1. randomly assign a date within reporting date for symptom onset date -------
nivel_daily_data <- nivel_wkly_data %>%
  uncount(cases) %>% # Repeat rows based on the number of cases
  mutate(
    iso_week = paste0(year, "-W", sprintf("%02d", as.numeric(week_num))),
    first_day = ISOweek2date(paste0(iso_week, "-1")),
    random_day = sample(0:6, n(), replace = TRUE),
    onset_date = first_day + days(random_day)
  ) %>%
  select(-iso_week, -first_day, -random_day)

# 2. We will implement the method from Wallinga and Teunis 2004 ----------------

# Determine likelihood of transmission pair for every possible pair
# we assume the serial interval distribution is N(95.57, (15.17)^2)
# we will loop over each observation and calculate the likelihood of the
# transmission pair.
n <- dim(nivel_daily_data)[1]

# instead of calculating the probability of every pair of persons, we will calculate
# the likelihood of an event occurring for every pair of time points. This will
# decrease the dimensions of our likelihood matrix.

# Group data by onset_date and calculate the daily incidence
df <- nivel_daily_data %>%
  group_by(onset_date) %>%
  mutate(count = n()) %>%
  distinct(onset_date, count) %>%
  arrange(onset_date) %>%
  mutate(num_date = as.numeric(onset_date)) %>%
  ungroup() %>%
  rename(inc = count) %>%
  select(onset_date, inc)

# si dist parameters
si_mean <- 123.24
si_sd <- 31.55

# Calculate the initial R_t
rt_estimates <- rt_estim(df, mean_si = si_mean, sd_si = si_sd, dist_si = "normal")

# quick visual check!
rt_estimates1 <- rt_estimates %>%
  filter(!is.nan(rt),
         onset_date > min(onset_date) + 123, # trim first serial interval
         onset_date < max(onset_date) - 123  # trim last serial interval
         ) %>%
  mutate(
    rt_rollmean = rollmean(rt, 35, fill = NA),
    rt_adj_rollmean = rollmean(rt_adjusted, 35, fill = NA)
    )

ggplot(data = rt_estimates1, aes(x = onset_date, y = rt_rollmean)) +
  geom_line() +
  geom_hline(yintercept = 1.1, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = 1, linetype = "solid", color = "black") +
  theme_minimal()

# bootstrap to get 95% CIs
rt_bootstrap <- rt_estim_w_boot(df, mean_si = si_mean, sd_si = si_sd, dist_si = "normal",
                               n_bootstrap = 100)

# plot
rt_df_for_plot <- left_join(rt_estimates, rt_bootstrap$results, by = "onset_date") %>%
  # filter out first serial interval
  filter(onset_date > min(onset_date) + 123,
         onset_date < max(onset_date) - 21)
saveRDS(rt_df_for_plot, "vignettes/data/rt_df_for_plot.rds")

# compare rt unadjusted to adjusted
rt_plot_comparison <- ggplot(rt_df_for_plot, aes(x = onset_date)) +
  geom_line(aes(y = rollmean(rt, 35, fill = NA), color = "Rt unadjusted")) +
  #geom_line(aes(y = rollmean(rt_adjusted, 35, fill = NA), color = "Rt adjusted")) +
  #geom_ribbon(aes(ymin = rollmean(lower_rt, 35, fill = NA), ymax = rollmean(upper_rt, 35, fill = NA)), alpha = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  # Customizing the plot
  labs(
    x = "Date of symptom onset",
    y = "Reproduction number (Rt)",
    color = "Legend"
  ) +
  theme_minimal()

# plot rt_adjusted with confidence bounds (and season bars)
# Breaks for background rectangles
start_breaks <- c(as.Date("2011-01-01"),
                  seq(as.Date("2011-03-01"), tail(rt_df_for_plot$onset_date, 1), by = "quarter"))
end_breaks <- c(start_breaks[-1] - 1, tail(rt_df_for_plot$onset_date,1))
rects <- data.frame(xstart = start_breaks,
                    xend = end_breaks,
                    season = factor(c(rep(c("winter", "spring", "summer", "autumn"),
                                          length.out = length(start_breaks))),
                                    levels = c("winter", "spring", "summer", "autumn")))

p <- ggplot() +
  geom_rect(data = rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf,
                              fill = season), alpha = 0.4) +
  scale_fill_viridis_d() +
  geom_line(data = rt_df_for_plot, aes(x = onset_date,
                                       y = rollmean(rt_adjusted, 35, fill = NA))) +
  geom_ribbon(data = rt_df_for_plot, aes(x = onset_date,
                                         ymin = rollmean(lower_rt_adjusted, 35, fill = NA),
                                         ymax = rollmean(upper_rt_adjusted, 35, fill = NA)),
              alpha = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = "Date of Symptom Onset", y = "Reproduction Number") +
  theme(
    panel.background = element_blank()
  )

p

# plot epidemic curve
epidemic_curve <- ggplot(df, aes(x = onset_date, y = inc)) +
  geom_col(fill = "steelblue") +
  #geom_smooth(method = "loess", color = "red", se = FALSE) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    x = "Date of Diagnosis",
    y = "Incidence of scabies infections per 1000 people"
  ) +
  theme_minimal()
epidemic_curve
ggsave("vignettes/figures/nivel_epidemic_curve.png", plot = epidemic_curve, width = 7, height = 5, dpi = 300)

# Check rt estimates -----------------------------------------------------------
# 1. Calculate b(t+u)/b(t) where b() is the incidence, t is the current time
#    point, and u is the mean SI
seq1 <- 1:(nrow(rt_estimates1))
seq2 <- seq1 + 123
rt_check1 <- data.frame(
  onset_date = (seq2 - 1) + min(rt_estimates1$onset_date),
  rt_check1 = rt_estimates1$rt[seq2]/rt_estimates1$rt[seq1]
)

# 2. Calculate geometric mean
# This should be approx R0 = 1.1
# First remove last serial interval
rt_check2 <- exp(mean(log(rt_estimates2$rt), na.rm = TRUE))

# 3. Calculate rt using EpiEstim
# use discritized normal distribution
discretize_normal <- function(mean, sd, lower, upper) {
  bins <- lower:upper
  probs <- diff(pnorm(c(bins, Inf), mean, sd))  # Probabilities for each bin
  return(probs)
}

# Example usage
distr_si <- discretize_normal(mean = si_mean, sd = si_sd, lower = 0, upper = 250)

window <- 1 # time window (in days) over which to calculate Rt
ts <- 1:nrow(df)
ts <- ts[ts > 1 & ts <= (max(ts)-window+1)]
te <- ts+(window-1)

config <- make_config(incid = df$inc,
                      si_distr = c(0,distr_si),
                      t_start = ts,
                      t_end = te)
config$n_sim <- 100

rt_epiestim <- wallinga_teunis(
  incid = df$inc,
  method = "non_parametric_si",
  config = config)

rt_check3 <- rt_epiestim$R %>%
  mutate(onset_date = t_start + min(df$onset_date)) %>%
  select(onset_date, `Mean(R)`)

# plot checks
df_check <- rt_estimates1 %>%
  select(onset_date, rt) %>%
  left_join(., rt_check1, by = "onset_date") %>%
  left_join(., rt_check3, by = "onset_date") %>%
  pivot_longer(-onset_date) %>%
  mutate(Method = factor(
    case_when(
      name == "rt" ~ "WT (mitey)",
      name == "Mean(R)" ~ "WT (EpiEstim)",
      name == "rt_check1" ~ "Growth Rate"
    ), levels = c("Growth Rate", "WT (mitey)", "WT (EpiEstim)")
  )
  )

p_check <- ggplot(data = df_check,
                  aes(x = onset_date, y = value, color = Method)) +
  geom_line() +
  geom_hline(yintercept = rt_check2, linetype = "dashed", col = "grey") +
  scale_color_viridis_d(option = "D") +
  labs(title = "Compare Methods To Estimate Rt",
       x = "Onset Date",
       y = "Daily Reproduction Number",
       color = "Method") +
  theme_minimal()
p_check

# write results to file
#saveRDS(rt_epiestim$R, "vignettes/data/rt_epiestim_gamma_si.rds")
saveRDS(rt_epiestim$R, "vignettes/data/rt_epiestim_normal_si.rds")

### Sensitivity analysis -------------------------------------------------------
# assume Gamma distributed SI distribution
# Calculate the initial R_t
rt_estimates_sa <- rt_estim(df, mean_si = si_mean, sd_si = si_sd, dist_si = "gamma")

# bootstrap to get 95% CIs
rt_bootstrap_sa <- rt_estim_w_boot(df, mean_si = si_mean, sd_si = si_sd,
                               dist_si = "gamma", n_bootstrap = 100)

rt_df_for_plot_sa <- left_join(rt_estimates_sa, rt_bootstrap_sa$results, by = "onset_date") %>%
  # filter out first serial interval
  filter(onset_date > min(onset_date) + 123,
         onset_date < as.Date("2024-01-01", format = "%Y-%m-%d"))

# calculate geometric mean
geom_mean_rt_sa <- exp(mean(log(rt_df_for_plot_sa$rt_adjusted[rt_df_for_plot_sa$rt_adjusted != 0])))

# SA plot -------
start_breaks <- c(as.Date("2011-01-01"),
                  seq(as.Date("2011-03-01"), tail(rt_df_for_plot_sa$onset_date, 1), by = "quarter"))
end_breaks <- c(start_breaks[-1] - 1, tail(rt_df_for_plot_sa$onset_date,1))
rects <- data.frame(xstart = start_breaks,
                    xend = end_breaks,
                    season = factor(c(rep(c("winter", "spring", "summer", "autumn"),
                                          length.out = length(start_breaks))),
                                    levels = c("winter", "spring", "summer", "autumn")))

p_sa <- ggplot() +
  geom_rect(data = rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf,
                              fill = season), alpha = 0.4) +
  scale_fill_viridis_d() +
  geom_line(data = rt_df_for_plot_sa, aes(x = onset_date,
                                       y = rollmean(rt_adjusted, 35, fill = NA))) +
  geom_ribbon(data = rt_df_for_plot_sa, aes(x = onset_date,
                                         ymin = rollmean(lower_rt_adjusted, 35, fill = NA),
                                         ymax = rollmean(upper_rt_adjusted, 35, fill = NA)),
              alpha = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = "Date of Symptom Onset", y = "Reproduction Number") +
  theme(
    panel.background = element_blank()
  )

p_sa
### Checks ---------------------------------------------------------------------

# estimate using EpiEstim (specifying the Wallinga and Teunis method) as a check

# EpiEstim SI distribution discretization:
# assumes shifted gamma
si_mean <- 123.24
si_sd <- 31.55
si_max <- 250
distr_si <- discr_si(seq(0, si_max), mu = si_mean, sigma = si_sd)
sum(distr_si*(seq(0, si_max)))


# ------------------------
# load R estimates
rt_epiestim_gamma <- readRDS("vignettes/data/rt_epiestim_gamma_si.rds")
rt_epiestim_normal <- readRDS("vignettes/data/rt_epiestim_normal_si.rds")

# do some data wrangling
# df_rt <- rt_epiestim_gamma %>%
#   mutate(onset_date = )

# adjust
onset_times <- 1:nrow(df_rt_check)
correction_factors <- sapply(onset_times, function(t) {
  right_truncation_correction(t, nt, si_mean, si_sd, distribution = "gamma")
})

# Adjusted Rt values for right truncation
df_rt_check$adjusted_rt <- df_rt_check$`Mean(R)` * correction_factors
df_rt_check$adjusted_rt_lower <- df_rt_check$`Quantile.0.025(R)` * correction_factors
df_rt_check$adjusted_rt_upper <- df_rt_check$`Quantile.0.975(R)` * correction_factors

# calculate mean
mean_rt <- mean(df_rt_check$adjusted_rt, na.rm = TRUE)
mean_rt

# plot
# trim data set
df_plot <- df_rt_check %>%
  #filter()
  filter(#onset_date < max(onset_date) - 123,
         onset_date < as.Date("2024-01-01", format = "%Y-%m-%d"),
         onset_date > min(onset_date) + 123)

start_breaks <- c(as.Date("2011-01-01"),
                  seq(as.Date("2011-03-01"), tail(df_plot$onset_date, 1), by = "quarter"))
end_breaks <- c(start_breaks[-1] - 1, tail(df_plot$onset_date,1))
rects <- data.frame(xstart = start_breaks,
                    xend = end_breaks,
                    season = factor(c(rep(c("winter", "spring", "summer", "autumn"),
                                          length.out = length(start_breaks))),
                                    levels = c("winter", "spring", "summer", "autumn")))

p_epiestim <- ggplot() +
  geom_rect(data = rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf,
                              fill = season), alpha = 0.4) +
  scale_fill_viridis_d() +
  geom_line(data = df_plot,
            aes(x = onset_date,y = adjusted_rt)) +
  geom_ribbon(data = df_plot,
              aes(x = onset_date, ymin = adjusted_rt_lower, ymax = adjusted_rt_upper),
              alpha = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = "Date of Symptom Onset", y = "Reproduction Number") +
  theme(
    panel.background = element_blank()
  )

p_epiestim

# Compare estimates graphically
# ggplot(df %>%
#          filter(onset_date > min(onset_date) + 123),
#        aes(x = onset_date)) +
#   # R_t with confidence intervals
#   #geom_ribbon(aes(ymin = rt_lower, ymax = rt_upper), fill = "blue", alpha = 0.2) +
#   geom_line(aes(y = rt_mean, color = "R_t"), size = 1) +
#   # Mean(R) with confidence intervals
#   # geom_ribbon(aes(ymin = `Quantile.0.025(R)`, ymax = `Quantile.0.975(R)`), fill = "red", alpha = 0.2) +
#   #geom_line(aes(y = `Mean(R)`, color = "Mean(R)`"), size = 1) +
#   # Horizontal line at Rt = 1
#   geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 0.8) +
#   # Customizing the plot
#   labs(
#     title = "R_t and Mean(R) with Confidence Intervals",
#     x = "Onset Date",
#     y = "Rt",
#     color = "Legend"
#   ) +
#   # Setting colors
#   #scale_color_manual(values = c("R_t" = "blue", "Mean(R)" = "red")) +
#   # Theme adjustments
#   theme_minimal()

# save r_t data set, so that we don't have to re-run bootstrapping
#saveRDS(df, "inst/extdata/data/rt_df.rds")


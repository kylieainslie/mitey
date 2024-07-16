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
# To begin, we will try a smaller example and only use data from the first week.
week1 <- nivel_daily_data %>%
  filter(year == 2011, week_num == 1) %>%
  # order by onset date
  arrange(onset_date)

# Determine likelihood of transmission pair for every possible pair
# we assume the serial interval distribution is N(95.57, (15.17)^2)
# we will loop over each observation and calculate the likelihood of the
# transmission pair.
n <- dim(nivel_daily_data)[1]

# instead of calculating the probability of every pair of persons, we will calculate
# the likelihood of an event occurring for every pair of time points. This will
# decrease the dimensions of our likelihood matrix.

# Step 1: Group data by onset_date and calculate the daily incidence
df <- nivel_daily_data %>%
  group_by(onset_date) %>%
  mutate(count = n()) %>%
  distinct(onset_date, count) %>%
  arrange(onset_date) %>%
  mutate(num_date = as.numeric(onset_date)) %>%
  ungroup()


# Step 2: Compute the likelihood matrix
nt <- nrow(df)

# Initialize the likelihood matrix
likelihood_mat <- matrix(0, nrow = nt, ncol = nt)

# Pre-compute the likelihood values for all possible serial intervals
serial_intervals <- outer(df$num_date, df$num_date, "-")

# Apply get_likelihood only to positive serial_intervals
likelihood_values <- apply(serial_intervals, 1:2, get_likelihood,
      mu = 95.57,
      sigma = 15.17,
      distn = "normal",
      tail_cut = 180,
      positive_only = TRUE
    )

# Step 2 & 3: Combine the calculation of the likelihood matrix and incorporate incidence data
for (i in 1:nt) {
  for (j in 1:nt) {
    if (i > j) {  # Only consider pairs where i > j
      likelihood_mat[i, j] <- likelihood_values[i, j] * df$count[j]
    }
  }
}


# Step 4: Get the marginal likelihood by summing each row
marginal_likelihood <- rowSums(likelihood_mat)

# Step 5: Calculate the probability matrix
prob_mat <- sweep(likelihood_mat, 1, marginal_likelihood, FUN = "/")
prob_mat[is.nan(prob_mat)] <- NA

# Step 6: Calculate the expected reproduction number per day (R_t)
expected_rt <- colSums(prob_mat, na.rm = TRUE)

# Add the R_t values to the data frame
df$R_t <- expected_rt

# plot -------------------------------------------------------------------------
df_plot <- df %>%
  # create new variable for plotting Rt where the first 95 days (equivalent to
  # 1 mean serial interval) are NA
  mutate(rt_plot = if_else(onset_date < min(onset_date) + 95, NA, R_t),
         rt_plot_rolling_avg = rollmean(rt_plot, k = 14, fill = NA, align = "right")
         ) %>%
  filter(onset_date < as.Date("2024-01-01"))

# Breaks for background rectangles
start_breaks <- c(as.Date("2011-01-01"),
                  seq(as.Date("2011-03-01"), tail(df_plot$onset_date, 1), by = "quarter"))
end_breaks <- c(start_breaks[-1] - 1, tail(df_plot$onset_date,1))
rects <- data.frame(xstart = start_breaks,
                    xend = end_breaks,
                    season = factor(c(rep(c("winter", "spring", "summer", "autumn"),
                                          length.out = length(start_breaks))),
                                    levels = c("winter", "spring", "summer", "autumn")))

p <- ggplot() +
  geom_rect(data = rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf,
                              fill = season), alpha = 0.4) +
  scale_fill_viridis_d() +
  geom_line(data = df_plot, aes(x = onset_date, y = rt_plot_rolling_avg)) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = "Date of Symptom Onset", y = "Reproduction Number") +
  theme(
    panel.background = element_blank()
  )

p

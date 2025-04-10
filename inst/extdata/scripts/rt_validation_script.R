### Load required packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(zoo)

### Generate Rt values that follow a sine wave

# Set seed
set.seed(123)

# Parameters
time <- seq(1, 730, by = 1) # Days in a year
amplitude <- 0.5           # Amplitude of the sine wave
cycles <- 4                # 4 cycles
frequency <- cycles / 365  # Cycles per year
phase <- 0                 # No phase shift
baseline <- 1.5            # Average R_t value
noise_sd <- 0.1            # Standard deviation of noise

# Generate sine wave
rt <- baseline + amplitude * sin(2 * pi * frequency * time + phase)

# Add noise
rt_noisy <- rt + rnorm(length(rt), mean = 0, sd = noise_sd)
rt_noisy <- pmax(rt_noisy, 0)  # Ensure R_t values are positive

# use EpiLPS to create a fake dataset
library(EpiLPS)

# set seed
set.seed(1234)

# specify serial interval distribution
si_spec <- Idist(mean = 8, sd = 5.66, dist = "gamma")
plot(si_spec)

# generate data
datasim <- episim(si = si_spec$pvec, Rpattern = 5, endepi = 40, dist = "negbin", overdisp = 15, plotsim = TRUE)

# Plot
t_vec <- seq(0, length(datasim$y)-1, by = 1)
rt_true <- datasim$Rtrue(t_vec)
data <- data.frame(time = t_vec, inc = datasim$y, Rt_True = rt_true) %>%
  mutate(onset_date = as.Date("2024-01-01") + time)

ggplot(data, aes(x = time, y = Rt_True)) +
  geom_line() +
  labs(x = "Days", y = "R") +
  theme_minimal()

# use mitey::rt_estim()
library(devtools)
load_all()

rt_estimated <- rt_estim(
  inc_dat = data,
  mean_si = 8,
  sd_si = 5.66,
  dist_si = "gamma"
  )

data$Rt_Estimated <- rt_estimated$rt_adjusted

data %>%
  select(time, onset_date, inc, Rt_True, Rt_Estimated) %>%
  pivot_longer(Rt_True:Rt_Estimated) %>%
  rename(Quantity = name) %>%
  ggplot(aes(x = time, y = value, color = Quantity)) +
  geom_line() +
  labs(x = "Days", y = "R") +
  theme_minimal()

# repeat with bootstrap to get confidence bounds
rt_boot <- rt_estim_w_boot(
  inc_dat = data,
  mean_si = 8,
  sd_si = 5.66,
  dist_si = "gamma",
  n_bootstrap = 100
)

data$Rt_Lower <- rt_boot$results$lower_rt_adjusted
data$Rt_Upper <- rt_boot$results$upper_rt_adjusted

# smooth Rt_Estimated
data_smoothed <- data %>%
  filter(!is.infinite(Rt_Estimated),
         !is.na(Rt_Estimated),
         !is.nan(Rt_Estimated),
         Rt_Estimated < 20) %>%
  mutate(
    Rt_Estimated_Smoothed = rollmean(Rt_Estimated, k = 7, fill = NA, align = "center"),
    Rt_Lower_Smoothed = rollmean(Rt_Lower, k = 7, fill = NA, align = "center"),
    Rt_Upper_Smoothed = rollmean(Rt_Upper, k = 7, fill = NA, align = "center"),
    )

# shift Rt_Estimated
data_shifted <- data %>%
  mutate(Rt_Estimated_Shifted = lag(Rt_Estimated, n = 8),
         Rt_Lower_Shifted = lag(Rt_Lower, n = 8),
         Rt_Upper_Shifted = lag(Rt_Upper, n = 8)
         )

data_smoothed %>%
  ggplot(aes(x = time)) +
  # Add confidence bounds as ribbon
  geom_ribbon(aes(ymin = Rt_Lower_Smoothed, ymax = Rt_Upper_Smoothed),
              fill = "blue", alpha = 0.2) +
  # Add shifted Rt_Estimated line
  # geom_line(aes(y = Rt_Estimated, color = "Rt_Estimated"), linewidth = 1, linetype = "solid") +
  # Add original Rt_Estimated line
  geom_line(aes(y = Rt_Estimated_Smoothed, color = "Rt_Estimated_Smoothed"), linewidth = 1) +
  # Add true Rt line
  geom_line(aes(y = Rt_True, color = "Rt_True"), linewidth = 1) +
  # Customize labels and theme
  labs(x = "Days", y = "R", color = "Legend") +
  scale_color_manual(values = c("Rt_True" = "red", "Rt_Estimated_Smoothed" = "blue"
                                # , "Rt_Estimated" = "green"
                                )) +
  theme_minimal()

### Create symptom onset data using the above Rt values

# Function to simulate new cases based on Rt
simulate_new_cases <- function(rt_values,
                                 generation_time,
                                 initial_cases = 1,
                                 days = 180,
                                 population = 1e6,
                                 waning_rate = 0.001) {

  # Normalize generation time, so it sums to 1
  generation_time <- generation_time / sum(generation_time)

  # Initialize variables
  cases <- numeric(days)
  cumulative_infections <- 0  # Start with no infections
  cases[1] <- initial_cases
  susceptible_pool <- population

  # Track onset dates
  onset_dates <- rep(1, initial_cases)

  for (day in 2:days) {
    # Compute remaining susceptibles
    susceptible_pool <- susceptible_pool + waning_rate * (population - susceptible_pool)

    # Compute expected new cases using convolution with generation interval
    lambda <- sum(cases[max(1, day - length(generation_time)):(day - 1)] *
                    rev(generation_time[1:min(day - 1, length(generation_time))]))

    # Adjust Rt for susceptible depletion
    effective_rt <- rt_values[day] * (susceptible_pool / population)

    # Compute expected new cases
    expected_new_cases <- effective_rt * lambda

    # Simulate new cases
    new_cases <- rpois(1, expected_new_cases)

    # Update cases, cumulative infections, and onset dates
    cases[day] <- new_cases
    cumulative_infections <- cumulative_infections + new_cases

    # Update susceptible pool
    susceptible_pool <- susceptible_pool - new_cases
    susceptible_pool <- max(susceptible_pool, 0)  # Ensure non-negative

    if (new_cases > 0) {
      onset_dates <- c(onset_dates, rep(day, new_cases))
    }
  }

  # Return results
  list(cases = cases, onset_dates = onset_dates, cumulative_infections = cumulative_infections)
}
###

### Run simulation
# Define parameters
rt_values <- rt_noisy #1.5 + 0.5 * sin(2 * pi * (2 / 365) * seq(1, 365))  # Rt values (sine wave)
generation_time <- dgamma(1:30, shape = 2, scale = 2)  # Generation interval

# Run simulation
simulation_result <- simulate_new_cases(rt_values = rt_values,
                                        generation_time = generation_time,
                                        initial_cases = 25,
                                        days = 730, #365,
                                        waning_rate = 0.025)

# Extract cases and scaled Rt
cases <- simulation_result$cases
days <- seq_along(cases)  # Day numbers
scaled_rt <- rt_values * max(cases) / max(rt_values)

# Combine into a data frame
my_data <- data.frame(
  Day = days,
  New_Cases = cases,
  Scaled_Rt = scaled_rt
)

# plot input (scaled Rt) with cases
data_long <- my_data %>%
  pivot_longer(cols = c("New_Cases", "Scaled_Rt"),
               names_to = "Metric", values_to = "Value")

# Plot using ggplot
ggplot(data_long, aes(x = Day, y = Value, color = Metric, linetype = Metric)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("New_Cases" = "blue",
                                "Scaled_Rt" = "red")) +
  scale_linetype_manual(values = c("New_Cases" = "solid",
                                   "Scaled_Rt" = "solid")) +
  labs(
    x = "Day",
    y = "Value",
    color = "Legend",
    linetype = "Legend"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )


### Estimate Rt from simulated cases using mitey::rt_estim()
# Create incidence data frame for rt_estim
inc_dat <- data.frame(
  onset_date = seq_along(cases),
  inc = cases
) %>%
  mutate(onset_date = as.Date("2023-01-01") + (onset_date - 1))

# Estimate Rt using rt_estim
gi_mean <- sum(seq_along(generation_time) * generation_time)
gi_sd <- sqrt(sum(generation_time * (seq_along(generation_time) - gi_mean)^2))

rt_estimated <- rt_estim_w_boot(
  inc_dat = inc_dat,
  mean_si = gi_mean,
  sd_si = gi_sd,
  dist_si = "gamma",
  n_bootstrap = 100
)

# Merge the estimated Rt into the main data frame
my_data$Estimated_Rt <- rt_estimated$results$median_rt_adjusted
my_data$Estimated_Rt_Lower <- rt_estimated$results$lower_rt_adjusted
my_data$Estimated_Rt_Upper <- rt_estimated$results$upper_rt_adjusted
my_data$True_Rt <- rt_values

# Reshape data for ggplot
data_long <- my_data %>%
  pivot_longer(cols = c("New_Cases", "Scaled_Rt", "Estimated_Rt",
                        "Estimated_Rt_Lower", "Estimated_Rt_Upper", "True_Rt"),
               names_to = "Metric", values_to = "Value")

# Plot using ggplot
ggplot(data_long %>%
         filter(Metric %in% c("Estimated_Rt", "Estimated_Rt_Lower",
                              "Estimated_Rt_Upper", "True_Rt"),
                Day > 45, # remove first period before epidemic has sinusoidal pattern
                Day < max(Day) - (2 * gi_mean)),
       aes(x = Day, y = Value, color = Metric, linetype = Metric)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("Estimated_Rt" = "blue",
                                "Estimated_Rt_Lower" = "blue",
                                "Estimated_Rt_Upper" = "blue",
                                "True_Rt" = "red")) +
  scale_linetype_manual(values = c("Estimated_Rt" = "solid",
                                   "Estimated_Rt_Lower" = "dotted",
                                   "Estimated_Rt_Upper" = "dotted",
                                   "True_Rt" = "solid")) +
  labs(
    x = "Day",
    y = "Value",
    color = "Legend",
    linetype = "Legend"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

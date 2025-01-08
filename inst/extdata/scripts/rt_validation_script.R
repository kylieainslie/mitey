### Load required packages
library(ggplot2)
library(dplyr)

### Generate Rt values that follow a sine wave

# Parameters
time <- seq(0, 365, by = 1) # Days in a year
amplitude <- 0.5           # Amplitude of the sine wave
frequency <- 2 / 365       # Two cycles per year
phase <- 0                 # No phase shift
baseline <- 1.5            # Average R_t value
noise_sd <- 0.1            # Standard deviation of noise

# Generate sine wave
rt <- baseline + amplitude * sin(2 * pi * frequency * time + phase)

# Add noise
set.seed(123) # For reproducibility
rt_noisy <- rt + rnorm(length(rt), mean = 0, sd = noise_sd)

# Ensure R_t values are positive
rt_noisy <- pmax(rt_noisy, 0)

# Plot
data <- data.frame(time = time, R_t = rt_noisy)
ggplot(data, aes(x = time, y = R_t)) +
  geom_line() +
  labs(title = "Synthetic R_t Data", x = "Time", y = "R_t") +
  theme_minimal()


### Create symptom onset data using the above Rt values

# Define parameters
set.seed(123)
rt_values <- rt_noisy #1.5 + 0.5 * sin(2 * pi * (2 / 365) * seq(1, 365))  # Rt values (sine wave)
generation_time <- dgamma(1:30, shape = 2, scale = 2)  # Generation interval

# Function to simulate new cases based on Rt
simulate_onset_dates <- function(rt_values, generation_time, initial_cases = 5, days = 180) {
  cases <- numeric(days)
  cases[1] <- initial_cases  # Initial cases on day 1

  onset_dates <- rep(1, initial_cases)  # Initial symptom onset dates

  for (day in 2:days) {
    # Compute expected new cases using convolution with generation interval
    lambda <- sum(cases[max(1, day - length(generation_time)):(day - 1)] *
                    rev(generation_time[1:min(day - 1, length(generation_time))]))

    # Simulate new cases
    new_cases <- rpois(1, rt_values[day] * lambda)

    cases[day] <- new_cases
    if (new_cases > 0) {
      onset_dates <- c(onset_dates, rep(day, new_cases))  # Append new onset dates
    }
  }

  list(cases = cases, onset_dates = onset_dates)
}


# Run simulation
simulation_result <- simulate_onset_dates(rt_values, generation_time, days = 180)

# Convert onset dates to actual dates
onset_dates <- as.Date("2024-01-01") + simulation_result$onset_dates - 1

# Create a data frame for plotting
onset_df <- data.frame(onset_date = onset_dates)

# Plot the simulated epidemic curve
ggplot(onset_df, aes(x = onset_date)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Simulated Symptom Onset Dates", x = "Date", y = "Number of Cases") +
  theme_minimal()

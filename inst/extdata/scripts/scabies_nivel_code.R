# scabies

# estimate exponential growth rate from yearly NIVEL data
# https://www.nivel.nl/nl/resultaten-van-onderzoek/nivel-cijfers-ziekten-op-jaarbasis

# load required packages
library(readxl)
library(tidyverse)
library(broom)
library(nls2)

# read in data set
scabies_inc_total <- read_xlsx("./inst/extdata/data/scabies_data_yearly.xlsx",
                               sheet = "total") %>%
  rename(inc = `Inc per 1.000`) %>%
  mutate(cases = as.numeric(inc),
         time = Year) %>%
  select("time", "cases")

# Define the exponential growth model
exponential_model <- function(x, a, b) {
  a * exp(b * x)
}

# get starting values
# Fit linear regression to log-transformed data
lm_fit <- lm(log(cases) ~ time, data = scabies_inc_total)

# Extract coefficients
a_init <- exp(coef(lm_fit)[1])
b_init <- coef(lm_fit)[2]

# Fitting the model with initial values
fit <- nls(cases ~ exponential_model(time, a, b),
           data = scabies_inc_total,
           start = list(a = a_init, b = b_init),
           control = list(minFactor = 0.00001, maxiter = 100000))

# Generate future time points
future_time <- seq(max(scabies_inc_total$time), max(scabies_inc_total$time) + 10,
                   by = 1)  # Example: extrapolating 10 time points into the future

# Create a data frame with both original and future time points
all_time <- data.frame(time = c(scabies_inc_total$time, future_time))

# Predict cases for all time points with confidence intervals using bootstrapping
set.seed(123)  # For reproducibility
bootstraps <- 1000
bootstrap_results <- replicate(bootstraps, {
  boot_indices <- sample(1:nrow(scabies_inc_total), replace = TRUE)
  boot_data <- scabies_inc_total[boot_indices, ]
  boot_fit <- nls(cases ~ exponential_model(time, a, b),
                  data = boot_data,
                  start = list(a = a_init, b = b_init),
                  control = list(minFactor = 0.000001, maxiter = 100000))
  predict(boot_fit, newdata = all_time)
})

# Calculate mean and confidence intervals
predicted_mean <- apply(bootstrap_results, 1, mean)
predicted_lower <- apply(bootstrap_results, 1, quantile, probs = 0.025)
predicted_upper <- apply(bootstrap_results, 1, quantile, probs = 0.975)

# Combine into a data frame
predicted_df <- data.frame(
  time = all_time$time,
  cases = predicted_mean,
  lower = predicted_lower,
  upper = predicted_upper
)

# Combine original data with predictions
scabies_plot_df <- bind_rows(
  scabies_inc_total %>% mutate(type = "Observed"),
  predicted_df %>% mutate(type = "Predicted")
)

# Plot using ggplot2
ggplot(scabies_plot_df, aes(x = time, y = cases, color = type)) +
  geom_point(data = scabies_inc_total, aes(x = time, y = cases), size = 2) +
  geom_line(size = 1.2) +
  geom_ribbon(data = predicted_df, aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "blue") +
  labs(title = "Scabies Incidence with Exponential Growth Model",
       x = "Year", y = "Incidence per 1000",
       color = "Data Type") +
  theme_minimal()

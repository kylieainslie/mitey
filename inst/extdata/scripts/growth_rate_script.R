# scabies

# estimate exponential growth rate from yearly NIVEL data
# https://www.nivel.nl/nl/resultaten-van-onderzoek/nivel-cijfers-ziekten-op-jaarbasis

# load required packages
library(readxl)
library(tidyverse)
library(broom)
library(nls2)
library(viridis)

# read in data set
scabies_inc_total <- read_xlsx("./inst/extdata/data/scabies_data_yearly.xlsx",
                               sheet = "total") %>%
  rename(inc = `Inc per 1.000`) %>%
  mutate(cases = as.numeric(inc),
         time = Year) %>%
  select("time", "cases")

scabies_inc_total <- scabies_inc_total %>%
  filter(!is.na(cases) & !is.infinite(cases) & !is.na(time) & !is.infinite(time))

# Define the exponential growth model
exponential_model <- function(x, a, b) {
  a * exp(b * x)
}

### get starting values
# Center the time variable around its mean for numerical stability
scabies_inc_total <- scabies_inc_total %>%
  mutate(time_centered = time - mean(time))

# Log-transform the cases for linear regression
log_cases <- log(scabies_inc_total$cases)

# Perform linear regression on log-transformed data
lm_fit <- lm(log_cases ~ time_centered, data = scabies_inc_total)

# Extract initial values from linear regression
a_init <- exp(coef(lm_fit)[1])
b_init <- coef(lm_fit)[2]

### Fitting the model with initial values
fit <- nls(cases ~ exponential_model(time_centered, a, b),
           data = scabies_inc_total,
           start = list(a = a_init, b = b_init),
           control = list(minFactor = 0.00001, maxiter = 100000))

# Generate future time points
future_time <- seq(max(scabies_inc_total$time) + 1, max(scabies_inc_total$time) + 10, by = 1)
future_time_centered <- future_time - mean(scabies_inc_total$time)

# predict future points
predicted_values <- predict(fit, newdata = data.frame(time_centered = future_time_centered))

### Bootstrap

# Number of bootstrap samples
bootstraps <- 1000

# Initialize the counter and results storage
null_fit_count <- 0
successful_fit_count <- 0
total_iterations <- 0

# Initialize a matrix to store the bootstrap results
bootstrap_results <- matrix(NA, nrow = length(future_time), ncol = bootstraps)
boot_a <- c(rep(NA, bootstraps))
boot_b <- c(rep(NA, bootstraps))
# Bootstrap procedure
while (successful_fit_count < bootstraps) {
  total_iterations <- total_iterations + 1
  boot_indices <- sample(1:nrow(scabies_inc_total), replace = TRUE)
  boot_data <- scabies_inc_total[boot_indices, ]
  boot_fit <- tryCatch({
    nls(cases ~ exponential_model(time_centered, a, b),
        data = boot_data,
        start = list(a = a_init, b = b_init),
        control = list(minFactor = 0.00001, maxiter = 100000))
  }, error = function(e) NULL)

  if (is.null(boot_fit)) {
    null_fit_count <- null_fit_count + 1  # Increment the NULL fit counter
  } else {
    boot_a[successful_fit_count + 1] <- coef(boot_fit)[1]
    boot_b[successful_fit_count + 1] <- coef(boot_fit)[2]
    bootstrap_results[, (successful_fit_count + 1)] <- predict(boot_fit, newdata = data.frame(time_centered = future_time_centered))
    successful_fit_count <- successful_fit_count + 1  # Increment the successful fit counter
  }
}

# Calculate growth rate confidence intervals
#gr <- mean(boot_b, na.rm = TRUE)
gr_lower <- quantile(boot_b, probs = 0.025, na.rm = TRUE)
gr_upper <- quantile(boot_b, probs = 0.975, na.rm = TRUE)

# Calculate mean and confidence intervals for predictions
predicted_means <- apply(bootstrap_results, 1, mean, na.rm = TRUE)
predicted_lowers <- apply(bootstrap_results, 1, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
predicted_uppers <- apply(bootstrap_results, 1, function(x) quantile(x, probs = 0.975, na.rm = TRUE))

# Combine original data with predicted values
extrapolated_df <- data.frame(time = c(tail(scabies_inc_total$time, 1), future_time),
                              cases_mean = c(tail(scabies_inc_total$cases, 1), predicted_means),
                              cases_lower = c(tail(scabies_inc_total$cases, 1), predicted_lowers),
                              cases_upper = c(tail(scabies_inc_total$cases, 1), predicted_uppers))

# Combine original and extrapolated data for plotting
scabies_plot_df <- bind_rows(
  scabies_inc_total %>% mutate(type = "Original Data", cases_mean = cases, cases_lower = NA, cases_upper = NA),
  extrapolated_df %>% mutate(type = "Extrapolated Data")
)

# Plot the results with ggplot2
saveRDS(scabies_plot_df, file = "vignettes/scabies_gr_proj_df.rds")

# Define the custom palette
custom_palette <- c("#21908C", "#440154") # Custom colors

# Plot with custom palette and specified line and point styles
ggplot(scabies_plot_df, aes(x = time, y = cases_mean, color = type, fill = type)) +
  # Points with no line through them
  geom_point(data = filter(scabies_plot_df, type == "Original Data"), aes(y = cases), size = 2, shape = 16) +

  # Line with points for other types
  geom_line(data = filter(scabies_plot_df, type != "Original Data"), aes(y = cases_mean), linewidth = 1.2) +
  geom_point(data = filter(scabies_plot_df, type != "Original Data"), aes(y = cases), size = 2, shape = 16) +

  # Ribbon for the confidence interval
  geom_ribbon(aes(ymin = cases_lower, ymax = cases_upper, fill = type), alpha = 0.2, color = NA) +

  # Labels and theme
  labs(x = "Year", y = "Incidence per 1000",
       color = "Data Type", fill = "Data Type") +
  scale_color_manual(values = custom_palette) +
  scale_fill_manual(values = custom_palette) +
  theme_minimal() +
  theme(legend.position = "right")

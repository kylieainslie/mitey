# scabies

# estimate exponential growth rate from yearly NIVEL data
# https://www.nivel.nl/nl/resultaten-van-onderzoek/nivel-cijfers-ziekten-op-jaarbasis

# load required packages
library(readxl)
library(tidyverse)

# read in data set
setwd("~/Dropbox/Kylie/Projects/RIVM/Projects/scabies/data")
scabies_inc_total <- read_xlsx("scabies_data_yearly.xlsx", sheet = "total") %>%
  rename(inc = `Inc per 1.000`) %>%
  mutate(cases = as.numeric(inc),
         time = Year) %>%
  select("time", "cases")
#scabies_in_by_age_group <- read_xlsx("scabies_data.xlsx", sheet = "Scabies inc per 1.000 by age gr")

# Plot the data
plot(scabies_inc_total$time, scabies_inc_total$cases,
     xlab = "Time", ylab = "Incidence", main = "Exponential Growth Data")


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

# Summary of the model
summary(fit)

# Plot 
plot(scabies_inc_total$time, scabies_inc_total$cases, 
     ylim = c(0, max(scabies_inc_total$cases) * 1.1), 
     xlim = c(min(scabies_inc_total$time), max(scabies_inc_total$time)), 
     xlab = "Year", ylab = "Cases per 1000", main = "Scabies Exponential Growth Model Fit")
lines(scabies_inc_total$time, predict(fit), col = "red")
legend("topleft", legend = c("Data", "Model Fit"), col = c("black", "red"), lty = 1)

# Generate future time points
future_time <- seq(max(scabies_inc_total$time), max(scabies_inc_total$time) + 10, by = 1)  # Example: extrapolating 10 time points into the future

# Predict future values
predicted_values <- predict(fit, newdata = data.frame(time = future_time))

# Combine original data with predicted values
extrapolated_df <- data.frame(time = future_time, cases = predicted_values)

# plot fitted data with extrapolated data
plot(scabies_inc_total$time, scabies_inc_total$cases, 
     ylim = c(0, max(scabies_inc_total$cases, extrapolated_df$cases) * 1.1), 
     xlim = c(min(scabies_inc_total$time), max(extrapolated_df$time)), 
     xlab = "Year", ylab = "Incidence per 1000", 
     main = "Scabies Incidence")
lines(scabies_inc_total$time, predict(fit), col = "blue", lwd = 2) # Fitted line
lines(extrapolated_df$time, extrapolated_df$cases, col = "red", lwd = 2, lty = 2) # Extrapolated line
legend("topleft", legend = c("Original Data", "Extrapolated Data", "Model Fit"), 
       col = c("black", "red", "blue"), lty = c(NA, 2, 1), lwd = 2, pch = c(1, NA, NA))

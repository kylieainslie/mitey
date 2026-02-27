## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7, 
  fig.height = 5,
  dpi = 300,
  warning = FALSE,
  message = FALSE
)


## ----setup--------------------------------------------------------------------
library(mitey)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(cowplot) 
library(here)
library(outbreaks)


## ----simulate-si-data---------------------------------------------------------
set.seed(1234)

# Parameters for simulation
N <- 500           # Number of observations
true_mean <- 15    # True mean serial interval (days)
true_sd <- 3       # True standard deviation (days)
route_weights <- c(0.2, 0.5, 0.2, 0.1)  # Weights for transmission routes

# Generate data for different transmission routes
CP <- fdrtool::rhalfnorm((route_weights[1]*N), theta=sqrt(pi/2)/(sqrt(2)*true_sd))  # Co-Primary
PS <- rnorm(route_weights[2]*N, mean=true_mean, sd=true_sd)                # Primary-Secondary
PT <- rnorm(route_weights[3]*N, mean=2*true_mean, sd=sqrt(2)*true_sd)      # Primary-Tertiary
PQ <- rnorm(route_weights[4]*N, mean=3*true_mean, sd=sqrt(3)*true_sd)      # Primary-Quaternary

# Combine and round to days
sim_icc_intervals <- round(c(CP, PS, PT, PQ))

# Visualize the simulated data
hist(sim_icc_intervals, 
     breaks = seq(min(sim_icc_intervals)-0.5, max(sim_icc_intervals)+0.5, by=1),
     main = "Simulated ICC Intervals", 
     xlab = "Days since index case onset",
     col = "lightblue")


## ----si-estimate-sim----------------------------------------------------------
# Estimate serial interval assuming Normal distribution
si_results <- si_estim(sim_icc_intervals, dist = "normal")
si_results


## ----si_estimate-sim-res, echo=FALSE------------------------------------------
# Display results
cat("True parameters:\n")
cat("Mean:", true_mean, "days\n")
cat("SD:", true_sd, "days\n\n")

cat("Estimated parameters:\n")
cat("Mean:", round(si_results$mean[1], 2), "days\n")
cat("SD:", round(si_results$sd[1], 2), "days\n")


## ----plot-sim-fit, fig.height=5, fig.width=7----------------------------------
# Extract weights
weights <- c(si_results$wts[1], 
            si_results$wts[2] + si_results$wts[3],
            si_results$wts[4] + si_results$wts[5], 
            si_results$wts[6] + si_results$wts[7])

# Plot the fitted distribution
plot_si_fit(
  dat = sim_icc_intervals,
  mean = si_results$mean[1],
  sd = si_results$sd[1],
  weights = weights,
  dist = "normal"
) +
  ggtitle("Fitted Serial Interval Distribution (Simulated Data)") +
  theme(plot.title = element_text(hjust = 0.5))


## ----convergence-diagnostics--------------------------------------------------
# Check convergence information from our previous result
cat("Converged:", si_results$converged, "\n")
cat("Iterations used:", si_results$iterations, "out of 50 (default max)\n")
cat("Log-likelihood:", round(si_results$loglik, 2), "\n")


## ----convergence-comparison---------------------------------------------------
# Compare with no early stopping (tol = 0)
si_no_early_stop <- si_estim(sim_icc_intervals, dist = "normal", n = 50, tol = 0)

cat("With early stopping: ", si_results$iterations, " iterations\n")
cat("Without early stopping:", si_no_early_stop$iterations, " iterations\n")
cat("Results are identical:",
    all.equal(si_results$mean, si_no_early_stop$mean, tolerance = 1e-5), "\n")


## ----multiple-restarts--------------------------------------------------------
set.seed(456)

# Single restart with poor initial values
result_single <- si_estim(sim_icc_intervals, init = c(30, 8), n_starts = 1)

# Multiple restarts - algorithm explores parameter space more thoroughly
result_multi <- si_estim(sim_icc_intervals, init = c(30, 8), n_starts = 5)

cat("Single restart:\n")
cat("  Mean:", round(result_single$mean, 2), ", SD:", round(result_single$sd, 2), "\n")
cat("  Log-likelihood:", round(result_single$loglik, 2), "\n\n")

cat("Multiple restarts (n_starts = 5):\n")
cat("  Mean:", round(result_multi$mean, 2), ", SD:", round(result_multi$sd, 2), "\n")
cat("  Log-likelihood:", round(result_multi$loglik, 2), "\n")
cat("  (Best result selected from", result_multi$n_restarts, "restarts)\n")


## ----load-real-si-data--------------------------------------------------------
# Load scabies ICC interval data
file_path <-  here("vignettes", "data", "si_data.rds")
scabies_si_data <- readRDS(file_path)


## ----si-estimate-real---------------------------------------------------------
# Estimate serial interval for each study
result_by_study <- scabies_si_data %>%
  group_by(study) %>%
  summarise(result = list(si_estim(icc_interval))) %>%
  mutate(
    mean = map_dbl(result, "mean"),
    sd = map_dbl(result, "sd"),
    wts = map(result, "wts")
  ) %>%
  select(-result)

# Display results
result_by_study %>%
  select(study, mean, sd) %>%
  mutate(across(c(mean, sd), round, 2)) %>%
  arrange(mean) %>%
  knitr::kable(caption = "Estimated mean and standard deviation of serial interval (days) by study")


## ----process-si-results-------------------------------------------------------
# Process weights for plotting
result_wide <- result_by_study %>%
  unnest(wts) %>%
  pivot_longer(
    cols = c(mean, sd, wts),
    names_to = "statistic",
    values_to = "value"
  ) %>%
  group_by(study, statistic) %>%
  mutate(
    occurrence = row_number(),
    statistic = if_else(statistic == "wts", paste0("weight_", occurrence), statistic)
  ) %>%
  filter(statistic != "mean" | occurrence == 1) %>%
  filter(statistic != "sd" | occurrence == 1) %>%
  select(-occurrence) %>%
  ungroup() %>%
  pivot_wider(
    names_from = statistic,
    values_from = value
  )

# Merge with original data for plotting
df_merged <- scabies_si_data %>%
  left_join(result_wide, by = "study", relationship = "many-to-many")


## ----plot-si-multi, fig.height=8, fig.width=10--------------------------------
# Create a function to generate plot for each study
plot_study <- function(study_data) {
  study_name <- unique(study_data$study)
  
  plot_si_fit(
    dat = study_data$icc_interval,
    mean = study_data$mean[1],
    sd = study_data$sd[1],
    weights = c(study_data$weight_1[1], 
                study_data$weight_2[1] + study_data$weight_3[1],
                study_data$weight_4[1] + study_data$weight_5[1], 
                study_data$weight_6[1] + study_data$weight_7[1]),
    dist = "normal",
    scaling_factor = 0.25
  ) +
    ggtitle(study_name) +
    theme(plot.title = element_text(hjust = 0.5, size = 11))
}

# Generate plots for each study
study_plots <- df_merged %>%
  group_by(study) %>%
  group_split() %>%
  map(plot_study)

# Combine plots
combined_plot <- plot_grid(
  plotlist = study_plots,
  labels = "AUTO",
  ncol = 2
)

# Display combined plot
combined_plot


## ----func_arguments, eval=FALSE-----------------------------------------------
# wallinga_lipsitch(
#   incidence,          # Vector of case counts
#   dates,              # Vector of dates corresponding to incidence
#   si_mean,            # Mean of serial interval distribution (days)
#   si_sd,              # Standard deviation of serial interval
#   si_dist = "gamma",  # Distribution type ("gamma" or "normal")
#   smoothing = 0,      # Window size for smoothing (0 for no smoothing)
#   bootstrap = FALSE,  # Whether to compute bootstrap CIs
#   n_bootstrap = 1000, # Number of bootstrap samples
#   conf_level = 0.95,  # Confidence level for intervals
#   shift = FALSE       # Whether to shift estimates by one serial interval
# )


## ----sim-rt-------------------------------------------------------------------
set.seed(42)

# Simulation parameters
t_end <- 100        # Simulation duration (days)
si_mean <- 7        # Mean serial interval (days)
si_sd <- 2          # SD of serial interval (days)

# Define a time-varying reproduction number function
# This creates a pattern where Rt starts high, decreases below 1, then increases again
true_rt <- function(t) {
  if(t < 20) return(2.5)                  # Initial high Rt
  if(t < 40) return(2.5 - 0.075 * (t - 20)) # Linear decrease
  if(t < 60) return(1.0)                  # Stable period at Rt=1
  if(t < 80) return(1.0 + 0.05 * (t - 60))  # Linear increase
  return(2.0)                             # Final high Rt
}

# Create vector of true Rt values for plotting
true_rt_values <- sapply(1:t_end, true_rt)

# Initialize with some seed cases
cases <- numeric(t_end)
cases[1:5] <- c(1, 2, 3, 5, 8)

# Serial interval distribution (discretized normal)
si_pmf <- dnorm(0:30, mean = si_mean, sd = si_sd)
si_pmf <- si_pmf / sum(si_pmf)  # Normalize to sum to 1

# Generate incidence using renewal equation model
for(t in 6:t_end) {
  # Calculate expected new cases
  lambda <- 0
  for(s in 1:min(t-1, length(si_pmf))) {
    lambda <- lambda + cases[t-s] * true_rt(t-s) * si_pmf[s]
  }
  
  # Add randomness (negative binomial distribution)
  cases[t] <- rnbinom(1, mu = lambda, size = 10)
}

# Create dates sequence
sim_dates <- seq.Date(as.Date("2023-01-01"), by = "day", length.out = t_end)

# Create data frame
sim_epidemic <- data.frame(
  date = sim_dates,
  cases = cases,
  true_rt = true_rt_values
)

# Plot the simulated epidemic curve with true Rt
p1 <- ggplot(sim_epidemic, aes(x = date, y = cases)) +
  geom_col(fill = "steelblue") +
  labs(
    x = "Date",
    y = "Daily Cases",
    title = "Simulated Epidemic Curve"
  ) +
  theme_minimal()

p2 <- ggplot(sim_epidemic, aes(x = date, y = true_rt)) +
  geom_line(color = "red", linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  labs(
    x = "Date",
    y = "Reproduction Number",
    title = "True Time-Varying Rt"
  ) +
  ylim(0, 3) +
  theme_minimal()

# Combine plots
plot_grid(p1, p2, ncol = 1)


## ----estimate-rt--------------------------------------------------------------
# Estimate Rt using wallinga_lipsitch
rt_estimates <- wallinga_lipsitch(
  incidence = sim_epidemic$cases,
  dates = sim_epidemic$date,
  si_mean = si_mean,
  si_sd = si_sd,
  si_dist = "normal",
  smoothing = 7,        # 7-day smoothing window
  bootstrap = TRUE,
  n_bootstrap = 100     # Use more in practice
)

# Convert to data frame for plotting
rt_est_df <- as.data.frame(rt_estimates)

# Compare estimated vs true Rt
ggplot() +
  # True Rt
  geom_line(
    data = sim_epidemic, 
    aes(x = date, y = true_rt, color = "True Rt"),
    linewidth = 1
  ) +
  
  # Estimated Rt (with right-truncation correction)
  geom_line(
    data = rt_est_df, 
    aes(x = date, y = R_corrected, color = "Estimated Rt"),
    linewidth = 1
  ) +
  
  # Confidence intervals
  geom_ribbon(
    data = rt_est_df,
    aes(x = date, ymin = R_corrected_lower, ymax = R_corrected_upper),
    fill = "blue", alpha = 0.2
  ) +
  
  # R=1 threshold
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  
  # Styling
  scale_color_manual(
    values = c("True Rt" = "red", "Estimated Rt" = "blue"),
    name = ""
  ) +
  labs(
    x = "Date",
    y = "Reproduction Number (Rt)",
    title = "Comparison of True vs. Estimated Rt",
    subtitle = "Using Wallinga-Lipsitch method with 7-day smoothing"
  ) +
  ylim(0, 3) +
  theme_minimal() +
  theme(legend.position = "bottom")


## ----rt-shift-----------------------------------------------------------------
# Re-estimate with shift=TRUE
rt_shifted <- wallinga_lipsitch(
  incidence = sim_epidemic$cases,
  dates = sim_epidemic$date,
  si_mean = si_mean,
  si_sd = si_sd,
  si_dist = "normal",
  smoothing = 7,
  bootstrap = FALSE,  # Skip bootstrap for speed
  shift = TRUE        # Use shift parameter
)

# Compare true Rt with both original and shifted estimates
ggplot() +
  # True Rt
  geom_line(
    data = sim_epidemic, 
    aes(x = date, y = true_rt, color = "True Rt"),
    linewidth = 1
  ) +
  
  # Original estimated Rt
  geom_line(
    data = rt_est_df, 
    aes(x = date, y = R_corrected, color = "Estimated Rt (original)"),
    linewidth = 1
  ) +
  
  # Shifted estimated Rt
  geom_line(
    data = as.data.frame(rt_shifted), 
    aes(x = shifted_date, y = R_corrected, color = "Estimated Rt (shifted)"),
    linewidth = 1
  ) +
  
  # R=1 threshold
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  
  # Styling
  scale_color_manual(
    values = c("True Rt" = "red", 
               "Estimated Rt (original)" = "blue",
               "Estimated Rt (shifted)" = "green"),
    name = ""
  ) +
  labs(
    x = "Date",
    y = "Reproduction Number (Rt)",
    title = "Effect of the 'shift' parameter on Rt estimates",
    subtitle = "Shifting forward by one serial interval (7 days)"
  ) +
  ylim(0, 3) +
  theme_minimal() +
  theme(legend.position = "bottom")


## ----get_zika_data------------------------------------------------------------
data(zika_girardot_2015)

# Examine the data structure
str(zika_girardot_2015)


## ----zike_epicurve, echo=FALSE------------------------------------------------
zika_epicurve <- ggplot(zika_girardot_2015, aes(x = date, y = cases)) +
  geom_col(fill = "steelblue") +
  scale_x_date(date_breaks = "14 days", date_labels = "%d-%m-%Y") +
  labs(
    x = "Date of onset",
    y = "Number of cases",
    title = "Daily Zika virus cases in Girardot, Colombia (2015-2016)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

zika_epicurve


## ----rt-estimate-zika---------------------------------------------------------
# Parameters for Zika
zika_si_mean <- 7  # days
zika_si_sd <- 1.5     # days

# Estimate Rt using the Wallinga-Lipsitch method
rt_zika <- wallinga_lipsitch(
  incidence = zika_girardot_2015$cases,
  dates = zika_girardot_2015$date,
  si_mean = zika_si_mean,
  si_sd = zika_si_sd,
  si_dist = "gamma",
  smoothing = 0,
  bootstrap = TRUE,
  n_bootstrap = 100,
  conf_level = 0.95
)

head(rt_zika)


## ----rt_viz-------------------------------------------------------------------
# Prepare data for visualization
rt_plot_data <- rt_zika %>%
  filter(!is.na(R_corrected)) %>%
  # Skip the first 7 days of unstable estimates
  filter(date > min(date) + 7)

# Plot Rt over time
rt_plot <- ggplot(rt_plot_data, aes(x = date)) +
  geom_ribbon(aes(ymin = R_corrected_lower, ymax = R_corrected_upper), 
              fill = "#21908C", alpha = 0.2) +
  geom_line(aes(y = R_corrected), color = "#21908C", size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  #coord_cartesian(ylim = c(0, 3)) +
  scale_x_date(date_breaks = "14 days", date_labels = "%d-%m-%Y") +
  labs(
    x = "Time", 
    y = "Reproduction number (Rt)",
    title = "Estimated reproduction number for Zika in Girardot, Colombia",
    subtitle = paste0("Serial interval: ", zika_si_mean, " days (SD: ", zika_si_sd, " days)")
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Show Rt plot with epicurve
plot_grid(zika_epicurve, rt_plot, ncol = 1, rel_heights = c(1, 1.5))


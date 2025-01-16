# script to validate Rt estimation using code from Gostic et al.

# To generate the synthetic data we borrow code from [Gostic et al.](https://github.com/cobeylab/Rt_estimation) [@gostic2020]. Briefly, the synthetic data were generated from a deterministic SEIR model in which the transmission rate changes abruptly.

# load packages
library(devtools)
load_all()
library(ggplot2)
library(tidyr)
library(dplyr)
library(viridis)


# set seed for reproducibility
set.seed(2984)

# define helper functions
# Define beta(t) from SEIR simulation
rt <- function(time) {
  baseline + amplitude * sin(2 * pi * frequency * time + phase)
}

# arnaught <- function(time) {
#   c(rt(time)[1], rt(time))
# }

beta_t <- function(time) {
  rt(time) * gamma
}


# simulate data
n <- 1e6 # total population size
n_I <- 10
n_t <- 365

# Parameters for Rt
time <- seq(0, n_t, by = 1)    # time_vector
amplitude <- 0.5               # Amplitude of the sine wave
cycles <- 4                    # Number of cycles
frequency <- cycles / 365      # Cycles per year
phase <- 0                     # No phase shift
baseline <- 1.5                # Average R_t value (R_0)

# Input parameters
rt_input <- rt(time)
t_E <- 4
t_I <- 4

sim_data <- simulate_seir_ode(
  arnaught = rt_input,
  t_E = t_E,
  t_I = t_I,
  N = n,
  S_init = n-n_I,
  E_init = 0,
  I_init = n_I,
  n_t = 365,
  n_steps_per_t = 1
)

# show sim data
head(sim_data)

#
# Figure 1. Dynamics of an SEIR model showing the number of individuals in each state (Susceptible, Exposed, Infectious, and Recovered) over time.
sim_data %>%
  pivot_longer(S:dIR) %>%
  rename(State = name) %>%
  filter(State %in% c("S", "E", "I", "R")) %>%
  mutate(State = factor(State, levels = c("S", "E", "I", "R"))) %>%
  ggplot()+
  geom_line(aes(x = time, y = value, color = State)) +
  labs(x = "Time (days)", y = "Number of Individuals") +
  #facet_wrap(.~name, scales = 'free_y') +
  #xlim(c(0, 200))
  theme_minimal()

# Figure 2. Plot of time-varying reproduction number.
ggplot(sim_data, aes(x = time, y = Rt)) +
  geom_line(color = "blue", linewidth = 1) +
  labs(
    title = "Effective Reproduction Number (R_t) Over Time",
    x = "Time (days)",
    y = expression(R[t])
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5)
  )

# Calculate true Rt
#beta_t <- approxfun(sim_data$time, c(r0,rt) / t_I)

# arnaught <- 2.0 #rt_noisy
# beta_t <- function(t) { arnaught / t_I }
#  function(t) {
#   approx(1:length(arnaught), arnaught, xout = t)$y  # Interpolates Rt values
# }

# Define parameters for integrate_Rt
sigma <- 1 / t_E
gamma <- 1 / t_I

# Calculate true case reproduction number
Rt_output <- integrate_Rt(
  beta_t = beta_t,
  sigma = sigma,
  gamma = gamma,
  N = 1,
  T_final = 365,
  E0 = 0.001,
  I0 = 0)

# Extract and plot results
case_rt <- Rt_output$case_rt_df
case_rt$R_instantaneous <- beta_t(sim_data$time) * (sim_data$S/n)

plot(case_rt$time, case_rt$R_case, type = "l", xlab = "Time", ylab = "R", ylim = c(0,4))
lines(case_rt$time, case_rt$R_instantaneous, col = "blue")

## Calculate Rt with WL method using mitey and symptom onset data from simulated data
df_sim <- sim_data %>%
  select(time, dEI) %>%
  rename(onset_date = time,
         inc = dEI) %>%
  mutate(inc = if_else(is.na(inc), 0, inc))

rt_estimates <- rt_estim(df_sim, mean_si = 8, sd_si = 5.66, dist_si = "gamma")
rt_boot <- rt_estim_w_boot(df_sim, mean_si = 8, sd_si = 5.66, dist_si = "gamma", n_bootstrap = 100)

case_rt$R_WL <- rt_boot$results$median_rt_adjusted
case_rt$R_WL_lower <- rt_boot$results$lower_rt_adjusted
case_rt$R_WL_upper <- rt_boot$results$upper_rt_adjusted

# Compare estimates
case_rt %>%
  filter(time > t_E + t_I) %>%
  ggplot(aes(x = time)) +
  # Add confidence bounds for R_WL
  geom_ribbon(aes(ymin = R_WL_lower, ymax = R_WL_upper),
              fill = viridis(1, option = "viridis"), alpha = 0.2) +
  # Add the R_WL line
  geom_line(aes(y = R_WL, color = "R_WL"), linewidth = 1) +
  # Add the R_case line
  geom_line(aes(y = R_case, color = "R_case"), linewidth = 1) +
  # Customize the plot
  labs(x = "Time", y = "R", color = "Legend") +
  scale_color_viridis_d(option = "viridis", begin = 0.2, end = 0.8) +
  coord_cartesian(ylim = c(0, 5)) +
  theme_minimal()

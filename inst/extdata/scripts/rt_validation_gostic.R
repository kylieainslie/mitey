# script to validate Rt estimation using code from Gostic et al.

# To generate the synthetic data we borrow code from [Gostic et al.](https://github.com/cobeylab/Rt_estimation) [@gostic2020]. Briefly, the synthetic data were generated from a deterministic SEIR model in which the transmission rate changes abruptly.

# simulate data
n <- 1e6 # total population size
n_I <- 10

sim_data <- simulate_seir_ode(
  arnaught = 2.0,
  t_E = 4,
  t_I = 4,
  N = n,
  S_init = n-n_I,
  E_init = 0,
  I_init = n_I,
  n_t = 365,
  n_steps_per_t = 1
)

# show sim data
head(sim_data)

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

# Rt comparison
arnaught <- 2.0 #rt_noisy
beta_t <- function(t) { arnaught / t_I }
#  function(t) {
#   approx(1:length(arnaught), arnaught, xout = t)$y  # Interpolates Rt values
# }
t_E <- 4
t_I <- 4

## Calculate true case reproduction number
Rt_output <- integrate_Rt(
  beta_t = beta_t,
  sigma = 1/t_E,
  gamma = 1/t_I,
  N = 1,
  T_final = 365,
  E0 = 0.001,
  I0 = 0)

## Calculate Rt with WT method using mitey and symptom onset data from simulated data
df_sim <- sim_data %>%
  select(time, dEI) %>%
  rename(onset_date = time,
         inc = dEI)

rt_estimates <- rt_estim(df_sim[-1,], mean_si = 8, sd_si = 5.66, dist_si = "gamma")

# Compare estimates
Rt_output$case_rt_df[-1,] %>%
  mutate(R_WL = rt_estimates$rt_adjusted) %>%
  pivot_longer(R_case:R_WL) %>%
  rename(Method = name) %>%
  filter(time > 16) %>%
  ggplot() +
  geom_line(aes(x = time, y = value, color = Method)) +
  labs(x = "Time (days)", y = "Case Reproduction Number") +
  theme_minimal()

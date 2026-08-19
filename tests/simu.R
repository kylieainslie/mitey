rm(list= ls())

library(dplyr)
library(devtools)
library(fdrtool)


pkg_root <- r"(C:\Users\Administrator.DESKTOP-JQG2K7N\桌面\2025 melb\master_project\mitey\mitey)"
setwd(pkg_root)

load_all()
his_dat_path <- r"(C:\Users\Administrator.DESKTOP-JQG2K7N\桌面\2025 melb\master_project\mitey\mitey\inst\extdata\data\si_data.rds)"

#####################################
# Simulation using historical dataset
si_data <- readRDS(his_dat_path)

# Extract normal model results 
normal_results <- data.frame(
  study = c("Akunzirwe et al.", "Ariza et al.", "Kaburi et al.", "Tjon-Kon-Fat et al"),
  normal_mean = c(122.9, 98.4, 167.3, 110.7),
  normal_sd = c(26.9, 8.5, 9.7, 16.1)
)

# Use lapply to batch process the four studies and run the lognormal model
lognormal_results_list <- lapply(normal_results$study, function(study_name) {
  
  # Extract data for the current study
  dat <- si_data %>% 
    filter(study == study_name) %>% 
    pull(icc_interval)
  
  # Run the lognormal model
  fit <- si_estim(dat = dat, dist = "lognormal")
  
  data.frame(
    study = study_name,
    n_sample = length(dat),
    zero_icc = sum(dat == 0),
    lognormal_mean = fit$mean,
    lognormal_sd = fit$sd,
    converged = fit$converged,
    loglik = fit$loglik
  )
})

lognormal_results <- bind_rows(lognormal_results_list)

#normal model
dat <- si_data %>% 
    filter(study == "Akunzirwe et al.") %>% 
    pull(icc_interval)

fit <- si_estim(dat = dat, dist = "normal")




# Merge lognormal and normal results for comparison
comparison_table <- normal_results %>%
  left_join(lognormal_results, by = "study") %>%
  select(
    study, 
    normal_mean, lognormal_mean,  
    normal_sd, lognormal_sd, 
  )

# Print the final comparison table
print(comparison_table)


##########################################
# Simulation using synthetic dataset
rlnorm_natural <- function(n, mean, sd) {
  sdlog <- sqrt(log(1 + (sd / mean)^2))
  meanlog <- log(mean) - 0.5 * sdlog^2
  rlnorm(n, meanlog, sdlog)
}

set.seed(123)

# Parameters for simulation
N <- 200           
true_mean <- 15    # True mean serial interval (days)
true_sd <- 3      
route_weights <- c(0.2, 0.5, 0.2, 0.1)  

# Generate data for different transmission routes 
generate_sim_data <- function(N, true_mean, true_sd, route_weights) {
  CP <- rhalfnorm(route_weights[1] * N, theta = sqrt(pi / 2) / (sqrt(2) * true_sd))
  PS <- rlnorm_natural(route_weights[2] * N, mean = true_mean, sd = true_sd)
  PT <- rlnorm_natural(route_weights[3] * N, mean = true_mean, sd = true_sd) +
        rlnorm_natural(route_weights[3] * N, mean = true_mean, sd = true_sd)
  PQ <- rlnorm_natural(route_weights[4] * N, mean = true_mean, sd = true_sd) +
        rlnorm_natural(route_weights[4] * N, mean = true_mean, sd = true_sd) +
        rlnorm_natural(route_weights[4] * N, mean = true_mean, sd = true_sd)
  
  return(round(c(CP, PS, PT, PQ)))
}
sim_icc_intervals <- generate_sim_data(N, true_mean, true_sd, route_weights)

# Fit all three distributions and collect results
dists <- c("lognormal", "normal", "gamma")
si_results <- lapply(dists, function(d) {
  si_estim(sim_icc_intervals, dist = d, init = c(true_mean, true_sd))
})
names(si_results) <- dists

# Summary table comparing estimates to true values
do.call(rbind, lapply(dists, function(d) {
  fit <- si_results[[d]]
  data.frame(
    dist       = d,
    true_mean  = true_mean,
    est_mean   = fit$mean,
    true_sd    = true_sd,
    est_sd     = fit$sd,
    converged  = fit$converged,
    iterations = fit$iterations
  )
}))

# Plot fitted curves (estimated mean/sd/weights) against the data histogram
for (d in dists) {
  print(plot_si_fit_result(si_results[[d]], sim_icc_intervals, dist = d))
}

# ---------------------------------------------------------------------
# Compare change in estiamtes between weighted and unweighted models

# source unweightd version for comparison
source(r"(C:\Users\Administrator.DESKTOP-JQG2K7N\桌面\2025 melb\master_project\mitey\mitey\R\old\si_estim_unweighted.R)")

n_reps <- 50
# Run batch simulation
sim_results_list <- lapply(1:n_reps, function(i) {
  
  # Generate new synthetic data for each iteration
  sim_icc_intervals <- generate_sim_data(N, true_mean, true_sd, route_weights)
  
  # Fit both models separately
  fit_w <- si_estim(sim_icc_intervals, dist = "lognormal", init = c(true_mean, true_sd))
  fit_u <- si_estim_unweighted(sim_icc_intervals, dist = "lognormal", init = c(true_mean, true_sd))
  
  data.frame(
    rep = i,
    model = c("weighted", "unweighted"),
    est_mean = c(fit_w$mean, fit_u$mean),
    est_sd = c(fit_w$sd, fit_u$sd),
    converged = c(fit_w$converged, fit_u$converged)
  )
})

# Combine all repetition results
sim_results <- bind_rows(sim_results_list)

# 3. Calculate evaluation metrics
mse_comparison <- sim_results %>%
  # Evaluate only successfully converged results
  filter(converged == TRUE) %>%
  group_by(model) %>%
  summarise(
    convergence_rate = round(n() / n_reps, 3),
    mean_bias = mean(est_mean) - true_mean,
    mean_mse = mean((est_mean - true_mean)^2),
    .groups = 'drop'
  )

print(mse_comparison)

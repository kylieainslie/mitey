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
N <- 500           
true_mean <- 15    # True mean serial interval (days)
true_sd <- 3      
route_weights <- c(0.2, 0.5, 0.2, 0.1)  

# Generate data for different transmission routes 
CP <- rhalfnorm(route_weights[1] * N, theta = sqrt(pi / 2) / (sqrt(2) * true_sd))
PS <- rlnorm_natural(route_weights[2] * N, mean = true_mean, sd = true_sd)
PT <- rlnorm_natural(route_weights[3] * N, mean = true_mean, sd = true_sd) +
  rlnorm_natural(route_weights[3] * N, mean = true_mean, sd = true_sd)
PQ <- rlnorm_natural(route_weights[4] * N, mean = true_mean, sd = true_sd) +
  rlnorm_natural(route_weights[4] * N, mean = true_mean, sd = true_sd) +
  rlnorm_natural(route_weights[4] * N, mean = true_mean, sd = true_sd)

sim_icc_intervals <- round(c(CP, PS, PT, PQ))

si_results <- si_estim(
  sim_icc_intervals,
  dist = "lognormal",
  init = c(true_mean, true_sd)
)

# results
data.frame(
  true_mean  = true_mean,
  est_mean   = si_results$mean,
  true_sd    = true_sd,
  est_sd     = si_results$sd,
  converged  = si_results$converged,
  iterations = si_results$iterations
)








# This file is used for full simulation senario.

rm(list= ls())

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

pacman::p_load(
  dplyr,
  purrr,
  tidyr,
  fdrtool,
  devtools,
  openxlsx,
  future,
  furrr
)


pkg_root <- r"(C:\Users\Administrator.DESKTOP-JQG2K7N\桌面\2025 melb\master_project\mitey\mitey)"
results_path <- file.path(pkg_root, "notes", "Simulation_results")
setwd(pkg_root)
devtools::load_all(pkg_root, quiet = TRUE)

#  Data generation function 
rlnorm_natural <- function(n, mean, sd) {
  sdlog <- sqrt(log(1 + (sd / mean)^2))
  meanlog <- log(mean) - 0.5 * sdlog^2
  rlnorm(n, meanlog, sdlog)
}

generate_sim_data <- function(N, true_mean, true_sd, route_weights) {
  n_routes <- round(N * route_weights) 
  
  CP <- rhalfnorm(n_routes[1], theta = sqrt(pi / 2) / (sqrt(2) * true_sd))
  PS <- rlnorm_natural(n_routes[2], mean = true_mean, sd = true_sd)
  PT <- rlnorm_natural(n_routes[3], mean = true_mean, sd = true_sd) +
        rlnorm_natural(n_routes[3], mean = true_mean, sd = true_sd)
  PQ <- rlnorm_natural(n_routes[4], mean = true_mean, sd = true_sd) +
        rlnorm_natural(n_routes[4], mean = true_mean, sd = true_sd) +
        rlnorm_natural(n_routes[4], mean = true_mean, sd = true_sd)
  
  return(round(c(CP, PS, PT, PQ)))
}

# Configure the parameter matrix for the 9 scenarios

scenarios_grid <- expand_grid(
  mean_sd = list(c(2, 2), c(15, 3), c(5, 10)),
  weight_set = list(
    baseline   = c(0.20, 0.50, 0.20, 0.10),
    early_dominant   = c(0.40, 0.35, 0.15, 0.10),
    heavy_tailed = c(0.10, 0.35, 0.35, 0.20)
  )
) %>%

  mutate(
    scenario_id = row_number(),
    true_mean = map_dbl(mean_sd, 1),
    true_sd   = map_dbl(mean_sd, 2),
    weight_name = names(weight_set)
  ) %>%
  select(scenario_id, true_mean, true_sd, weight_name, weight_set)

##############################################
# Main function for simulation execution
run_batch_simulations <- function(grid, N_sample = 200, reps = 100, n_iter = 50, dist_type = "lognormal") {
  
  #Iterate over each scenario in the grid
  all_scenarios <- pmap(grid, function(scenario_id, true_mean, true_sd, weight_name, weight_set) {

  total_scenarios <- nrow(grid)
  # Print a  header when a new Scenario starts
    message(sprintf("\n [Scenario %d/%d] Mean: %g | SD: %g | Weights: %s", 
                    scenario_id, total_scenarios, true_mean, true_sd, weight_name))
    
    # Iterate over the specified number of repetitions
    rep_runs <- future_map(1:reps, function(i) {
      # Use local mitey package for every CPU worker
      devtools::load_all(pkg_root, quiet = TRUE) 
      cat(sprintf("\r  - Running repetition: %d / %d ...", i, reps))
      sim_data <- generate_sim_data(N_sample, true_mean, true_sd, weight_set)
      fit <- si_estim(sim_data, n = n_iter, dist = dist_type, init = c(true_mean, true_sd))
      
      result_row <- tibble(
        scenario_id = scenario_id,
        true_mean   = true_mean,
        true_sd     = true_sd,
        weight_name = weight_name,
        rep         = i,
        est_mean    = fit$mean,
        est_sd      = fit$sd,
        converged   = fit$converged,
        iterations  = fit$iterations,
        loglik      = fit$loglik
      )
      
      history_df <- fit$history %>%
        mutate(
          scenario_id = scenario_id,
          true_mean   = true_mean,
          true_sd     = true_sd,
          weight_name = weight_name,
          rep         = i,
          .before     = 1
        )
      
      # Return both the final result and the iteration history as a list
      list(result = result_row, history = history_df)
    }, .options = furrr_options(seed = TRUE))
    cat("\n  - Completed!\n")
    transposed_reps <- transpose(rep_runs)
    
    list(
      result  = list_rbind(transposed_reps$result),
      history = list_rbind(transposed_reps$history)
    )
  })
  
  #Combine the grouped data frames from all scenarios into the final datasets
  transposed_all <- transpose(all_scenarios)
  
  list(
    results = list_rbind(transposed_all$result),
    history = list_rbind(transposed_all$history)
  )
}

# Using multiple CPU for faster simulation
plan(multisession, workers = max(1, availableCores() - 1))

set.seed(2026)
##############################################
# Execute the simulation and extract outputs
sim_out <- run_batch_simulations(scenarios_grid, N_sample = 200, reps = 100)

final_sim_results <- sim_out$results
iteration_history <- sim_out$history

summary_metrics <- final_sim_results %>%
  group_by(scenario_id, true_mean, true_sd, weight_name) %>%
  summarise(
    n_reps = n(),
    convergence_rate = mean(converged),
    mean_bias = mean(est_mean - true_mean, na.rm = TRUE),
    mean_mse = mean((est_mean - true_mean)^2, na.rm = TRUE),
    sd_bias = mean(est_sd - true_sd, na.rm = TRUE),
    sd_mse = mean((est_sd - true_sd)^2, na.rm = TRUE),
    .groups = "drop"
  )
print(summary_metrics)

#########################################
# Save Results
excel_data <- list(
  "Per_Senario"   = summary_metrics,
  "Per_Repitition" = final_sim_results,
  "Per_Iteration" = iteration_history
)

output_file <- file.path(results_path, "Simulation_All_Results.xlsx")
write.xlsx(excel_data, file = output_file, rowNames = FALSE)


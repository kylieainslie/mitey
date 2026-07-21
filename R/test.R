## Repeated simulation test for mitey::si_estim(dist = "lognormal")
## Output: test_repetition_all_results.csv

library(dplyr)
library(purrr)
library(readr)

# 1. Locate source directory and output path

script_dir <- r"(C:\Users\Administrator.DESKTOP-JQG2K7N\桌面\2025 melb\master_project\mitey\mitey\R)"
output_dir <- r"(C:\Users\Administrator.DESKTOP-JQG2K7N\桌面\2025 melb\master_project\meeting_summarys\simulation_outputs\test_repetition)"


# 2. Simulation parameters
reps <- 50
n_iter <- 50
seed <- 20260720

# Load local R source files
r_files <- setdiff(list.files(script_dir, pattern = "\\.R$", full.names = TRUE), normalizePath(file.path(script_dir, "test.R"), winslash = "/", mustWork = TRUE))
invisible(walk(r_files, source))

# 3. Scenarios setup
scenario_grid <- tibble(
  scenario_id = c("A1_n100", "A2_PS_dominant", "A2_CP_heavy", "A2_higher_generation_heavy", "B1_gamma", "B2_normal"),
  scenario_type = c(rep("correct_specification", 4), rep("misspecification", 2)),
  true_dist = c(rep("lognormal", 4), "gamma", "normal"),
  n_obs = 100,
  true_mean = 6,
  true_sd = 3,
  w_cp = c(0.20, 0.10, 0.40, 0.10, 0.20, 0.20),
  w_ps = c(0.50, 0.70, 0.35, 0.35, 0.50, 0.50),
  w_pt = c(0.20, 0.15, 0.15, 0.35, 0.20, 0.20),
  w_pq = c(0.10, 0.05, 0.10, 0.20, 0.10, 0.10)
) %>% mutate(weights = sprintf("(%.2f, %.2f, %.2f, %.2f)", w_cp, w_ps, w_pt, w_pq))

#Convert Natural-Scale Mean and SD to Lognormal Parameters (meanlog and sdlog)
natural_to_lognormal <- function(mean, sd) {
  sdlog <- sqrt(log(1 + sd^2 / mean^2))
  c(meanlog = log(mean) - 0.5 * sdlog^2, sdlog = sdlog)
}


# Draw Serial Interval Samples from Specified Distribution
draw_si <- function(n, dist, mean, sd) {
  if (n == 0) return(numeric(0))
  switch(dist,
         "lognormal" = { pars <- natural_to_lognormal(mean, sd); rlnorm(n, pars["meanlog"], pars["sdlog"]) },
         "gamma"     = rgamma(n, shape = mean^2 / sd^2, scale = sd^2 / mean),
         "normal"    = rnorm(n, mean = mean, sd = sd),
         stop("Unsupported SI distribution: ", dist)
  )
}

# Simulate Index Case-to-Case Intervals
simulate_icc <- function(n, dist, mean, sd, wts) {
  route <- sample(c("CP", "PS", "PT", "PQ"), size = n, replace = TRUE, prob = wts)
  counts <- table(factor(route, levels = c("CP", "PS", "PT", "PQ")))
  
  val <- map(c("CP1","CP2","PS","PT1","PT2","PQ1","PQ2","PQ3"), ~draw_si(counts[substr(.x,1,2)], dist, mean, sd))
  names(val) <- c("cp1", "cp2", "ps", "pt1", "pt2", "pq1", "pq2", "pq3")
  
  icc <- numeric(n)
  if (dist == "normal") {
    icc[route == "CP"] <- abs(val$cp1 - val$cp2)
    icc[route == "PS"] <- abs(val$ps)
    icc[route == "PT"] <- abs(val$pt1 + val$pt2)
    icc[route == "PQ"] <- abs(val$pq1 + val$pq2 + val$pq3)
  } else {
    icc[route == "CP"] <- abs(val$cp1 - val$cp2)
    icc[route == "PS"] <- val$ps
    icc[route == "PT"] <- val$pt1 + val$pt2
    icc[route == "PQ"] <- val$pq1 + val$pq2 + val$pq3
  }
  pmax(round(icc), 0)
}

# Fit Model for a Single Repetition
# Simulates an ICC dataset for a given scenario
# using `si_estim()`, and records performance metrics and execution time.
fit_one_repetition <- function(scenario_row, repetition_id, seed_value) {
  set.seed(seed_value)
  weights <- c(CP = scenario_row$w_cp, PS = scenario_row$w_ps, PT = scenario_row$w_pt, PQ = scenario_row$w_pq)
  icc <- simulate_icc(scenario_row$n_obs, scenario_row$true_dist, scenario_row$true_mean, scenario_row$true_sd, weights)
  
  start_time <- Sys.time()
  fit <- tryCatch(
    si_estim(dat = icc, dist = "lognormal", n = n_iter, init = c(scenario_row$true_mean, scenario_row$true_sd), tol = 1e-5, n_starts = 1),
    error = function(e) e
  )
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  is_err <- inherits(fit, "error")
  
  tibble(
    scenario_id = scenario_row$scenario_id,
    repetition_id = repetition_id,
    scenario_type = scenario_row$scenario_type,
    true_dist = scenario_row$true_dist,
    fit_dist = "lognormal",
    n_obs = scenario_row$n_obs,
    weights = scenario_row$weights,
    n_iter = n_iter,
    true_mean = scenario_row$true_mean,
    true_sd = scenario_row$true_sd,
    est_mean  = fit$mean,
    est_sd    = fit$sd,
    bias_mean = fit$mean - scenario_row$true_mean,
    bias_sd   = fit$sd - scenario_row$true_sd,
    converged = isTRUE(fit$converged),
    iterations= fit$iterations,
    loglik    = fit$loglik,
    elapsed_seconds = elapsed,
    seed = seed_value,
    error = if (is_err) conditionMessage(fit) else NA_character_
  )
}

# 4. Run simulations across all scenarios and repetitions
all_results <- map_dfr(seq_len(nrow(scenario_grid)), function(scenario_index) {
  scenario_row <- scenario_grid[scenario_index, ]
  
  map_dfr(seq_len(reps), function(rep_id) {
    rep_seed <- seed + scenario_index * 100000L + rep_id
    fit_one_repetition(scenario_row, rep_id, rep_seed)
  })
})

# 5. Output 
output_file <- file.path(output_dir, "test_repetition_all_results.csv")
write_csv(all_results, output_file)
cat(sprintf("Simulation complete. Results saved to: %s\n", output_file))
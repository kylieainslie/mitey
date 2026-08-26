rm(list = ls())

library(dplyr)
library(purrr)
library(fdrtool)
library(devtools)
library(openxlsx)

pkg_root <- r"(C:\Users\Administrator.DESKTOP-JQG2K7N\桌面\2025 melb\master_project\mitey\mitey)"
setwd(pkg_root)
load_all()

generate_normal_icc <- function(N, true_mean, true_sd, route_weights) {
  n_routes <- round(N * route_weights)

  CP <- rhalfnorm(n_routes[1], theta = sqrt(pi / 2) / (sqrt(2) * true_sd))
  PS <- rnorm(n_routes[2], mean = true_mean, sd = true_sd)
  PT <- rnorm(n_routes[3], mean = 2 * true_mean, sd = sqrt(2) * true_sd)
  PQ <- rnorm(n_routes[4], mean = 3 * true_mean, sd = sqrt(3) * true_sd)

  round(c(CP, PS, PT, PQ))
}
# For normal, there are 8 free parameters
get_k <- function(dist) {
  if (dist == "normal") return(8)
  if (dist %in% c("gamma", "lognormal")) return(5)
}

fit_three_models <- function(dat, true_mean, true_sd) {
  dists <- c("normal", "gamma", "lognormal")

  map_dfr(dists, function(d) {
    fit <- si_estim(dat, dist = d, init = c(true_mean, true_sd))
    k <- get_k(d)

    tibble(
      fit_dist = d,
      est_mean = fit$mean,
      est_sd = fit$sd,
      loglik = fit$loglik,
      AIC = -2 * fit$loglik + 2 * k,
      BIC = -2 * fit$loglik + log(length(dat)) * k,
      converged = fit$converged
    )
  })
}

set.seed(2026)

N <- 200
n_reps <- 10
true_mean <- 5
true_sd <- 5
route_weights <- c(0.20, 0.50, 0.20, 0.10)

model_selection_results <- map_dfr(1:n_reps, function(i) {
  dat <- generate_normal_icc(N, true_mean, true_sd, route_weights)
  dat <- dat[dat >= 0] # Gamma and lognormal can't handle negative ICC values
  fit_three_models(dat, true_mean, true_sd) %>%
    mutate(
      rep = i,
      true_dist = "normal",
      true_mean = true_mean,
      true_sd = true_sd,
      .before = 1
    )
})

# Calculate the rate of selecting normal model under AIC/BIC criteria
selection_rates <- model_selection_results %>%
  group_by(rep) %>%
  summarise(
    best_AIC = fit_dist[which.min(AIC)],
    best_BIC = fit_dist[which.min(BIC)],
    .groups = "drop"
  ) %>%
  summarise(
    AIC_select_normal_rate = mean(best_AIC == "normal", na.rm = TRUE),
    BIC_select_normal_rate = mean(best_BIC == "normal", na.rm = TRUE)
  )

# Calculate mean AIC/BIC for each model
mean_metrics <- model_selection_results %>%
  group_by(fit_dist) %>%
  summarise(
    mean_AIC = mean(AIC, na.rm = TRUE),
    mean_BIC = mean(BIC, na.rm = TRUE),
    mean_loglik = mean(loglik, na.rm = TRUE),
    .groups = "drop"
  )

print(selection_rates)
print(mean_metrics)


# Save Results
excel_data <- list(
  "Raw_Results" = model_selection_results,
  "Selection_Rates" = selection_rates,
  "Mean_Metrics" = mean_metrics
)

write.xlsx(
  excel_data,
  file = file.path(pkg_root, "notes", "Simulation_results", "model_selection_normal.xlsx"),
  rowNames = FALSE
)
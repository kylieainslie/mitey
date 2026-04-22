compare_n_routes <- function(
    dat,
    n_routes_max = 8L,
    dist = "normal",
    n = 50,
    ...
) {
  results <- vector("list", n_routes_max - 1L)

  for (n_routes in 2:n_routes_max) {
    fit <- si_estim(dat, n = n, dist = dist, n_routes = n_routes, ...)
    print(plot_si_fit_result(fit,dat))

    if (dist == "normal") {
      n_params <- 2L + (2L * n_routes - 2L)
    } else {
      n_params <- 2L + (n_routes - 1L)
    }

    n_obs <- length(dat)
    aic <- -2 * fit$loglik + 2 * n_params
    bic <- -2 * fit$loglik + log(n_obs) * n_params

    results[[n_routes - 1L]] <- data.frame(
      n_routes   = n_routes,
      mu         = fit$mean,
      sigma      = fit$sd,
      loglik     = fit$loglik,
      aic        = aic,
      bic        = bic,
      n_params   = n_params,
      converged  = fit$converged,
      iterations = fit$iterations
    )
  }

  results_df <- do.call(rbind, results)

  # identifier le meilleur selon BIC
  best_bic <- results_df$n_routes[which.min(results_df$bic)]
  best_aic <- results_df$n_routes[which.min(results_df$aic)]

  # calculer les n_optimal
  n_opt1 <- n_optimal1(dat)
  n_opt2 <- n_optimal2(dat)
  n_opt3 <- n_optimal3(dat)

  cat("=== Comparaison des n_routes ===\n")
  print(results_df, digits = 4)
  cat("\nMeilleur n_routes selon BIC :", best_bic, "\n")
  cat("Meilleur n_routes selon AIC :", best_aic, "\n")

  invisible(results_df)
}

#' @export
compare_n_routes <- function(
    dat,
    n_routes_max = 8L,
    dist = "normal",
    n = 50,
    ...
) {
  results <- vector("list", n_routes_max - 1L)
  weights_list <- vector("list", n_routes_max - 1L)

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

    wts        <- fit$wts
    wt_names   <- paste0("w", seq_along(wts))   # w1, w2, w3, ...
    wts_df     <- as.data.frame(t(wts))
    names(wts_df) <- wt_names
    wts_df$n_routes <- n_routes
    weights_list[[n_routes - 1L]] <- wts_df
  }

  results_df <- do.call(rbind, results)

  weights_df <- do.call(
    dplyr::bind_rows,          # bind_rows gère les colonnes manquantes → NA
    weights_list
  )
  # Réordonner : n_routes en première colonne
  weights_df <- weights_df[, c("n_routes",
                               setdiff(names(weights_df), "n_routes"))]


  cat("\n")
  cat("\n")
  cat("\n=== Components weights by model ===\n")
  cat("(CP = w1 ; PS = w2/w3 ; PT = w4/w5 ; etc. — NA = no component)\n\n")
  print(weights_df, digits = 3, row.names = FALSE)

  # Identify the best using BIC
  best_bic <- results_df$n_routes[which.min(results_df$bic)]
  best_aic <- results_df$n_routes[which.min(results_df$aic)]

  # Compute the n_optimal
  n_opt1 <- n_optimal1(dat)
  n_opt2 <- n_optimal2(dat)

  cat("=== Comparison of the n_routes ===\n")
  print(results_df, digits = 4)
  cat("\nBest n_routes according to BIC :", best_bic, "\n")
  cat("Best n_routes according to AIC :", best_aic, "\n")

  invisible(results_df)
}

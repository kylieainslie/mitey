#' Compare Serial Interval Models Across Different Numbers of Transmission Routes
#'
#' Fits the serial interval mixture model for each number of transmission routes
#' from 2 up to \code{n_routes_max}, computes model selection criteria (AIC and BIC)
#' for each fit, and prints a summary table to help identify the optimal number
#' of routes.
#'
#' @param dat numeric vector; index case-to-case (ICC) intervals in days.
#' @param n_routes_max integer; maximum number of transmission routes to evaluate.
#'   Models are fitted for \code{n_routes} from 2 to \code{n_routes_max}.
#'   Defaults to \code{8L}.
#' @param dist character; assumed parametric family for the serial interval
#'   distribution. Must be either \code{"normal"} (default) or \code{"gamma"}.
#' @param n integer; number of EM algorithm iterations passed to
#'   \code{\link{si_estim}}. Defaults to 50.
#' @param ... additional arguments passed to \code{\link{si_estim}}.
#'
#' @details
#' For each value of \code{n_routes} in \code{2:n_routes_max}, the function:
#' \enumerate{
#'   \item Calls \code{\link{si_estim}} to fit the mixture model and retrieve
#'     the log-likelihood and estimated weights.
#'   \item Counts the number of free parameters \eqn{p}:
#'     \itemize{
#'       \item Normal distribution: \eqn{p = 2 + (2 \times n\_routes - 2)}
#'       \item Gamma distribution: \eqn{p = 2 + (n\_routes - 1)}
#'     }
#'   \item Computes AIC and BIC:
#'     \deqn{AIC = -2 \ell + 2p}
#'     \deqn{BIC = -2 \ell + p \log(n)}
#'   \item Prints the fitted mixture plot via \code{\link{plot_si_fit_result}}.
#' }
#' The model with the lowest BIC (or AIC) is reported as the best-fitting
#' number of routes.
#'
#' @return A data frame (returned invisibly) with one row per value of
#'   \code{n_routes} and the following columns:
#' \describe{
#'   \item{n_routes}{Number of transmission routes fitted.}
#'   \item{mu}{Estimated mean of the serial interval (days).}
#'   \item{sigma}{Estimated standard deviation of the serial interval (days).}
#'   \item{loglik}{Log-likelihood of the fitted model.}
#'   \item{aic}{Akaike Information Criterion.}
#'   \item{bic}{Bayesian Information Criterion.}
#'   \item{n_params}{Number of free parameters in the model.}
#'   \item{converged}{Logical; whether the EM algorithm converged.}
#'   \item{iterations}{Number of EM iterations performed.}
#' }
#' A summary table and the best \code{n_routes} according to BIC and AIC are
#' printed as side effects.
#'
#' @seealso \code{\link{si_estim}}, \code{\link{plot_si_fit_result}}
#'
#' @export
#' @examples
#' \donttest{
#' set.seed(42)
#' icc_data <- c(
#'   abs(rnorm(30, mean = 0,  sd = 2)),
#'   rnorm(80, mean = 12, sd = 3),
#'   rnorm(30, mean = 24, sd = 4)
#' )
#' icc_data <- round(pmax(icc_data, 0))
#'
#' # Compare 2 to 6 routes with normal distribution
#' results <- compare_n_routes(icc_data, n_routes_max = 6, dist = "normal", n = 50)
#'
#' # Compare using gamma distribution
#' results_gam <- compare_n_routes(icc_data, n_routes_max = 5, dist = "gamma", n = 50)
#' }

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
    print(plot_si_fit_result(fit,dat,dist=dist))

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
    dplyr::bind_rows,
    weights_list
  )

  weights_df <- weights_df[, c("n_routes",
                               setdiff(names(weights_df), "n_routes"))]


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

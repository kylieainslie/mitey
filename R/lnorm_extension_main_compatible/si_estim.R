#' Estimate Serial Interval Distribution Using the Vink Method
#'
#' Estimates the mean and standard deviation of the serial interval distribution from
#' outbreak data using the Expectation-Maximization (EM) algorithm developed by Vink
#' et al. (2014). The serial interval is defined as the time between symptom onset
#' in a primary case and symptom onset in a secondary case infected by that primary case.
#'
#' @param dat numeric vector; index case-to-case (ICC) intervals in days.
#' @param n integer; number of EM algorithm iterations to perform. Defaults to 50.
#' @param dist character; the assumed parametric family for the serial interval
#'   distribution. Must be "normal" (default), "gamma", or "lognormal".
#' @param init numeric vector of length 2; initial values for the mean and standard
#'   deviation. If NULL (default), uses the sample mean and sample standard deviation.
#' @param tol numeric; convergence tolerance for changes in the SI mean, SD, and component weights. Defaults to 1e-6. Set to 0 to disable early stopping.
#' @param n_starts integer; number of random restarts for the EM algorithm. Defaults to 1.
#' @param n_routes integer; number of transmission routes to model. Must be >= 2.
#'   Defaults to 4 (Co-Primary, Primary-Secondary, Primary-Tertiary, Primary-Quaternary).
#'   Increasing this allows modelling longer transmission chains.
#' @param wind The window censure interval
#' @param return_history logical; whether to return iteration history for the
#'   selected restart. Defaults to FALSE.
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{mean}: Estimated mean of the serial interval distribution (days)
#'   \item \code{sd}: Estimated standard deviation of the serial interval distribution (days)
#'   \item \code{wts}: Numeric vector of estimated component weights. Length is
#'     2*n_routes - 1 for normal distribution, n_routes for gamma and lognormal distributions.
#'   \item \code{converged}: Logical indicating whether the algorithm converged.
#'   \item \code{iterations}: Integer indicating the number of iterations performed.
#'   \item \code{loglik}: Log-likelihood of the fitted model.
#'   \item \code{history}: Optional data frame of iteration-wise SI estimates,
#'     parameter and weight changes, and component weights for the selected restart.
#'     Included only when \code{return_history = TRUE}.
#'   \item \code{n_restarts}: Number of restarts performed.
#'   \item \code{n_routes}: Number of transmission routes used.
#'   \item \code{wind}: Window censoring interval used.
#' }
#'
#' @references
#' Vink MA, Bootsma MCJ, Wallinga J (2014). Serial intervals of respiratory infectious
#' diseases: A systematic review and analysis. American Journal of Epidemiology,
#' 180(9), 865-875. \doi{10.1093/aje/kwu209}
#'
#' @export
#' @importFrom stats weighted.mean optim
#'
#' @examples
#' # Example 1: Basic usage with simulated data, default 4 routes
#' set.seed(123)
#' simulated_icc <- c(
#'   rep(1, 20), rep(2, 25), rep(3, 15), rep(4, 8)
#' )
#' result <- si_estim(simulated_icc)
#'
#' \donttest{
#' # Example 2: Using 5 routes
#' result_5routes <- si_estim(simulated_icc, n_routes = 5)
#'
#' # Example 3: Using gamma distribution with 3 routes
#' result_gamma <- si_estim(simulated_icc, dist = "gamma", n_routes = 3)
#'
#' # Example 4: Using lognormal distribution
#' result_lognormal <- si_estim(simulated_icc, dist = "lognormal")
#' }
#'
si_estim <- function(
    dat,
    n        = 50,
    dist     = "normal",
    init     = NULL,
    tol      = 1e-6,
    n_starts = 1,
    n_routes = 4L,
    wind = 1,
    return_history = FALSE
) {
  ## --- Input validation ---

  if (any(is.na(dat))) {
    stop("Data contains NA values. Please remove or impute NA values before analysis.")
  }
  if (!is.numeric(dat)) {
    stop("Input data must be numeric.")
  }
  if (any(dat < 0)) {
    warning(
      "Data contains negative values. While the Vink method can handle negative serial intervals ",
      "(e.g., co-primary infections), please ensure this is intended for your analysis."
    )
  }
  if (length(dat) < 2) {
    stop("Need at least 2 data points to estimate serial interval parameters.")
  }
  if (!dist %in% c("normal", "gamma", "lognormal")) {
    stop("Distribution must be either 'normal', 'gamma', or 'lognormal'.")
  }
  if (dist %in% c("gamma", "lognormal") && any(dat < 0)) {
    stop("Gamma and lognormal distributions cannot be used with negative values. Please use 'normal' distribution instead.")
  }
  if (is.null(init)) {
    init <- c(mean(dat), sd(dat))
  }
  if (length(init) != 2) {
    stop("Initial values must be a vector of length 2 (mean and sd).")
  }
  if (any(is.na(init))) {
    stop("Initial values cannot contain NA.")
  }
  if (init[2] <= 0) {
    stop("Initial standard deviation must be positive.")
  }
  if (dist == "lognormal" && (!is.finite(init[1]) || init[1] <= 0)) {
    stop("Initial mean must be positive and finite for lognormal distribution.")
  }
  if (!is.numeric(tol) || length(tol) != 1 || is.na(tol) ||
      !is.finite(tol) || tol < 0) {
    stop("Tolerance must be a single non-negative finite numeric value.")
  }
  if (!is.numeric(n_starts) || length(n_starts) != 1 || is.na(n_starts) ||
      !is.finite(n_starts) || n_starts < 1 || n_starts != floor(n_starts)) {
    stop("n_starts must be a positive integer.")
  }
  n_starts <- as.integer(n_starts)

  if (!is.numeric(n_routes) || length(n_routes) != 1 || is.na(n_routes) ||
      !is.finite(n_routes) || n_routes < 2 || n_routes != floor(n_routes)) {
    stop("n_routes must be an integer >= 2.")
  }
  n_routes <- as.integer(n_routes)

  if (!is.logical(return_history) || length(return_history) != 1L ||
      is.na(return_history)) {
    stop("return_history must be TRUE or FALSE.")
  }

  ## --- Setup ---

  j   <- length(dat)
  if (dist == "lognormal") {
    # Use a zero-safe proxy only for the lognormal M-step.
    lognormal_mstep_dat <- ifelse(dat == 0, 1, dat)
  }
  dat <- ifelse(dat == 0, 0.00001, dat)

  # Build component vector dynamically
  if (dist == "normal") {
    comp_vec <- c(1L, unlist(lapply(seq_len(n_routes - 1), function(i) c(2L*i, 2L*i + 1L))))
  } else {
    comp_vec <- c(1L, 2L * seq_len(n_routes - 1L))
  }

  n_comp <- length(comp_vec)

  # Data-derived bounds for random starting points
  data_sd    <- sd(dat)
  data_range <- range(dat)

  # Generate starting points
  starting_points <- vector("list", n_starts)
  starting_points[[1]] <- init
  if (n_starts > 1) {
    for (s in 2:n_starts) {
      mu_start    <- stats::runif(1, min = data_range[1], max = data_range[2])
      sigma_start <- stats::runif(1, min = data_sd * 0.5, max = data_sd * 2)
      starting_points[[s]] <- c(mu_start, sigma_start)
    }
  }

  ## --- EM helper ---

  run_single_em <- function(mu_init, sigma_init) {
    # E-step
    mu    <- mu_init
    sigma <- sigma_init
    converged      <- FALSE
    iterations_used <- n
    w <- rep(1 / n_comp, n_comp)
    if (return_history) history <- vector("list", n)

    for (k_iter in 1:n) {
      mu_prev    <- mu
      sigma_prev <- sigma
      w_prev     <- w

      tau <- matrix(0, nrow = n_comp, ncol = j)

      for (l in 1:j) {
        use_lower <- (dat[l] != 0.00001)
        for (comp_idx in seq_len(n_comp)) {
          tau[comp_idx, l] <- integrate_component(
            dat[l], mu, sigma,
            comp  = comp_vec[comp_idx],
            dist  = dist,
            lower = use_lower,
            wind = wind
          )
        }
      }

      # Weight component probabilities and normalize posterior responsibilities
      tau <- sweep(tau, 1, w, "*")
      denom <- colSums(tau)
      tau   <- sweep(tau, 2, denom, "/")


      # Component weights
      w <- rowSums(tau) / j

      # M-step: update mu and sigma
      if (dist == "normal") {
        # Use the first positive-route component (index 2 in comp_vec = component 2)

        pos_idx <- which(comp_vec == 2L)
        mu    <- weighted.mean(dat, tau[pos_idx, ])
        sigma <- sqrt(weighted_var(dat, tau[pos_idx, ]))
      } else if (dist == "gamma") {
        pos_idx <- which(comp_vec == 2L)
        opt <- optim(
          par     = c(mu, sigma),
          fn      = wt_loglik,
          tau2    = tau[pos_idx, ],
          dat     = dat,
          gr      = NULL,
          method  = "BFGS",
          hessian = FALSE
        )
        mu    <- opt$par[1]
        sigma <- opt$par[2]
      } else if (dist == "lognormal") {
        # Update PS parameters on the log scale, then return natural-scale estimates.
        pos_idx <- which(comp_vec == 2L)
        log_dat <- log(lognormal_mstep_dat)
        mu_log <- weighted.mean(log_dat, tau[pos_idx, ])
        sigma_log <- sqrt(weighted_var(log_dat, tau[pos_idx, ]))
        mu <- exp(mu_log + 0.5 * sigma_log^2)
        sigma <- sqrt((exp(sigma_log^2) - 1) * exp(2 * mu_log + sigma_log^2))
      }

      # Convergence check: SI parameters and component weights
      mu_change <- NA_real_
      sigma_change <- NA_real_
      w_change <- NA_real_
      if (k_iter > 1) {
        mu_change    <- abs(mu - mu_prev)    / (abs(mu_prev)    + .Machine$double.eps)
        sigma_change <- abs(sigma - sigma_prev) / (abs(sigma_prev) + .Machine$double.eps)
        w_change <- max(abs(w - w_prev))
      }

      if (return_history) {
        history[[k_iter]] <- data.frame(
          iteration = k_iter,
          mean = mu,
          sd = sigma,
          mu_change = mu_change,
          sigma_change = sigma_change,
          w_change = w_change,
          as.list(setNames(w, paste0("w_", comp_vec)))
        )
      }

      if (tol > 0 && k_iter > 1 &&
          mu_change < tol && sigma_change < tol && w_change < tol) {
        converged       <- TRUE
        iterations_used <- k_iter
        break
      }
    }

    loglik <- calculate_mixture_loglik(dat, mu, sigma, w, comp_vec, dist, wind)

    result <- list(
      mean       = mu,
      sd         = sigma,
      wts        = w,
      converged  = converged,
      iterations = iterations_used,
      loglik     = loglik
    )
    if (return_history) {
      result$history <- do.call(rbind, history[seq_len(iterations_used)])
    }
    return(result)

  }

  ## --- Run restarts ---

  best_result <- NULL
  best_loglik <- -Inf

  for (s in seq_len(n_starts)) {
    result <- run_single_em(starting_points[[s]][1], starting_points[[s]][2])
    if (is.finite(result$loglik) && result$loglik > best_loglik) {
      best_loglik <- result$loglik
      best_result <- result
    }
  }

  if (is.null(best_result)) {
    best_result <- result
  }

  best_result$n_restarts <- n_starts
  best_result$n_routes   <- n_routes
  best_result$wind   <- wind

  return(best_result)
}


#' Calculate Log-Likelihood for Mixture Model
#'
#' Internal function to calculate the log-likelihood of the fitted mixture model.
#'
#' @param dat numeric vector; the data
#' @param mu numeric; estimated mean
#' @param sigma numeric; estimated standard deviation
#' @param wts numeric vector; component weights
#' @param comp_vec integer vector; component indices
#' @param dist character; distribution type ("normal", "gamma", or "lognormal")
#' @param wind The window censure interval
#'
#' @return numeric; log-likelihood value
#' @keywords internal
calculate_mixture_loglik <- function(dat, mu, sigma, wts, comp_vec, dist,wind) {
  j      <- length(dat)
  loglik <- 0

  for (l in 1:j) {
    prob      <- 0
    use_lower <- (dat[l] != 0.00001)
    for (comp_idx in seq_along(comp_vec)) {
      comp      <- comp_vec[comp_idx]
      comp_prob <- integrate_component(
        dat[l], mu, sigma,
        comp  = comp,
        dist  = dist,
        lower = use_lower,
        wind = wind
      )
      prob <- prob + wts[comp_idx] * comp_prob
    }
    if (prob > 0) {
      loglik <- loglik + log(prob)
    } else {
      loglik <- loglik + log(.Machine$double.xmin)
    }
  }

  return(loglik)
}

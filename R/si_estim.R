#' Estimate Serial Interval Distribution Using the Vink Method
#'
#' Estimates the mean and standard deviation of the serial interval distribution from
#' outbreak data using the Expectation-Maximization (EM) algorithm developed by Vink
#' et al. (2014). The serial interval is defined as the time between symptom onset
#' in a primary case and symptom onset in a secondary case infected by that primary case.
#'
#' The Vink method addresses the challenge that individual transmission pairs are
#' typically unknown in outbreak investigations. Instead, it uses index case-to-case
#' (ICC) intervals - the time differences between the case with earliest symptom onset
#' (index case) and all other cases - to infer the underlying serial interval
#' distribution through a mixture modeling approach.
#'
#' \strong{Methodological Approach:}
#'
#' The method models ICC intervals as arising from a mixture of four transmission routes:
#' \itemize{
#'   \item \strong{Co-Primary (CP)}: Cases infected simultaneously from the same source
#'   \item \strong{Primary-Secondary (PS)}: Direct transmission from index case
#'   \item \strong{Primary-Tertiary (PT)}: Second-generation transmission
#'   \item \strong{Primary-Quaternary (PQ)}: Third-generation transmission
#' }
#'
#' The EM algorithm iteratively:
#' \enumerate{
#'   \item \strong{E-step}: Calculates the probability that each ICC interval belongs
#'         to each transmission route component
#'   \item \strong{M-step}: Updates the serial interval parameters (mean, standard deviation)
#'         and component weights based on these probabilities
#' }
#'
#' \strong{Distribution Choice:}
#'
#' \itemize{
#'   \item \strong{Normal distribution}: Allows negative serial intervals (useful for
#'         modeling co-primary infections) and uses 7 components (positive and negative
#'         pairs for PS, PT, PQ routes plus CP)
#'   \item \strong{Gamma distribution}: Restricts to positive values only, uses 4 components
#'         (CP, PS, PT, PQ without negative pairs). Recommended when negative serial
#'         intervals are epidemiologically implausible
#' }
#'
#' \strong{Key Assumptions:}
#'
#' \itemize{
#'   \item The case with earliest symptom onset is the index case
#'   \item Transmission occurs through at most 4 generations
#'   \item Serial intervals follow the specified parametric distribution
#'   \item Cases represent a single, homogeneously-mixing outbreak
#' }
#'
#' @param dat numeric vector; index case-to-case (ICC) intervals in days. These are
#'            calculated as the time difference between symptom onset in each case and
#'            symptom onset in the index case (case with earliest onset). Must contain
#'            at least 2 values. Values should be non-negative in most epidemiological
#'            contexts, though negative values are allowed for normal distribution
#' @param n integer; number of EM algorithm iterations to perform. More iterations
#'          generally improve convergence but increase computation time. Defaults to 50,
#'          which is typically sufficient for convergence
#' @param dist character; the assumed parametric family for the serial interval distribution. Must be either:
#' \itemize{
#'   \item \code{"normal"} (default): Allows negative serial intervals, uses 7 mixture components
#'   \item \code{"gamma"}: Restricts to positive serial intervals, uses 4 mixture components
#' }
#' @param init numeric vector of length 2; initial values for the mean and standard
#'             deviation to start the EM algorithm. If \code{NULL} (default), uses the
#'             sample mean and sample standard deviation of the input data. Providing
#'             good initial values can improve convergence, especially for challenging datasets
#' @param tol numeric; convergence tolerance for the EM algorithm. The algorithm stops
#'            early if the relative change in both mean and standard deviation between
#'            iterations is less than this value. Set to 0 to disable early stopping
#'            and always run all \code{n} iterations. Defaults to 1e-6
#' @param n_starts integer; number of random restarts for the EM algorithm. The algorithm
#'                 can converge to local optima, especially when initial values are far from
#'                 the true parameters. Using multiple restarts with different starting points
#'                 helps find the global optimum. The first restart uses the provided \code{init}
#'                 values (or data-derived values if \code{init = NULL}), and subsequent restarts
#'                 use random starting points sampled from the data range. The result with the
#'                 highest log-likelihood is returned. Defaults to 1 (no additional restarts)
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{mean}: Estimated mean of the serial interval distribution (days)
#'   \item \code{sd}: Estimated standard deviation of the serial interval distribution (days)
#'   \item \code{wts}: Numeric vector of estimated component weights representing the
#'         probability that cases belong to each transmission route. Length depends
#'         on distribution choice (7 for normal, 4 for gamma)
#'   \item \code{converged}: Logical indicating whether the algorithm converged before
#'         reaching the maximum number of iterations
#'   \item \code{iterations}: Integer indicating the number of iterations performed
#'   \item \code{loglik}: Log-likelihood of the fitted model (used for model comparison)
#'   \item \code{n_restarts}: Number of restarts performed
#' }
#'
#' @details
#' \strong{Input Data Preparation:}
#'
#' To prepare ICC intervals from outbreak data:
#' \enumerate{
#'   \item Identify the case with the earliest symptom onset date (index case)
#'   \item Calculate the time difference (in days) between each case's onset
#'         date and the index case onset date
#'   \item The resulting values are the ICC intervals for input to this function
#' }
#'
#' \strong{Convergence and Diagnostics:}
#'
#' The EM algorithm typically converges within 20-50 iterations. Users should:
#' \itemize{
#'   \item Examine the fitted distribution using \code{\link{plot_si_fit}}
#'   \item Consider alternative distribution choices if fit is poor
#'   \item Try different initial values if results seem unreasonable
#'   \item Ensure adequate sample size (generally >20 cases recommended)
#' }
#'
#'
#' @seealso \code{\link{plot_si_fit}} for diagnostic visualization,
#'          \code{\link{integrate_component}} for the underlying likelihood calculations
#'
#' @references
#' Vink MA, Bootsma MCJ, Wallinga J (2014). Serial intervals of respiratory infectious
#' diseases: A systematic review and analysis. American Journal of Epidemiology,
#' 180(9), 865-875. \doi{10.1093/aje/kwu209}
#'
#' @export
#' @importFrom stats weighted.mean
#' @importFrom stats optim
#'
#' @examples
#' # Example 1:Basic usage with simulated data
#' set.seed(123)
#' simulated_icc <- c(
#'   rep(1, 20),   # Short intervals (co-primary cases)
#'   rep(2, 25),   # Medium intervals (primary-secondary)
#'   rep(3, 15),   # Longer intervals (higher generation)
#'   rep(4, 8)
#' )
#'
#' result <- si_estim(simulated_icc)
#'
#' \donttest{
#' # Example 2: Larger simulated outbreak, specifying distribution
#' large_icc <- c(
#'   rep(1, 38),   # Short intervals (co-primary cases)
#'   rep(2, 39),   #
#'   rep(3, 30),   # Medium intervals (primary-secondary)
#'   rep(4, 17),   #
#'   rep(5, 7),    # Longer intervals (higher generation)
#'   rep(6, 4),
#'   rep(7, 2)
#' )
#'
#' result_normal <- si_estim(large_icc, dist = "normal")
#' result_gamma <- si_estim(large_icc, dist = "gamma")
#'
#' # Example 3: Using custom initial values
#' result_custom <- si_estim(large_icc, dist = "normal", init = c(3.0, 1.5))
#'
#' # Example 4: Specify iterations
#' result_iter <- si_estim(large_icc, n=100)
#'
#' # Example 5: Check convergence status
#' result_conv <- si_estim(large_icc)
#' if (result_conv$converged) {
#'   message("Converged in ", result_conv$iterations, " iterations")
#' }
#'
#' # Example 6: Use multiple restarts to avoid local optima
#' # Useful when initial values might be far from true parameters
#' result_restarts <- si_estim(large_icc, n_starts = 5)
#' message("Best result from ", result_restarts$n_restarts, " restarts")
#'
#' }
#'
si_estim <- function(
  dat,
  n = 50,
  dist = "normal",
  init = NULL,
  tol = 1e-6,
  n_starts = 1
) {
  ## Check inputs
  # Check inputs for NA values
  if (any(is.na(dat))) {
    stop(
      "Data contains NA values. Please remove or impute NA values before analysis."
    )
  }

  # Check for non-numeric input
  if (!is.numeric(dat)) {
    stop("Input data must be numeric.")
  }

  # Check for negative values and issue warning
  if (any(dat < 0)) {
    warning(
      "Data contains negative values. While the Vink method can handle negative serial intervals ",
      "(e.g., co-primary infections), please ensure this is intended for your analysis."
    )
  }

  # Check for single value
  if (length(dat) < 2) {
    stop("Need at least 2 data points to estimate serial interval parameters.")
  }

  # Check distribution parameter
  if (!dist %in% c("normal", "gamma")) {
    stop("Distribution must be either 'normal' or 'gamma'.")
  }

  # Check for gamma distribution with negative values
  if (dist == "gamma" && any(dat < 0)) {
    stop(
      "Gamma distribution cannot be used with negative values. Please use 'normal' distribution instead."
    )
  }

  # Set initial values if not provided
  if (is.null(init)) {
    init <- c(mean(dat), sd(dat))
  }

  # Check initial values
  if (length(init) != 2) {
    stop("Initial values must be a vector of length 2 (mean and sd).")
  }

  if (any(is.na(init))) {
    stop("Initial values cannot contain NA.")
  }

  if (init[2] <= 0) {
    stop("Initial standard deviation must be positive.")
  }

  # Check tolerance parameter
  if (!is.numeric(tol) || length(tol) != 1 || is.na(tol) || !is.finite(tol) || tol < 0) {
    stop("Tolerance must be a single non-negative finite numeric value.")
  }

  # Check n_starts parameter
  if (!is.numeric(n_starts) || length(n_starts) != 1 || is.na(n_starts) ||
      !is.finite(n_starts) || n_starts < 1 || n_starts != floor(n_starts)) {
    stop("n_starts must be a positive integer.")
  }
  n_starts <- as.integer(n_starts)

  ## Vink et al. implementation code
  j <- length(dat)
  dat <- ifelse(dat == 0, 0.00001, dat)

  # Components depend on specified distribution
  if (dist == "normal") {
    comp_vec <- 1:7
  } else if (dist == "gamma") {
    comp_vec <- c(1, 2, 4, 6)
  }

  # Data-derived bounds for random starting points
  data_sd <- sd(dat)
  data_range <- range(dat)

  # Generate starting points for multiple restarts
  starting_points <- vector("list", n_starts)

  # First starting point: user-provided or data-derived
  starting_points[[1]] <- init  # Already set to c(mean(dat), sd(dat)) if NULL

  # Additional random starting points
  if (n_starts > 1) {
    for (s in 2:n_starts) {
      # Random mean within data range
      mu_start <- stats::runif(1, min = data_range[1], max = data_range[2])
      # Random SD between 0.5x and 2x the data SD
      sigma_start <- stats::runif(1, min = data_sd * 0.5, max = data_sd * 2)
      starting_points[[s]] <- c(mu_start, sigma_start)
    }
  }

  # Helper function to run a single EM
  run_single_em <- function(mu_init, sigma_init) {
    mu <- mu_init
    sigma <- sigma_init
    converged <- FALSE
    iterations_used <- n

    for (k in 1:n) {
      mu_prev <- mu
      sigma_prev <- sigma

      tau <- matrix(0, nrow = length(comp_vec), ncol = j)

      for (l in 1:j) {
        if (dat[l] == 0.00001) {
          for (comp in seq_along(comp_vec)) {
            tau[comp, l] <- integrate_component(
              dat[l], mu, sigma,
              comp = comp_vec[comp], dist = dist, lower = FALSE
            )
          }
        } else {
          for (comp in seq_along(comp_vec)) {
            tau[comp, l] <- integrate_component(
              dat[l], mu, sigma,
              comp = comp_vec[comp], dist = dist, lower = TRUE
            )
          }
        }
      }

      # Normalize tau
      denom <- colSums(tau)
      tau <- sweep(tau, 2, denom, "/")

      # Calculate the weights
      w <- rowSums(tau) / j

      # Update parameters
      if (dist == "normal") {
        mu <- weighted.mean(dat, tau[2, ])
        sigma <- sqrt(weighted_var(dat, tau[2, ]))
      } else if (dist == "gamma") {
        opt <- optim(
          par = c(mu, sigma),
          wt_loglik,
          tau2 = tau[2, ],
          dat = dat,
          gr = NULL,
          method = "BFGS",
          hessian = FALSE
        )
        mu <- opt$par[1]
        sigma <- opt$par[2]
      }

      # Check for convergence
      if (tol > 0 && k > 1) {
        mu_change <- abs(mu - mu_prev) / (abs(mu_prev) + .Machine$double.eps)
        sigma_change <- abs(sigma - sigma_prev) / (abs(sigma_prev) + .Machine$double.eps)

        if (mu_change < tol && sigma_change < tol) {
          converged <- TRUE
          iterations_used <- k
          break
        }
      }
    }

    # Calculate log-likelihood for model comparison
    loglik <- calculate_mixture_loglik(dat, mu, sigma, w, comp_vec, dist)

    list(
      mean = mu,
      sd = sigma,
      wts = w,
      converged = converged,
      iterations = iterations_used,
      loglik = loglik
    )
  }

  # Run EM for each starting point and keep track of results
  best_result <- NULL
  best_loglik <- -Inf

  for (s in seq_len(n_starts)) {
    result <- run_single_em(starting_points[[s]][1], starting_points[[s]][2])

    if (is.finite(result$loglik) && result$loglik > best_loglik) {
      best_loglik <- result$loglik
      best_result <- result
    }
  }

  # If no valid result found (shouldn't happen), use the last result

  if (is.null(best_result)) {
    best_result <- result
  }

  # Add n_restarts to output
  best_result$n_restarts <- n_starts

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
#' @param dist character; distribution type ("normal" or "gamma")
#'
#' @return numeric; log-likelihood value
#' @keywords internal
calculate_mixture_loglik <- function(dat, mu, sigma, wts, comp_vec, dist) {
  j <- length(dat)
  loglik <- 0

  for (l in 1:j) {
    # Calculate probability for each component
    prob <- 0
    for (comp_idx in seq_along(comp_vec)) {
      comp <- comp_vec[comp_idx]
      if (dat[l] == 0.00001) {
        comp_prob <- integrate_component(
          dat[l], mu, sigma,
          comp = comp, dist = dist, lower = FALSE
        )
      } else {
        comp_prob <- integrate_component(
          dat[l], mu, sigma,
          comp = comp, dist = dist, lower = TRUE
        )
      }
      prob <- prob + wts[comp_idx] * comp_prob
    }

    if (prob > 0) {
      loglik <- loglik + log(prob)
    } else {
      loglik <- loglik + log(.Machine$double.xmin)  # Avoid -Inf
    }
  }

  return(loglik)
}

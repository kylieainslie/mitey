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
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{mean}: Estimated mean of the serial interval distribution (days)
#'   \item \code{sd}: Estimated standard deviation of the serial interval distribution (days)
#'   \item \code{wts}: Numeric vector of estimated component weights representing the
#'         probability that cases belong to each transmission route. Length depends
#'         on distribution choice (7 for normal, 4 for gamma)
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
#' # Example 1: Basic usage with simulated data
#' # Simulate ICC intervals from a mixed outbreak
#' set.seed(123)
#' simulated_icc <- c(
#'   rep(1, 38),   # Short intervals (co-primary cases)
#'   rep(2, 39),   #
#'   rep(3, 30),   # Medium intervals (primary-secondary)
#'   rep(4, 17),   #
#'   rep(5, 7)     # Longer intervals (higher generation)
#' )
#'
#' # Estimate serial interval parameters
#' result <- si_estim(simulated_icc, dist = "normal")
#' print(result)
#'
#' # Example 2: Using gamma distribution
#' result_gamma <- si_estim(simulated_icc, dist = "gamma")
#' print(result_gamma)
#'
#' # Example 3: Providing custom initial values
#' result_custom <- si_estim(
#'   simulated_icc,
#'   dist = "normal",
#'   init = c(10, 3),  # Start with mean=10, sd=3
#'   n = 100           # More iterations
#' )
#'
si_estim <- function(
  dat,
  n = 50,
  dist = "normal",
  init = NULL
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

  ## Vink et al. implementation code
  j <- length(dat)
  dat <- ifelse(dat == 0, 0.00001, dat)

  # Initial guesses
  if (is.null(init)) {
    mu <- mean(dat)
    sigma <- sd(dat)
  } else {
    mu <- init[1]
    sigma <- init[2]
  }
  # components depend on specified distribution
  if (dist == "normal") {
    comp_vec <- 1:7
  } else if (dist == "gamma") {
    comp_vec <- c(1, 2, 4, 6)
  }

  # E-step
  # calculate the absolute probability of interval belonging to a component

  # Iterations
  for (k in 1:n) {
    tau <- matrix(0, nrow = length(comp_vec), ncol = j)

    for (l in 1:j) {
      if (dat[l] == 0.00001) {
        for (comp in 1:length(comp_vec)) {
          tau[comp, l] <- integrate_component(
            dat[l],
            mu,
            sigma,
            comp = comp_vec[comp],
            dist = dist,
            lower = FALSE
          )
        }
      } else {
        for (comp in 1:length(comp_vec)) {
          tau[comp, l] <- integrate_component(
            dat[l],
            mu,
            sigma,
            comp = comp_vec[comp],
            dist = dist,
            lower = TRUE
          )
        }
      }
    }

    # Normalize tau
    denom <- colSums(tau)
    tau <- sweep(tau, 2, denom, "/")

    # Calculate the weights
    w <- rowSums(tau) / j

    # update parameters
    if (dist == "normal") {
      mu <- weighted.mean(dat, tau[2, ])
      sigma <- sqrt(weighted_var(dat, tau[2, ]))
    } else if (dist == "gamma") {
      # estimates for the mean and standard deviation of the primary-secondary
      # infection component
      opt <- optim(
        par = c(mu, sigma),
        wt_loglik,
        tau2 = tau[2, ],
        dat = dat,
        gr = NULL,
        method = c("BFGS"),
        hessian = FALSE
      )
      mu <- opt$par[1]
      sigma <- opt$par[2]
    }

    rtn <- list(mean = mu, sd = sigma, wts = w)
  }

  return(rtn)
}

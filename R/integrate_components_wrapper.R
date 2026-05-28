#' Compute Serial Interval Component Integrals for All Transmission Routes
#'
#' This wrapper function efficiently computes the likelihood contributions for all
#' relevant transmission route components for a given index case-to-case (ICC) interval.
#' It is a key component of the Vink method's Expectation-Maximization algorithm for
#' estimating serial interval parameters from outbreak data.
#'
#' @param d numeric; the index case-to-case (ICC) interval in days. Must be non-negative.
#' @param mu numeric; the mean of the serial interval distribution in days.
#' @param sigma numeric; the standard deviation of the serial interval distribution in days.
#' @param dist character; the assumed underlying distribution family for the serial
#'             interval. Must be either "normal" or "gamma". Defaults to "normal".
#' @param n_routes integer; the number of transmission routes to model. Must be >= 2.
#'   Defaults to 4 (Co-Primary, Primary-Secondary, Primary-Tertiary, Primary-Quaternary).
#'   For normal distribution, generates 2*n_routes - 1 components.
#'   For gamma distribution, generates n_routes components.
#'
#' @return numeric vector; integrated likelihood values for each relevant transmission
#'         route component. The length depends on the distribution and n_routes:
#' \itemize{
#'   \item Normal distribution: 2*n_routes - 1 values
#'   \item Gamma distribution: n_routes values
#' }
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' # Default 4 routes
#' integrate_components_wrapper(d = 10, mu = 15, sigma = 3, dist = "normal")
#' integrate_components_wrapper(d = 10, mu = 15, sigma = 3, dist = "gamma")
#'
#' # 5 routes
#' integrate_components_wrapper(d = 10, mu = 15, sigma = 3, dist = "normal", n_routes = 5)
#' integrate_components_wrapper(d = 10, mu = 15, sigma = 3, dist = "gamma", n_routes = 5)
#' }
#'
integrate_components_wrapper <- function(
    d,
    mu,
    sigma,
    dist    = "normal",
    n_routes = 4L
) {
  dist <- match.arg(dist, c("normal", "gamma"))

  if (!is.numeric(n_routes) || length(n_routes) != 1 ||
      is.na(n_routes) || n_routes < 2 || n_routes != floor(n_routes)) {
    stop("n_routes must be an integer >= 2.")
  }
  n_routes <- as.integer(n_routes)

  if (dist == "normal") {
    # Component 1 (CP) + pairs (2i, 2i+1) for i in 1:(n_routes - 1)
    comp_vec <- c(1L, unlist(lapply(seq_len(n_routes - 1L), function(i) c(2L*i, 2L*i + 1L))))
  } else {
    # Component 1 (CP) + even components 2i for i in 1:(n_routes - 1)
    comp_vec <- c(1L, 2L * seq_len(n_routes - 1L))
  }

  result <- sapply(comp_vec, function(comp) {
    if (d == 0) {
      integrate_component(d, mu, sigma, comp, dist = dist, lower = FALSE)
    } else {
      integrate_component(d, mu, sigma, comp, dist = dist, lower = TRUE)
    }
  })

  return(result)
}

#' Wrapper function to for integrate_components()
#'
#' This function computes the integral for all components for a given data point.
#'
#' @param d numeric; the data point.
#' @param mu numeric; the mean value.
#' @param sigma numeric; the standard deviation.
#'
#' @return A vector containing the integrated values for each component.
#' The components are:
#' \itemize{
#'   \item Component 1: Integration for CP route.
#'   \item Component 2: Integration for PS route (mean = mu, sd = sigma).
#'   \item Component 3: Integration for PS route (mean = -mu, sd = sigma).
#'   \item Component 4: Integration for PT route (mean = 2*mu, sd = sqrt(2)*sigma).
#'   \item Component 5: Integration for PT route (mean = -2*mu, sd = sqrt(2)*sigma).
#'   \item Component 6: Integration for PQ route (mean = 3*mu, sd = sqrt(3)*sigma).
#'   \item Component 7: Integration for PQ route (mean = -3*mu, sd = sqrt(3)*sigma).
#' }
#' @examples
#' # Example 1: Compute integrations for a specific data point
#' integrate_components_wrapper(1, 10, 2)
#' @export
integrate_components_wrapper <- function(d, mu, sigma) {
  result <- sapply(1:7, function(comp) {
    if (d == 0) {
      integrate_component(d, mu, sigma, comp, lower = FALSE)
    } else {
      integrate_component(d, mu, sigma, comp, lower = TRUE)
    }
  })
  return(result)
}

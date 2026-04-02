#' Compute Serial Interval Component Integrals for All Transmission Routes
#'
#' This wrapper function efficiently computes the likelihood contributions for all
#' relevant transmission route components for a given index case-to-case (ICC) interval.
#' It is a key component of the Vink method's Expectation-Maximization algorithm for
#' estimating serial interval parameters from outbreak data.
#'
#' The function handles different integration scenarios based on the distribution type
#' and ICC interval value:
#' \itemize{
#'   \item For \strong{normal distribution}: Uses all 7 components representing the full
#'         mixture of transmission routes (co-primary, primary-secondary with positive
#'         and negative components, primary-tertiary, and primary-quaternary routes)
#'   \item For \strong{gamma distribution}: Uses components 1, 2, 4, and 6 only, as
#'         the gamma distribution naturally handles only positive serial intervals,
#'         eliminating the need for negative component pairs
#'   \item For \strong{ICC interval = 0}: Uses upper integration (\code{lower = FALSE})
#'         representing the special case of simultaneous symptom onset
#'   \item For \strong{ICC interval > 0}: Uses lower integration (\code{lower = TRUE})
#'         representing the standard transmission likelihood calculation
#' }
#'
#' @param d numeric; the index case-to-case (ICC) interval in days. Represents the
#'          time difference between the symptom onset of the index case (latest case)
#'          and the current case being evaluated. Must be non-negative
#' @param mu numeric; the mean of the serial interval distribution in days. Must be
#'           positive for meaningful epidemiological interpretation
#' @param sigma numeric; the standard deviation of the serial interval distribution
#'              in days. Must be positive
#' @param dist character; the assumed underlying distribution family for the serial
#'             interval. Must be either "normal" or "gamma". Defaults to "normal".
#'             Gamma distribution is often preferred for serial intervals as it
#'             naturally restricts to positive values
#'
#' @return numeric vector; integrated likelihood values for each relevant transmission
#'         route component. The length depends on the distribution:
#' \itemize{
#'   \item Normal distribution: 7 values (components 1-7)
#'   \item Gamma distribution: 4 values (components 1, 2, 4, 6)
#' }
#'
#' @details
#' This function is primarily used internally by \code{si_estim()} as part of the
#' E-step in the EM algorithm. Each component represents a different hypothesis
#' about the transmission route:
#' \itemize{
#'   \item Component 1: Co-primary transmission (simultaneous exposure)
#'   \item Components 2-3: Primary-secondary transmission (direct transmission)
#'   \item Components 4-5: Primary-tertiary transmission (second generation)
#'   \item Components 6-7: Primary-quaternary transmission (third generation)
#' }
#'
#' For gamma distributions, components 3, 5, and 7 are omitted because the gamma
#' distribution naturally handles the asymmetry that these components would otherwise
#' model in the normal distribution case.
#'
#' @seealso \code{\link{integrate_component}}, \code{\link{si_estim}}, \code{\link{flower}},
#'          \code{\link{fupper}}, \code{\link{f0}}
#'
#' @references
#' Vink MA, Bootsma MCJ, Wallinga J (2014). Serial intervals of respiratory infectious
#' diseases: A systematic review and analysis. American Journal of Epidemiology,
#' 180(9), 865-875.
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' integrate_components_wrapper(d = 10, mu = 15, sigma = 3, dist = "normal")
#' integrate_components_wrapper(d = 10, mu = 15, sigma = 3, dist = "gamma")
#' }
#'
integrate_components_wrapper <- function(
  d,
  mu,
  sigma,
  dist = "normal"
) {
  dist <- match.arg(dist, c("normal", "gamma"))

  if (dist == "normal") {
    comp_vec <- 1:7
  } else if (dist == "gamma") {
    comp_vec <- c(1, 2, 4, 6)
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

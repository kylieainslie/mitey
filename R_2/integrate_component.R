#' Integrate Serial Interval Component Functions for Likelihood Calculation
#'
#' This function performs numerical integration of serial interval component functions
#' used in the Vink method for estimating serial interval distributions. It integrates
#' the probability density functions for different transmission routes over specified
#' intervals as part of the Expectation-Maximization algorithm.
#'
#' The function supports two integration modes:
#' \itemize{
#'   \item \code{lower = TRUE}: Integrates using \code{flower} and \code{fupper} functions
#'         over intervals \code{[d-1, d]} and \code{[d, d+1]} respectively, representing the likelihood
#'         contribution when case occurs at day d
#'   \item \code{lower = FALSE}: Integrates using \code{f0} function over interval \code{[d, d+1]},
#'         representing an alternative likelihood formulation
#' }
#'
#' The components represent different transmission routes in outbreak analysis:
#' \itemize{
#'   \item Component 1: Co-Primary (CP) transmission
#'   \item Components 2+3: Primary-Secondary (PS) transmission
#'   \item Components 4+5: Primary-Tertiary (PT) transmission
#'   \item Components 6+7: Primary-Quaternary (PQ) transmission
#' }
#'
#' @param d numeric; the index case-to-case (ICC) interval in days for which to calculate the likelihood contribution
#' @param mu numeric; the mean of the serial interval distribution in days
#' @param sigma numeric; the standard deviation of the serial interval distribution in days
#' @param comp integer; the transmission route component number (1 to 7). See Details for component definitions
#' @param dist character; the assumed underlying distribution of the serial interval.
#'             Must be either "normal" or "gamma". Defaults to "normal"
#' @param lower logical; if \code{TRUE} (default), performs integration using \code{flower} and \code{fupper} functions. If \code{FALSE}, uses \code{f0} function
#'
#' @return numeric; the integrated likelihood value for the specified component and data point. Used in the EM algorithm for serial interval estimation
#'
#' @details This function is primarily used internally by \code{si_estim()} as part of the Vink method for estimating serial interval parameters from outbreak data.
#'
#' @seealso \code{\link{flower}}, \code{\link{fupper}}, \code{\link{f0}}, \code{\link{si_estim}}
#'
#' @references
#' Vink MA, Bootsma MCJ, Wallinga J (2014). Serial intervals of respiratory infectious
#' diseases: A systematic review and analysis. American Journal of Epidemiology,
#' 180(9), 865-875.
#' @keywords internal
#' @examples
#' \dontrun{
#' integrate_component(d = 15, mu = 12, sigma = 3, comp = 2, dist = "normal", lower = TRUE)
#' integrate_component(d = 15, mu = 12, sigma = 3, comp = 2, dist = "normal", lower = FALSE)
#' }
integrate_component <- function(
  d,
  mu,
  sigma,
  comp,
  dist = c("normal", "gamma"),
  lower = TRUE
) {
  if (lower) {
    return(
      integrate(
        f = flower,
        lower = (d - 1),
        upper = d,
        r = d,
        mu = mu,
        sigma = sigma,
        comp = comp,
        dist = dist
      )[[1]] +

        integrate(
          f = fupper,
          lower = d,
          upper = (d + 1),
          r = d,
          mu = mu,
          sigma = sigma,
          comp = comp,
          dist = dist
        )[[1]]
    )
  } else {
    return(
      integrate(
        f = f0,
        lower = d,
        upper = (d + 1),
        mu = mu,
        sigma = sigma,
        comp = comp,
        dist = dist
      )[[1]]
    )
  }
}

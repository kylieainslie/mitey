#' Wrapper function to for integrate_components()
#'
#' This function computes the integral for all components for a given data point.
#'
#' @param d numeric; the data point.
#' @param mu numeric; the mean value.
#' @param sigma numeric; the standard deviation.
#' @param dist string; assumed distribution of the serial interval; take "normal" or "gamma"
#'
#' @return A vector containing the integrated values for each component.
#'
#' @examples
#' # Example 1: Compute integrations for a specific data point
#' integrate_components_wrapper(1, 10, 2)
#' @export
integrate_components_wrapper <- function(d, mu, sigma, dist= c("normal", "gamma")) {

  if(dist == "normal"){ comp_vec <- 1:7
  } else if (dist == "gamma") { comp_vec <- c(1,2,4,6)
  }

  result <- sapply(1:7, function(comp) {
    if (d == 0) {
      integrate_component(d, mu, sigma, comp, dist = dist, lower = FALSE)
    } else {
      integrate_component(d, mu, sigma, comp, dist = dist, lower = TRUE)
    }
  })
  return(result)
}

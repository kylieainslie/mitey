#' Integrate Component Function
#'
#' This function integrates the component function over a specified interval.
#'
#' @param d numeric; the data point.
#' @param mu numeric; the mean value.
#' @param sigma numeric, the standard deviation.
#' @param comp integer; the component number, accepts values 1 to 7.
#' @param dist string; assumed distribution of the serial interval; take "normal" or "gamma".
#' @param lower logical; indicates whether to integrate the lower part (TRUE) or the upper part (FALSE).
#'
#' @return The integrated value.
#' @export
integrate_component <- function(d, mu, sigma, comp, dist = c("normal", "gamma"),
                                lower = TRUE) {

  if (lower) {
    return(integrate(f = flower,
                     lower = (d - 1),
                     upper = d,
                     r = d,
                     mu = mu,
                     sigma = sigma,
                     comp = comp,
                     dist = dist
                     )[[1]] +

          integrate(f = fupper,
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
    return(integrate(f = f0,
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

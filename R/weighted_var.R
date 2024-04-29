#' Calculate weighted variance
#'
#' This function calculates the weighted variance which is then used in the estimation step (E-step) of the EM algorithm
#'
#' @param x numeric vector
#' @param w numeric vecotr of weights to be used in the caluclation of the weighted variance of x
#' @param na.rm logical; remove NA values?
#'
#' @return numeric estimate of the weighted variance of x
#' @export
#' @importFrom stats weighted.mean
#'
#' @examples
#' a <- 1:10
#' b <- seq(0.1, 1, length.out = 10)
#'
#' weighted.var(a, b)

weighted_var <- function(x, w, na.rm = FALSE) {
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  sum.w <- sum(w)
  rtn <- (sum.w*(sum(w*(x-weighted.mean(x,w))^2))) / (sum.w^2 - sum(w^2))
  return(rtn)
}

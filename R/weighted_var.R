# function for weighted variance, to be used in E-step
# weighted variance
weighted.var <- function(x, w, na.rm = FALSE) {
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  sum.w <- sum(w)
  (sum.w*(sum(w*(x-weighted.mean(x,w))^2))) / (sum.w^2 - sum(w^2))
}

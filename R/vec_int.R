#' Vectorized integration
#'
#' @param x a (numeric) vector of index case to case intervals
#' @param l lower bound for integration
#' @param u upper bound for integration
#' @param rt route of transmission; integer from 1 to 7
#' @param q character string; quantity to estimate, "zero", "lower", "upper"
#'
#' @return a vector of solutions for the integral for every value of x
#' @export
#' @importFrom stats integrate
#'
#' @examples
#' N <- 1000; hmu<-15; hsigma<-3
#' set.seed(1234)
#' CP <- fdrtool::rhalfnorm((0.2*N),theta=sqrt(pi/2)/(sqrt(2)*hsigma))
#' PS <- stats::rnorm(0.5*N,mean=hmu,sd=hsigma)
#' PT <- stats::rnorm(0.3*N,mean=2*hmu,sd=sqrt(2)*hsigma)
#' PQ <- stats::rnorm(0.1*N,mean=3*hmu,sd=sqrt(3)*hsigma)
#' icc_intervals <- round(c(CP,PS,PT,PQ))
#'
#' tmp <- vec_int(icc_intervals, l = icc_intervals, u = (icc_intervals + 1), rt = 1, q = "zero")

# vectorized method
vec_int <- Vectorize(function(x, l, u, rt, q) integrate(f = conv_tri_dist, lower = l, upper = u,
                                      sigma = hsigma, r = x, mu = hmu, route = rt, quantity = q)[[1]])





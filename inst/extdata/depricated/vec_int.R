#' Vectorized integration
#'
#' @param x a (numeric) vector of index case to case intervals
#' @param l lower bound for integration
#' @param u upper bound for integration
#' @param s standard deviation of distribution function to be integrated
#' @param m mean of distribution function to be integrated
#' @param rt route of transmission; integer from 1 to 7
#' @param q character string; quantity to estimate, "zero", "lower", "upper"
#'
#' @return a vector of solutions for the integral for every value of x
#' @export
#' @importFrom stats integrate
#'
#' @examples
#' df <-c(rep(1,38),rep(2,39),rep(3,30),rep(4,17),rep(5,7))
#'
#' vec_int(x = df, l = df, u = (df + 1), s = sd(df), m = mean(df), rt = 1, q = "zero")

# vectorized method
vec_int <- Vectorize(function(x, l, u, s, m, rt, q) integrate(f = conv_tri_dist, lower = l, upper = u,
                                      sigma = s, r = x, mu = m, route = rt, quantity = q)[[1]])







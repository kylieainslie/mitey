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
#' my_data<-c(rep(1,38),rep(2,39),rep(3,30),rep(4,17),rep(5,7))
#'
#' vec_int(x = my_data, l = my_data, u = (my_data + 1), rt = 1, q = "zero")

# vectorized method
vec_int <- Vectorize(function(x, l, u, rt, q) integrate(f = conv_tri_dist, lower = l, upper = u,
                                      sigma = hsigma, r = x, mu = hmu, route = rt, quantity = q)[[1]])





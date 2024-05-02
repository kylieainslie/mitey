#' Estimate serial interval using the method from Vink et al. (2014)
#'
#' @param dat a (numeric) vector of index case to case intervals
#' @param n integer; number of iterations for EM algorithm; defaults to n = 50
#'
#' @return vector with estimates for the mean and standard deviation of the primary-secondary infection component
#' @export
#' @importFrom stats weighted.mean
#'
#' @examples
#' my_data<-c(rep(1,38),rep(2,39),rep(3,30),rep(4,17),rep(5,7))
#'
#' si_estim(my_data)

si_estim <- function(dat, n = 50){

  j<-length(dat)

  # EM algorithm:
  # specify plausible starting values/educated guess
  my_mu<-mean(dat)
  my_sigma<-sd(dat)

  # E-step
  for(k in 1:n){

    # calculate the absolute probability of interval belonging to a component
    # the output of get_abs_prob is a 7 x length(dat matrix); each row represents
    # a route of transmission
    tmp_mat <- get_abs_prob(dat, sigma = my_sigma, mu = my_mu)

    # divide each row by the column sums
    denom <- colSums(tmp_mat)
    tau <- sweep(tmp_mat, 2, denom, "/")
    # calculate the weights for each of the components
    w <- rowSums(tau)/j

    # M-step
    # estimates for the mean and standard deviation of the primary-secondary
    # infection component can be calculated directly
    my_mu <- weighted.mean(dat,tau[2,])
    my_sigma <- sqrt(weighted_var(dat, tau[2,]))
    rtn <- c(my_mu,my_sigma)
  }

  return(rtn)
}

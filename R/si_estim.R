#' Estimate serial interval using the EM Algorithm as developed by Vink et al. (2014)
#'
#' This function estimates the serial interval using the Expectation-Maximization (EM) algorithm.
#'
#' @param dat vector; a numeric vector of index case to case intervals
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

si_estim <- function(dat, n = 50) {
  j <- length(dat)

  # Initial guesses
  mu <- mean(dat)
  sigma <- sd(dat)

  # Iterations
  for (k in 1:n) {
    tau <- matrix(0, nrow = 7, ncol = j)

    for (l in 1:j) {
      if (dat[l] == 0) {
        for (comp in 1:7) {
          tau[comp, l] <- integrate_component(dat[l], mu, sigma, comp, lower = FALSE)
        }
      } else {
        for (comp in 1:7) {
          tau[comp, l] <- integrate_component(dat[l], mu, sigma, comp, lower = TRUE)
        }
      }
    }

    # Normalize tau
    denom <- colSums(tau)
    tau <- sweep(tau, 2, denom, "/")

    # Update parameters
    w <- rowSums(tau) / j
    mu <- weighted.mean(dat, tau[2, ])
    sigma <- sqrt(weighted_var(dat, tau[2, ]))
    rtn <- c(mu, sigma)
    print(rtn)
  }

  return(rtn)
}

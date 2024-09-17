#' Estimate serial interval using the EM Algorithm as developed by Vink et al. (2014)
#'
#' This function estimates the serial interval using the Expectation-Maximization (EM) algorithm.
#'
#' @param dat vector; a numeric vector of index case to case intervals
#' @param n integer; number of iterations for EM algorithm; defaults to n = 50
#' @param dist string; assumed distribution of the serial interval; takes "normal" or "gamma"; defaults to "normal".
#' @param init numeric vector of length 2 specifying the initial values to use for the mean and standard deviation. If init= NULL, then the sample mean and sample standard deviation each divided by 4 is used.
#' @return vector with estimates for the mean and standard deviation of the primary-secondary infection component
#' @export
#' @importFrom stats weighted.mean
#' @importFrom stats optim
#'
#' @examples
#' my_data<-c(rep(1,38),rep(2,39),rep(3,30),rep(4,17),rep(5,7))
#'
#' si_estim(my_data)

si_estim <- function(dat, n = 50, dist = "normal", init = NULL) {

  j <- length(dat)
  dat <- ifelse(dat == 0, 0.00001, dat)

  # Initial guesses
  if(is.null(init)){
    mu <- mean(dat)/4
    sigma <- sd(dat)/4
  } else {
    mu <- init[1]
    sigma <- init[2]
  }
  # components depend on specified distribution
  if(dist == "normal"){ comp_vec <- 1:7
  } else if (dist == "gamma") { comp_vec <- c(1,2,4,6)
  }

  # E-step
  # calculate the absolute probability of interval belonging to a component

  # Iterations
  for (k in 1:n) {

    tau <- matrix(0, nrow = length(comp_vec), ncol = j)

    for (l in 1:j) {

      if (dat[l] == 0.00001) {
        for (comp in 1:length(comp_vec)) {
          tau[comp, l] <- integrate_component(dat[l], mu, sigma, comp = comp_vec[comp],
                                              dist = dist, lower = FALSE)
        }
      } else {
        for (comp in 1:length(comp_vec)) {
          tau[comp, l] <- integrate_component(dat[l], mu, sigma, comp = comp_vec[comp],
                                              dist = dist, lower = TRUE)
        }
      }
    }

    # Normalize tau
    denom <- colSums(tau)
    tau <- sweep(tau, 2, denom, "/")

    # Calculate the weights
    w <- rowSums(tau) / j

    # update parameters
    if (dist == "normal"){
      mu <- weighted.mean(dat, tau[2, ])
      sigma <- sqrt(weighted_var(dat, tau[2, ]))

    } else if (dist == "gamma"){
      # estimates for the mean and standard deviation of the primary-secondary
      # infection component
      opt <- optim(par = c(mu, sigma), wt_loglik, tau2 = tau[2,], dat = dat,
                   gr=NULL, method = c("BFGS"), hessian=FALSE)
      mu <- opt$par[1]
      sigma <- opt$par[2]
    }

    rtn <- list(mean = mu,
                sd = sigma,
                wts = w)
  }

  return(rtn)
}

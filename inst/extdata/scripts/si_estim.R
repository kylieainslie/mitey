#' Estimate serial interval using the method from Vink et al. (2014)
#'
#' @param dat a (numeric) vector of index case to case intervals
#'
#' @return vector with estimates for the mean and standard deviation of the primary-secondary infection component
#' @export
#' @import fdrtool
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
#' si_estim(icc_intervals)

# [1] 19.05389  9.62618

si_estim <- function(dat){
  # mixture with 7 components
  # we split the folded normal distribution for the PS, PT and PQ route into two parts
  # component 1: CP route
  # component 2+3: PS route
  # component 4+5: PT route
  # component 6+7: PQ route

  j<-length(dat)

  # EM algorithm:
  tau1<-numeric(j)    #mixture weights for (coprimary, coprimary) pairs
  tau2<-numeric(j)    #mixture weights for (primary, secondary) pairs
  tau3<-numeric(j)    #mixture weights for (secondary, primary) pairs
  tau4<-numeric(j)    #mixture weights for (primary, tertiary) pairs
  tau5<-numeric(j)    #mixture weights for (tertiary, primary) pairs
  tau6<-numeric(j)    #mixture weights for (primary, quaternary) pairs
  tau7<-numeric(j)    #mixture weights for (quaternary, primary) pairs

  # plausible starting values/educated guess
  mu<-mean(dat)
  sigma<-sd(dat)

  # number of iterations
  N<-50

  # E-step

  # calculate the absolute probability of interval belonging to a component
  # TODO vectorise
  for(k in 1:N){



      # then calculate the relative probability of a data point to belong to one of
      # the components
      dummy<-d1+d2+d3+d4+d5+d6+d7
      tau1[l]<-d1/dummy
      tau2[l]<-d2/dummy
      tau3[l]<-d3/dummy
      tau4[l]<-d4/dummy
      tau5[l]<-d5/dummy
      tau6[l]<-d6/dummy
      tau7[l]<-d7/dummy
    }

    # now calculate the weights for each of the components
    w1<-sum(tau1)/j
    w2<-sum(tau2)/j
    w3<-sum(tau3)/j
    w4<-sum(tau4)/j
    w5<-sum(tau5)/j
    w6<-sum(tau6)/j
    w7<-sum(tau7)/j

    # M-step
    # estimates for the mean and standard deviation of the primary-secondary
    # infection component can be calculated directly
    mu<-weighted.mean(dat,tau2)
    sigma<-sqrt(weighted_var(dat, tau2))
    rtn <- c(mu,sigma)
    return(rtn)
  }
}

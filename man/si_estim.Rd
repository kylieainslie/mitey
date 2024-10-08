% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/si_estim.R
\name{si_estim}
\alias{si_estim}
\title{Estimate serial interval using the EM Algorithm as developed by Vink et al. (2014)}
\usage{
si_estim(dat, n = 50, dist = "normal", init = NULL)
}
\arguments{
\item{dat}{vector; a numeric vector of index case to case intervals}

\item{n}{integer; number of iterations for EM algorithm; defaults to n = 50}

\item{dist}{string; assumed distribution of the serial interval; takes "normal" or "gamma"; defaults to "normal".}

\item{init}{numeric vector of length 2 specifying the initial values to use for the mean and standard deviation. If init= NULL, then the sample mean and sample standard deviation each divided by 4 is used.}
}
\value{
vector with estimates for the mean and standard deviation of the primary-secondary infection component
}
\description{
This function estimates the serial interval using the Expectation-Maximization (EM) algorithm.
}
\examples{
my_data<-c(rep(1,38),rep(2,39),rep(3,30),rep(4,17),rep(5,7))

si_estim(my_data)
}

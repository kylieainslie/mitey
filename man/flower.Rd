% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flower.R
\name{flower}
\alias{flower}
\title{Calculate flower for Different Components}
\usage{
flower(x, r, mu, sigma, comp, dist = "normal")
}
\arguments{
\item{x}{The value at which to evaluate the function.}

\item{r}{The value of r.}

\item{mu}{The mean value.}

\item{sigma}{The standard deviation.}

\item{comp}{The component number (1 to 7).}

\item{dist}{string; assumed distribution of the serial interval; accepts "normal" or "gamma"; defaults to "normal"}
}
\value{
The calculated value of flower.
}
\description{
This function calculates the value of flower based on the component.
}

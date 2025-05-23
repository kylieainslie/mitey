% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wallinga_lipsitch.R
\name{wallinga_lipsitch}
\alias{wallinga_lipsitch}
\title{Estimate time-varying reproduction number using Wallinga-Lipsitch method with bootstrap confidence intervals}
\usage{
wallinga_lipsitch(
  incidence,
  dates,
  si_mean,
  si_sd,
  si_dist = "gamma",
  smoothing = 0,
  bootstrap = FALSE,
  n_bootstrap = 1000,
  conf_level = 0.95,
  shift = FALSE
)
}
\arguments{
\item{incidence}{Numeric vector of daily case counts}

\item{dates}{Vector of dates corresponding to the incidence data}

\item{si_mean}{Mean of the serial interval distribution}

\item{si_sd}{Standard deviation of the serial interval distribution}

\item{si_dist}{Distribution to use for serial interval ("gamma" or "normal")}

\item{smoothing}{Window size for smoothing estimates (0 for no smoothing)}

\item{bootstrap}{Logical; whether to compute bootstrap confidence intervals}

\item{n_bootstrap}{Number of bootstrap samples to generate}

\item{conf_level}{Confidence level for intervals (0.95 = 95\% CI)}

\item{shift}{Logical; whether to shift estimates by one mean serial interval (TRUE) or not (FALSE).
When TRUE, adds a 'shifted_date' column to align case reproduction number estimates with
instantaneous reproduction number estimates for proper comparison.}
}
\value{
Data frame with columns:
\itemize{
\item date: Original dates from input
\item incidence: Daily case counts
\item R: Estimated case reproduction number
\item R_corrected: Case reproduction number with right-truncation correction
\item If bootstrap=TRUE:
\itemize{
\item R_lower, R_upper: Confidence intervals for R
\item R_corrected_lower, R_corrected_upper: Confidence intervals for R_corrected
}
\item If shift=TRUE:
\itemize{
\item shifted_date: Dates shifted forward by one mean serial interval
}
}
}
\description{
This function implements the Wallinga-Lipsitch method (2007) for estimating the time-varying
case reproduction number (R_c) from incidence data. The case reproduction number represents
the average number of secondary infections generated by cases with symptom onset at time t.
}
\details{
The method works by:
\enumerate{
\item Creating a pairwise comparison between all possible infector-infectee pairs
\item Calculating the relative likelihood that case j infected case i based on their time difference
and the serial interval distribution
\item Aggregating these likelihoods to estimate the number of secondary cases generated by each case
\item Applying right-truncation correction to account for yet-unobserved secondary cases
\item Optionally calculating bootstrap confidence intervals by resampling individual cases
\item Optionally shifting the estimates forward by one serial interval (if shift=TRUE)
to align with instantaneous reproduction number estimates
}

The function supports both gamma and normal distributions for the serial interval,
smoothing of estimates, and bootstrap-based confidence intervals.
}
\references{
Wallinga J & Lipsitch M (2007). How generation intervals shape the relationship between
growth rates and reproductive numbers. Proceedings of the Royal Society B: Biological Sciences,
274(1609), 599-604.
}

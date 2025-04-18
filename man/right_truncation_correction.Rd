% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/right_truncation_correction.R
\name{right_truncation_correction}
\alias{right_truncation_correction}
\title{Correct for right truncation in estimates of time-varying reproduction number}
\usage{
right_truncation_correction(
  t,
  T_now,
  mean_serial_interval,
  sd_serial_interval,
  distribution = "normal"
)
}
\arguments{
\item{t}{integer; time point of infection or symptom onset}

\item{T_now}{interger; current time point}

\item{mean_serial_interval}{numeric; mean of serial interval distribution}

\item{sd_serial_interval}{numeric; standard deviation of serial interval distribution}

\item{distribution}{character string; assumed serial interval distribution, takes arguments "normal" or "gamma".}
}
\value{
vector of correction weights
}
\description{
This function calculates the correction factor for each observed time of symptom onset, t. It uses the cumulative distribution function (\code{pnorm}) of the normal distribution to calculate the probability that the serial interval is less than the time lag T_now - t, where T_now is the current time.
}

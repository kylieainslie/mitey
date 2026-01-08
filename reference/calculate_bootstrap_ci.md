# Calculate Bootstrap Confidence Intervals for R Estimates

Generates bootstrap confidence intervals for reproduction number
estimates by resampling the incidence data multiple times and
calculating quantiles of the resulting R distributions.

## Usage

``` r
calculate_bootstrap_ci(
  incidence,
  si_prob,
  dates,
  si_mean,
  si_sd,
  si_dist,
  smoothing,
  n_bootstrap,
  conf_level
)
```

## Arguments

- incidence:

  numeric vector; daily case counts

- si_prob:

  numeric matrix; serial interval probability matrix

- dates:

  vector; dates corresponding to incidence data

- si_mean:

  numeric; mean of the serial interval distribution

- si_sd:

  numeric; standard deviation of the serial interval distribution

- si_dist:

  character; distribution type, either "gamma" or "normal"

- smoothing:

  integer; window size for temporal smoothing

- n_bootstrap:

  integer; number of bootstrap samples to generate

- conf_level:

  numeric; confidence level (between 0 and 1)

## Value

named list with confidence interval bounds:

- `r_lower, r_upper`: Confidence intervals for raw R estimates

- `r_corrected_lower, r_corrected_upper`: Confidence intervals for
  corrected R estimates

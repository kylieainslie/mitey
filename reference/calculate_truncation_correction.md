# Calculate Right-Truncation Correction Factors

Computes correction factors to adjust reproduction number estimates for
right-truncation bias. This bias occurs because cases near the end of
the observation period may have generated secondary cases that are not
yet observed.

## Usage

``` r
calculate_truncation_correction(dates, si_mean, si_sd, si_dist)
```

## Arguments

- dates:

  vector; dates corresponding to each case

- si_mean:

  numeric; mean of the serial interval distribution in days

- si_sd:

  numeric; standard deviation of the serial interval distribution in
  days

- si_dist:

  character; distribution type, either "gamma" or "normal"

## Value

numeric vector; correction factors for each case. Values \> 1 indicate
upward adjustment needed. Returns NA when correction would be unreliable
(probability of observation \<= 0.5)

## Examples

``` r
if (FALSE) { # \dontrun{
case_dates <- seq(as.Date("2023-01-01"), as.Date("2023-01-20"), by = "day")
calculate_truncation_correction(case_dates, si_mean = 7, si_sd = 3, si_dist = "gamma")
} # }
```

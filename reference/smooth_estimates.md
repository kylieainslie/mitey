# Apply Moving Average Smoothing to R Estimates

Applies temporal smoothing to reproduction number estimates using a
centered moving average window. Handles missing and infinite values
appropriately.

## Usage

``` r
smooth_estimates(r_estimate, window)
```

## Arguments

- r_estimate:

  numeric vector; reproduction number estimates to smooth. Can contain
  NA or infinite values

- window:

  integer; size of the smoothing window in time units. Window is
  centered around each point

## Value

numeric vector; smoothed reproduction number estimates of the same
length as input. Returns NA for points with insufficient valid
neighboring values

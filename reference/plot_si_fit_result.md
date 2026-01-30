# Plot Serial Interval Fit from si_estim Result

A convenience wrapper for
[`plot_si_fit`](https://kylieainslie.github.io/mitey/reference/plot_si_fit.md)
that accepts the output from
[`si_estim`](https://kylieainslie.github.io/mitey/reference/si_estim.md)
directly, automatically handling the weight aggregation for different
distribution types.

## Usage

``` r
plot_si_fit_result(
  si_result,
  dat,
  dist = c("normal", "gamma"),
  scaling_factor = 1
)
```

## Arguments

- si_result:

  list; the output from
  [`si_estim`](https://kylieainslie.github.io/mitey/reference/si_estim.md)
  containing mean, sd, and wts (weights) components

- dat:

  numeric vector; the index case-to-case (ICC) intervals in days used
  for estimation

- dist:

  character; the distribution family used for estimation. Must be either
  "normal" (default) or "gamma". Should match the distribution used in
  the original
  [`si_estim()`](https://kylieainslie.github.io/mitey/reference/si_estim.md)
  call

- scaling_factor:

  numeric; multiplicative factor to adjust the height of the fitted
  density curve. Defaults to 1

## Value

A `ggplot2` object showing the fitted distribution overlaid on a
histogram of the observed data

## Details

This function simplifies the plotting workflow by automatically
aggregating the component weights from
[`si_estim()`](https://kylieainslie.github.io/mitey/reference/si_estim.md)
output:

- For **normal distribution**: Aggregates 7 weights into 4 transmission
  route weights (co-primary, primary-secondary, primary-tertiary,
  primary-quaternary)

- For **gamma distribution**: Passes the first 3 weights (co-primary,
  primary-secondary, primary-tertiary); the primary-quaternary weight is
  derived internally by
  [`f_gam()`](https://kylieainslie.github.io/mitey/reference/f_gam.md)
  as `1 - w1 - w2 - w3`

## See also

[`si_estim`](https://kylieainslie.github.io/mitey/reference/si_estim.md)
for serial interval estimation,
[`plot_si_fit`](https://kylieainslie.github.io/mitey/reference/plot_si_fit.md)
for the underlying plotting function

## Examples

``` r
# Simulate some ICC interval data
set.seed(123)
icc_data <- c(
  abs(rnorm(15, mean = 0, sd = 2)),
  rnorm(40, mean = 12, sd = 3),
  rnorm(15, mean = 24, sd = 4)
)
icc_data <- round(pmax(icc_data, 0))

# Estimate serial interval
# \donttest{
result <- si_estim(icc_data, n = 50)

# Plot using the convenience wrapper
plot_si_fit_result(result, icc_data, dist = "normal")

# }
```

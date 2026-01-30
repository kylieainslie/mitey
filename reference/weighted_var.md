# Calculate Sample Weighted Variance

Computes the sample weighted variance of a numeric vector using
precision weights. This function implements the standard unbiased
weighted variance estimator commonly used in statistical applications
where observations have different precisions or levels of confidence.

## Usage

``` r
weighted_var(x, w, na.rm = FALSE)
```

## Arguments

- x:

  numeric vector; the data values for which to calculate weighted
  variance. Missing values are allowed if `na.rm = TRUE`

- w:

  numeric vector; the precision weights corresponding to each
  observation in `x`. These represent how much confidence to place in
  each measurement (higher = more trusted). Must be the same length as
  `x`. Should be non-negative

- na.rm:

  logical; if `TRUE`, missing values (both `NA` and `NaN`) are removed
  before computation. If `FALSE` (default), missing values will cause
  the function to return `NA`

## Value

numeric; the weighted sample variance. Returns `NA` if insufficient data
or if `na.rm = FALSE` and missing values are present

## Details

The weighted variance is calculated using the formula: \$\$s^2_w =
\frac{\sum w_i}{(\sum w_i)^2 - \sum w_i^2} \sum w_i (x_i -
\bar{x}\_w)^2\$\$

where \\\bar{x}\_w\\ is the weighted mean and \\w_i\\ are the weights.

This function uses precision weights (also called reliability weights),
which represent how much confidence or trust to place in each
observation, rather than frequency weights that represent how many times
to count each observation.

The denominator correction (\\(\sum w_i)^2 - \sum w_i^2\\) provides an
unbiased estimator for precision weights. Examples of precision weights
include probabilities (0-1), measurement confidence scores, or inverse
error variances.

## See also

[`weighted.mean`](https://rdrr.io/r/stats/weighted.mean.html) for
weighted mean calculation, [`var`](https://rdrr.io/r/stats/cor.html) for
unweighted sample variance

## Examples

``` r
if (FALSE) { # \dontrun{
values <- c(2.1, 3.5, 1.8, 4.2, 2.9)
weights <- c(0.8, 0.3, 0.9, 0.5, 0.7)
weighted_var(values, weights)
} # }
```

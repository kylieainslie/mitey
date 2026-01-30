# Integrate Serial Interval Component Functions for Likelihood Calculation

This function performs numerical integration of serial interval
component functions used in the Vink method for estimating serial
interval distributions. It integrates the probability density functions
for different transmission routes over specified intervals as part of
the Expectation-Maximization algorithm.

## Usage

``` r
integrate_component(
  d,
  mu,
  sigma,
  comp,
  dist = c("normal", "gamma"),
  lower = TRUE
)
```

## Arguments

- d:

  numeric; the index case-to-case (ICC) interval in days for which to
  calculate the likelihood contribution

- mu:

  numeric; the mean of the serial interval distribution in days

- sigma:

  numeric; the standard deviation of the serial interval distribution in
  days

- comp:

  integer; the transmission route component number (1 to 7). See Details
  for component definitions

- dist:

  character; the assumed underlying distribution of the serial interval.
  Must be either "normal" or "gamma". Defaults to "normal"

- lower:

  logical; if `TRUE` (default), performs integration using `flower` and
  `fupper` functions. If `FALSE`, uses `f0` function

## Value

numeric; the integrated likelihood value for the specified component and
data point. Used in the EM algorithm for serial interval estimation

## Details

The function supports two integration modes:

- `lower = TRUE`: Integrates using `flower` and `fupper` functions over
  intervals `[d-1, d]` and `[d, d+1]` respectively, representing the
  likelihood contribution when case occurs at day d

- `lower = FALSE`: Integrates using `f0` function over interval
  `[d, d+1]`, representing an alternative likelihood formulation

The components represent different transmission routes in outbreak
analysis:

- Component 1: Co-Primary (CP) transmission

- Components 2+3: Primary-Secondary (PS) transmission

- Components 4+5: Primary-Tertiary (PT) transmission

- Components 6+7: Primary-Quaternary (PQ) transmission

This function is primarily used internally by
[`si_estim()`](https://kylieainslie.github.io/mitey/reference/si_estim.md)
as part of the Vink method for estimating serial interval parameters
from outbreak data.

## References

Vink MA, Bootsma MCJ, Wallinga J (2014). Serial intervals of respiratory
infectious diseases: A systematic review and analysis. American Journal
of Epidemiology, 180(9), 865-875.

## See also

[`flower`](https://kylieainslie.github.io/mitey/reference/flower.md),
[`fupper`](https://kylieainslie.github.io/mitey/reference/fupper.md),
[`f0`](https://kylieainslie.github.io/mitey/reference/f0.md),
[`si_estim`](https://kylieainslie.github.io/mitey/reference/si_estim.md)

## Examples

``` r
if (FALSE) { # \dontrun{
integrate_component(d = 15, mu = 12, sigma = 3, comp = 2, dist = "normal", lower = TRUE)
integrate_component(d = 15, mu = 12, sigma = 3, comp = 2, dist = "normal", lower = FALSE)
} # }
```

# Compute Serial Interval Component Integrals for All Transmission Routes

This wrapper function efficiently computes the likelihood contributions
for all relevant transmission route components for a given index
case-to-case (ICC) interval. It is a key component of the Vink method's
Expectation-Maximization algorithm for estimating serial interval
parameters from outbreak data.

## Usage

``` r
integrate_components_wrapper(d, mu, sigma, dist = "normal")
```

## Arguments

- d:

  numeric; the index case-to-case (ICC) interval in days. Represents the
  time difference between the symptom onset of the index case (latest
  case) and the current case being evaluated. Must be non-negative

- mu:

  numeric; the mean of the serial interval distribution in days. Must be
  positive for meaningful epidemiological interpretation

- sigma:

  numeric; the standard deviation of the serial interval distribution in
  days. Must be positive

- dist:

  character; the assumed underlying distribution family for the serial
  interval. Must be either "normal" or "gamma". Defaults to "normal".
  Gamma distribution is often preferred for serial intervals as it
  naturally restricts to positive values

## Value

numeric vector; integrated likelihood values for each relevant
transmission route component. The length depends on the distribution:

- Normal distribution: 7 values (components 1-7)

- Gamma distribution: 4 values (components 1, 2, 4, 6)

## Details

The function handles different integration scenarios based on the
distribution type and ICC interval value:

- For **normal distribution**: Uses all 7 components representing the
  full mixture of transmission routes (co-primary, primary-secondary
  with positive and negative components, primary-tertiary, and
  primary-quaternary routes)

- For **gamma distribution**: Uses components 1, 2, 4, and 6 only, as
  the gamma distribution naturally handles only positive serial
  intervals, eliminating the need for negative component pairs

- For **ICC interval = 0**: Uses upper integration (`lower = FALSE`)
  representing the special case of simultaneous symptom onset

- For **ICC interval \> 0**: Uses lower integration (`lower = TRUE`)
  representing the standard transmission likelihood calculation

This function is primarily used internally by
[`si_estim()`](https://kylieainslie.github.io/mitey/reference/si_estim.md)
as part of the E-step in the EM algorithm. Each component represents a
different hypothesis about the transmission route:

- Component 1: Co-primary transmission (simultaneous exposure)

- Components 2-3: Primary-secondary transmission (direct transmission)

- Components 4-5: Primary-tertiary transmission (second generation)

- Components 6-7: Primary-quaternary transmission (third generation)

For gamma distributions, components 3, 5, and 7 are omitted because the
gamma distribution naturally handles the asymmetry that these components
would otherwise model in the normal distribution case.

## References

Vink MA, Bootsma MCJ, Wallinga J (2014). Serial intervals of respiratory
infectious diseases: A systematic review and analysis. American Journal
of Epidemiology, 180(9), 865-875.

## See also

[`integrate_component`](https://kylieainslie.github.io/mitey/reference/integrate_component.md),
[`si_estim`](https://kylieainslie.github.io/mitey/reference/si_estim.md),
[`flower`](https://kylieainslie.github.io/mitey/reference/flower.md),
[`fupper`](https://kylieainslie.github.io/mitey/reference/fupper.md),
[`f0`](https://kylieainslie.github.io/mitey/reference/f0.md)

## Examples

``` r
# Basic example with normal distribution
# Returns 7 component values for ICC interval of 10 days
integrate_components_wrapper(d = 10, mu = 15, sigma = 3, dist = "normal")
#> [1] 1.193978e-02 3.369814e-02 1.997254e-16 1.547844e-06 6.921256e-21
#> [6] 1.234742e-11 5.012834e-26

# Same parameters with gamma distribution
# Returns 4 component values (components 1, 2, 4, 6)
integrate_components_wrapper(d = 10, mu = 15, sigma = 3, dist = "gamma")
#> [1] 1.196761e-02 3.355068e-02 2.444605e-10 6.004790e-24

# Special case: ICC interval of 0 (simultaneous onset)
integrate_components_wrapper(d = 0, mu = 12, sigma = 2, dist = "normal")
#> [1] 2.791926e-01 1.042142e-08 1.371410e-09 1.145498e-16 1.483522e-17
#> [6] 1.434277e-24 1.847568e-25
```

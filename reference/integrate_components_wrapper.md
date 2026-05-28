# Compute Serial Interval Component Integrals for All Transmission Routes

This wrapper function efficiently computes the likelihood contributions
for all relevant transmission route components for a given index
case-to-case (ICC) interval. It is a key component of the Vink method's
Expectation-Maximization algorithm for estimating serial interval
parameters from outbreak data.

## Usage

``` r
integrate_components_wrapper(d, mu, sigma, dist = "normal", n_routes = 4L)
```

## Arguments

- d:

  numeric; the index case-to-case (ICC) interval in days. Must be
  non-negative.

- mu:

  numeric; the mean of the serial interval distribution in days.

- sigma:

  numeric; the standard deviation of the serial interval distribution in
  days.

- dist:

  character; the assumed underlying distribution family for the serial
  interval. Must be either "normal" or "gamma". Defaults to "normal".

- n_routes:

  integer; the number of transmission routes to model. Must be \>= 2.
  Defaults to 4 (Co-Primary, Primary-Secondary, Primary-Tertiary,
  Primary-Quaternary). For normal distribution, generates 2\*n_routes -
  1 components. For gamma distribution, generates n_routes components.

## Value

numeric vector; integrated likelihood values for each relevant
transmission route component. The length depends on the distribution and
n_routes:

- Normal distribution: 2\*n_routes - 1 values

- Gamma distribution: n_routes values

## Examples

``` r
if (FALSE) { # \dontrun{
# Default 4 routes
integrate_components_wrapper(d = 10, mu = 15, sigma = 3, dist = "normal")
integrate_components_wrapper(d = 10, mu = 15, sigma = 3, dist = "gamma")

# 5 routes
integrate_components_wrapper(d = 10, mu = 15, sigma = 3, dist = "normal", n_routes = 5)
integrate_components_wrapper(d = 10, mu = 15, sigma = 3, dist = "gamma", n_routes = 5)
} # }
```

# Calculate serial interval mixture density assuming underlying normal distribution

This function computes the weighted mixture density for serial intervals
based on different transmission routes in an outbreak. It implements
part of the Vink et al. (2014) method for serial interval estimation,
assuming an underlying normal distribution for the serial interval.

## Usage

``` r
f_norm(x, weights, mu, sigma, n_routes)
```

## Arguments

- x:

  quantile or vector of quantiles (time in days since index case symptom
  onset)

- weights:

  numeric vector of length n_routes - 1; probability weights for each
  transmission route starting from Co-primary. The weight for the last
  route is derived as 1 - sum(weights).

- mu:

  mean serial interval in days (can be any real number)

- sigma:

  standard deviation of serial interval in days (must be positive)

## Value

Vector of weighted density values corresponding to input quantiles x.

## Details

The function models n_routes distinct transmission routes:

- Co-primary (CP): Cases infected simultaneously from the same source

- Primary-secondary (PS): Direct transmission from index case

- Primary-tertiary (PT): Transmission through one intermediate case

- And so on up to n_routes

The weights vector must have length n_routes - 1, with the remaining
probability (1 - sum(weights)) assigned to the last route. The
transmission route distributions are parameterized as:

- Co-primary: Half-normal with scale parameter derived from sigma

- Route i (i \>= 2): Normal(i \* mu, sqrt(i) \* sigma) and Normal(-i \*
  mu, sqrt(i) \* sigma)

## References

Vink, M. A., Bootsma, M. C. J., & Wallinga, J. (2014). Serial intervals
of respiratory infectious diseases: A systematic review and analysis.
American Journal of Epidemiology, 180(9), 865-875.

## See also

[`si_estim`](https://kylieainslie.github.io/mitey/reference/si_estim.md),
[`plot_si_fit`](https://kylieainslie.github.io/mitey/reference/plot_si_fit.md),
[`f_gam`](https://kylieainslie.github.io/mitey/reference/f_gam.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# 4 routes (default behaviour)
x <- seq(0, 400, by = 1)
density_values <- f_norm(x, weights = c(0.15, 0.50, 0.25), mu = 123, sigma = 32)
plot(x, density_values, type = "l")

# 5 routes
density_values5 <- f_norm(x, weights = c(0.15, 0.50, 0.20, 0.10), mu = 123, sigma = 32)
plot(x, density_values5, type = "l")
} # }
```

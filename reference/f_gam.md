# Calculate serial interval mixture density assuming underlying gamma distribution

This function computes the weighted mixture density for serial intervals
based on different transmission routes in an outbreak. It implements
part of the Vink et al. (2014) method for serial interval estimation,
assuming an underlying gamma distribution for the serial interval.

## Usage

``` r
f_gam(x, weights, mu, sigma)
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

  mean serial interval in days (must be positive)

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
probability (1 - sum(weights)) assigned to the last route. The function
converts mean and standard deviation to gamma distribution shape (k) and
scale (theta) parameters: \$\$k = \mu^2 / \sigma^2\$\$ \$\$\theta =
\sigma^2 / \mu\$\$ Route i uses shape i \* k and scale theta.

## References

Vink, M. A., Bootsma, M. C. J., & Wallinga, J. (2014). Serial intervals
of respiratory infectious diseases: A systematic review and analysis.
American Journal of Epidemiology, 180(9), 865-875.

## See also

[`si_estim`](https://kylieainslie.github.io/mitey/reference/si_estim.md),
[`plot_si_fit`](https://kylieainslie.github.io/mitey/reference/plot_si_fit.md),
[`f_norm`](https://kylieainslie.github.io/mitey/reference/f_norm.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# 4 routes (default behaviour)
x <- seq(0.1, 30, by = 0.1)
density_values <- f_gam(x, weights = c(0.1, 0.6, 0.2), mu = 6.5, sigma = 2.8)
plot(x, density_values, type = "l")

# 5 routes
density_values5 <- f_gam(x, weights = c(0.1, 0.6, 0.15, 0.10), mu = 6.5, sigma = 2.8)
plot(x, density_values5, type = "l")
} # }
```

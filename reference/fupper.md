# Calculate fupper for Different Components

This function calculates the value of fupper based on the component.

## Usage

``` r
fupper(x, r, mu, sigma, comp, dist = "normal")
```

## Arguments

- x:

  The value at which to evaluate the function.

- r:

  The value of r.

- mu:

  The mean value.

- sigma:

  The standard deviation.

- comp:

  The component number (1 to 7).

- dist:

  string; assumed distribution of the serial interval; accepts "normal"
  or "gamma"; defaults to "normal".

## Value

The calculated value of fupper.

## Examples

``` r
# Basic example with normal distribution
# Component 2 represents primary-secondary transmission
fupper(x = 15, r = 20, mu = 12, sigma = 3, comp = 2, dist = "normal")
#> [1] 0.4839414

# Same parameters with gamma distribution
fupper(x = 15, r = 20, mu = 12, sigma = 3, comp = 2, dist = "gamma")
#> [1] 0.4131908

# Component 1 represents co-primary transmission
fupper(x = 5, r = 25, mu = 8, sigma = 2, comp = 1, dist = "normal")
#> [1] 1.241736

# Calculate for all transmission route components
x_val <- 10
r_val <- 30
mu_val <- 12
sigma_val <- 3

# Components 1-7 represent different transmission routes
sapply(1:7, function(comp) {
  fupper(x_val, r_val, mu_val, sigma_val, comp, "normal")
})
#> [1] 2.455554e-01 2.236136e+00 5.865595e-12 8.531019e-03 2.237766e-14
#> [6] 5.898261e-06 1.547169e-17
```

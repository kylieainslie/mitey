# Calculate flower for Different Components

This function calculates the value of flower based on the component.

## Usage

``` r
flower(x, r, mu, sigma, comp, dist = "normal")
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
  or "gamma"; defaults to "normal"

## Value

The calculated value of flower.

## Examples

``` r
# Basic example with normal distribution
# Component 2 represents primary-secondary transmission
flower(x = 15, r = 10, mu = 12, sigma = 3, comp = 2, dist = "normal")
#> [1] 0.4839414

# Same parameters with gamma distribution
flower(x = 15, r = 10, mu = 12, sigma = 3, comp = 2, dist = "gamma")
#> [1] 0.4131908

# Component 1 represents co-primary transmission
flower(x = 5, r = 20, mu = 8, sigma = 2, comp = 1, dist = "normal")
#> [1] -0.8278239

# Calculate for all transmission route components
x_val <- 20
r_val <- 25
mu_val <- 10
sigma_val <- 3

# Components 1-7 represent different transmission routes
sapply(1:7, function(comp) {
  flower(x_val, r_val, mu_val, sigma_val, comp, "normal")
})
#> [1] -1.124267e-05 -2.056372e-03 -1.025946e-22 -3.761264e-01 -1.876536e-20
#> [6] -4.819912e-02 -2.404707e-21
```

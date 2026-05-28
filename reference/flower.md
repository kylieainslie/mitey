# Calculate flower for Different Components

This function calculates the value of flower based on the component.

## Usage

``` r
flower(x, r, mu, sigma, comp, dist = "normal", wind = 1)
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

  integer; the component number. Component 1 is Co-Primary. Even
  components 2i are positive routes, odd components 2i+1 are negative
  routes (normal only).

- dist:

  string; assumed distribution of the serial interval; accepts "normal"
  or "gamma"; defaults to "normal"

- wind:

  The window censure interval .

## Value

The calculated value of flower.

## Examples

``` r
if (FALSE) { # \dontrun{
flower(x = 15, r = 10, mu = 12, sigma = 3, comp = 2, dist = "normal")
flower(x = 15, r = 10, mu = 12, sigma = 3, comp = 2, dist = "gamma")
} # }
```

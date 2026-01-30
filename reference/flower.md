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
if (FALSE) { # \dontrun{
flower(x = 15, r = 10, mu = 12, sigma = 3, comp = 2, dist = "normal")
flower(x = 15, r = 10, mu = 12, sigma = 3, comp = 2, dist = "gamma")
} # }
```

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
if (FALSE) { # \dontrun{
fupper(x = 15, r = 20, mu = 12, sigma = 3, comp = 2, dist = "normal")
fupper(x = 15, r = 20, mu = 12, sigma = 3, comp = 2, dist = "gamma")
} # }
```

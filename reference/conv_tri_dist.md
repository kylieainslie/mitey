# Convolution of the triangular distribution with the mixture component density (continuous case)

We split the folded normal distribution for Primary-Secondary,
Primary-Tertiary and Primary-Quaternary routes into two parts

- component 1: Co-Primary route

- component 2+3: Primary-Secondary route

- component 4+5: Primary-Tertiary route

- component 6+7: Primary-Quaternary route

## Usage

``` r
conv_tri_dist(x, sigma = sd(x), r = x, mu = mean(x), route, quantity = "zero")
```

## Arguments

- x:

  vector of index case to case intervals

- sigma:

  standard deviation of density distribution

- r:

  numeric; the observed ICC interval value used for calculating
  integration bounds in "lower" and "upper" quantity modes

- mu:

  mean of density distribution

- route:

  integer; between 1 and 7 and indicates the route of transmission.

- quantity:

  character; "zero", "lower", "upper"

## Value

vector of density draws for each value of x

## Examples

``` r
if (FALSE) { # \dontrun{
iccs <- 1:30
conv_tri_dist(x = iccs, route = 1)
} # }
```

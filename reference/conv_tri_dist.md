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

  description??

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
iccs <- 1:30
conv_tri_dist(x = iccs, route = 1)
#>  [1]  0.0000000 -0.1265320 -0.2490151 -0.3651828 -0.4729775 -0.5706108
#>  [7] -0.6566121 -0.7298632 -0.7896184 -0.8355101 -0.8675395 -0.8860538
#> [13] -0.8917126 -0.8854446 -0.8683974 -0.8418841 -0.8073273 -0.7662056
#> [19] -0.7200021 -0.6701590 -0.6180381 -0.5648898 -0.5118284 -0.4598168
#> [25] -0.4096583 -0.3619946 -0.3173109 -0.2759454 -0.2381022 -0.2038672
```

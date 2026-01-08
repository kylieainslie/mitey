# Calculate Log-Likelihood for Mixture Model

Internal function to calculate the log-likelihood of the fitted mixture
model.

## Usage

``` r
calculate_mixture_loglik(dat, mu, sigma, wts, comp_vec, dist)
```

## Arguments

- dat:

  numeric vector; the data

- mu:

  numeric; estimated mean

- sigma:

  numeric; estimated standard deviation

- wts:

  numeric vector; component weights

- comp_vec:

  integer vector; component indices

- dist:

  character; distribution type ("normal" or "gamma")

## Value

numeric; log-likelihood value

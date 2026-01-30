# Calculate Weighted Negative Log-Likelihood for Gamma Distribution Parameters

Computes the weighted negative log-likelihood for gamma distribution
parameters in the M-step of the EM algorithm for serial interval
estimation. This function is used as the objective function for
numerical optimization when the serial interval distribution is assumed
to follow a gamma distribution.

## Usage

``` r
wt_loglik(par, dat, tau2)
```

## Arguments

- par:

  numeric vector of length 2; parameters to optimize where `par[1]` is
  the mean and `par[2]` is the standard deviation of the serial interval
  distribution

- dat:

  numeric vector; index case-to-case (ICC) intervals. Zero values are
  replaced with 0.00001 to avoid gamma distribution issues at zero

- tau2:

  numeric vector; posterior probabilities (weights) that each
  observation belongs to the primary-secondary transmission component.
  These are typically derived from the E-step of the EM algorithm

## Value

numeric; negative log-likelihood value for minimization. Returns a large
penalty value (1e10) if parameters result in invalid gamma distribution
parameters (non-positive shape/scale) or non-finite likelihood values

## Details

The function converts mean and standard deviation parameters to gamma
distribution shape and scale parameters, then calculates the weighted
log-likelihood based on the posterior probabilities from the E-step.
Returns the negative log-likelihood for minimization by
[`optim`](https://rdrr.io/r/stats/optim.html).

This function is used internally by
[`si_estim`](https://kylieainslie.github.io/mitey/reference/si_estim.md)
when `dist = "gamma"`. The gamma distribution is parameterized using
shape (k) and scale (theta) parameters derived from the mean and
standard deviation:

- Shape: \\k = \mu^2 / \sigma^2\\

- Scale: \\\theta = \sigma^2 / \mu\\

The weighted log-likelihood is calculated as: \$\$\sum\_{i} \tau\_{2,i}
\log f(x_i \| k, \theta)\$\$ where \\f(x_i \| k, \theta)\\ is the gamma
probability density function and \\\tau\_{2,i}\\ are the weights from
the E-step.

## Note

This function is primarily intended for internal use within the EM
algorithm. Users typically won't call this function directly but rather
use
[`si_estim`](https://kylieainslie.github.io/mitey/reference/si_estim.md)
with `dist = "gamma"`.

## See also

[`si_estim`](https://kylieainslie.github.io/mitey/reference/si_estim.md)
for the main serial interval estimation function,
[`optim`](https://rdrr.io/r/stats/optim.html) for the optimization
routine that uses this function

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
icc_intervals <- rgamma(50, shape = 2, scale = 3)
weights <- runif(50, 0.1, 1)
wt_loglik(c(5, 4), icc_intervals, weights)
} # }
```

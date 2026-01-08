# Estimate Serial Interval Distribution Using the Vink Method

Estimates the mean and standard deviation of the serial interval
distribution from outbreak data using the Expectation-Maximization (EM)
algorithm developed by Vink et al. (2014). The serial interval is
defined as the time between symptom onset in a primary case and symptom
onset in a secondary case infected by that primary case.

## Usage

``` r
si_estim(dat, n = 50, dist = "normal", init = NULL)
```

## Arguments

- dat:

  numeric vector; index case-to-case (ICC) intervals in days. These are
  calculated as the time difference between symptom onset in each case
  and symptom onset in the index case (case with earliest onset). Must
  contain at least 2 values. Values should be non-negative in most
  epidemiological contexts, though negative values are allowed for
  normal distribution

- n:

  integer; number of EM algorithm iterations to perform. More iterations
  generally improve convergence but increase computation time. Defaults
  to 50, which is typically sufficient for convergence

- dist:

  character; the assumed parametric family for the serial interval
  distribution. Must be either:

  - `"normal"` (default): Allows negative serial intervals, uses 7
    mixture components

  - `"gamma"`: Restricts to positive serial intervals, uses 4 mixture
    components

- init:

  numeric vector of length 2; initial values for the mean and standard
  deviation to start the EM algorithm. If `NULL` (default), uses the
  sample mean and sample standard deviation of the input data. Providing
  good initial values can improve convergence, especially for
  challenging datasets

## Value

A named list containing:

- `mean`: Estimated mean of the serial interval distribution (days)

- `sd`: Estimated standard deviation of the serial interval distribution
  (days)

- `wts`: Numeric vector of estimated component weights representing the
  probability that cases belong to each transmission route. Length
  depends on distribution choice (7 for normal, 4 for gamma)

## Details

The Vink method addresses the challenge that individual transmission
pairs are typically unknown in outbreak investigations. Instead, it uses
index case-to-case (ICC) intervals - the time differences between the
case with earliest symptom onset (index case) and all other cases - to
infer the underlying serial interval distribution through a mixture
modeling approach.

**Methodological Approach:**

The method models ICC intervals as arising from a mixture of four
transmission routes:

- **Co-Primary (CP)**: Cases infected simultaneously from the same
  source

- **Primary-Secondary (PS)**: Direct transmission from index case

- **Primary-Tertiary (PT)**: Second-generation transmission

- **Primary-Quaternary (PQ)**: Third-generation transmission

The EM algorithm iteratively:

1.  **E-step**: Calculates the probability that each ICC interval
    belongs to each transmission route component

2.  **M-step**: Updates the serial interval parameters (mean, standard
    deviation) and component weights based on these probabilities

**Distribution Choice:**

- **Normal distribution**: Allows negative serial intervals (useful for
  modeling co-primary infections) and uses 7 components (positive and
  negative pairs for PS, PT, PQ routes plus CP)

- **Gamma distribution**: Restricts to positive values only, uses 4
  components (CP, PS, PT, PQ without negative pairs). Recommended when
  negative serial intervals are epidemiologically implausible

**Key Assumptions:**

- The case with earliest symptom onset is the index case

- Transmission occurs through at most 4 generations

- Serial intervals follow the specified parametric distribution

- Cases represent a single, homogeneously-mixing outbreak

**Input Data Preparation:**

To prepare ICC intervals from outbreak data:

1.  Identify the case with the earliest symptom onset date (index case)

2.  Calculate the time difference (in days) between each case's onset
    date and the index case onset date

3.  The resulting values are the ICC intervals for input to this
    function

**Convergence and Diagnostics:**

The EM algorithm typically converges within 20-50 iterations. Users
should:

- Examine the fitted distribution using
  [`plot_si_fit`](https://kylieainslie.github.io/mitey/reference/plot_si_fit.md)

- Consider alternative distribution choices if fit is poor

- Try different initial values if results seem unreasonable

- Ensure adequate sample size (generally \>20 cases recommended)

## References

Vink MA, Bootsma MCJ, Wallinga J (2014). Serial intervals of respiratory
infectious diseases: A systematic review and analysis. American Journal
of Epidemiology, 180(9), 865-875.
[doi:10.1093/aje/kwu209](https://doi.org/10.1093/aje/kwu209)

## See also

[`plot_si_fit`](https://kylieainslie.github.io/mitey/reference/plot_si_fit.md)
for diagnostic visualization,
[`integrate_component`](https://kylieainslie.github.io/mitey/reference/integrate_component.md)
for the underlying likelihood calculations

## Examples

``` r
# Example 1:Basic usage with simulated data
set.seed(123)
simulated_icc <- c(
  rep(1, 20),   # Short intervals (co-primary cases)
  rep(2, 25),   # Medium intervals (primary-secondary)
  rep(3, 15),   # Longer intervals (higher generation)
  rep(4, 8)
)

result <- si_estim(simulated_icc)

# \donttest{
# Example 2: Larger simulated outbreak, specifying distribution
large_icc <- c(
  rep(1, 38),   # Short intervals (co-primary cases)
  rep(2, 39),   #
  rep(3, 30),   # Medium intervals (primary-secondary)
  rep(4, 17),   #
  rep(5, 7),    # Longer intervals (higher generation)
  rep(6, 4),
  rep(7, 2)
)

result_normal <- si_estim(large_icc, dist = "normal")
result_gamma <- si_estim(large_icc, dist = "gamma")

# Example 3: Using custom initial values
result_custom <- si_estim(large_icc, dist = "normal", init = c(3.0, 1.5))

# Example 4: Specify iterations
result_iter <- si_estim(large_icc, n=100)

# }
```

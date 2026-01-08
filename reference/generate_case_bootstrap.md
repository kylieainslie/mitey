# Generate Bootstrap Sample of Case Incidence

Creates a bootstrap sample by resampling individual cases with
replacement, then reconstructing daily incidence counts. This maintains
the temporal distribution while introducing sampling variation for
uncertainty estimation.

## Usage

``` r
generate_case_bootstrap(incidence)
```

## Arguments

- incidence:

  numeric vector; daily case counts (non-negative integers)

## Value

numeric vector; bootstrapped daily incidence of the same length as
input. Total number of cases remains the same but their temporal
distribution varies

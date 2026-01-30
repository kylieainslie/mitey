# Calculate Serial Interval Probability Matrix

Computes a matrix of transmission probabilities between all pairs of
cases based on their time differences and the specified serial interval
distribution. Only considers epidemiologically plausible transmission
pairs (earlier to later cases).

## Usage

``` r
calculate_si_probability_matrix(day_diffs, si_mean, si_sd, si_dist)
```

## Arguments

- day_diffs:

  numeric matrix; matrix of day differences between each pair of cases,
  where element `[i,j]` represents days between case i and case j

- si_mean:

  numeric; mean of the serial interval distribution in days

- si_sd:

  numeric; standard deviation of the serial interval distribution in
  days

- si_dist:

  character; distribution type, either "gamma" or "normal"

## Value

numeric matrix; matrix of transmission probabilities where element
`[i,j]` represents the probability that case j infected case i based on
their time difference and the serial interval distribution

## Examples

``` r
if (FALSE) { # \dontrun{
# Create sample day differences matrix
dates <- as.Date(c("2023-01-01", "2023-01-03", "2023-01-05"))
day_diffs <- create_day_diff_matrix(dates)

# Calculate probability matrix
prob_matrix <- calculate_si_probability_matrix(
  day_diffs, si_mean = 7, si_sd = 3, si_dist = "gamma"
)
} # }
```

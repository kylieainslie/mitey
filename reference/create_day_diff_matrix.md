# Create Day Difference Matrix

Creates a symmetric matrix containing the time differences (in days)
between all pairs of cases based on their symptom onset dates.

## Usage

``` r
create_day_diff_matrix(dates)
```

## Arguments

- dates:

  vector; dates of symptom onset for each case. Can be Date objects or
  any format coercible to dates

## Value

numeric matrix; symmetric matrix where element `[i,j]` represents the
number of days between case i and case j (positive if i occurs after j)

## Examples

``` r
# Create day difference matrix from onset dates
onset_dates <- as.Date(c("2023-01-01", "2023-01-04", "2023-01-07", "2023-01-10"))
day_differences <- create_day_diff_matrix(onset_dates)
print(day_differences)
#>      [,1] [,2] [,3] [,4]
#> [1,]    0   -3   -6   -9
#> [2,]    3    0   -3   -6
#> [3,]    6    3    0   -3
#> [4,]    9    6    3    0
```

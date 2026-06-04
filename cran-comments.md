## R CMD check results

0 errors | 0 warnings | 1 note

* NOTE: `unable to verify current time` — this is a transient network issue
  on the check machine and is unrelated to the package.

## Resubmission notes

This is a resubmission. Version 0.4.0 adds:

* `n_routes` parameter to `si_estim()` — user-configurable number of
  transmission routes (previously fixed at 4).
* `wind` parameter to `si_estim()` — window censure interval.
* New `compare_n_routes()` function for AIC/BIC model selection across
  different numbers of transmission routes.
* Two new vignettes: quick start guide and window parameter explanation.

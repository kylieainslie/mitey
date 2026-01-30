# Changelog

## mitey 0.3.1

### New Features

- **Convenience plotting function** - Added
  [`plot_si_fit_result()`](https://kylieainslie.github.io/mitey/reference/plot_si_fit_result.md)
  which accepts
  [`si_estim()`](https://kylieainslie.github.io/mitey/reference/si_estim.md)
  output directly, handling weight aggregation automatically for both
  normal and gamma distributions

### Bug Fixes

- Fixed logical operator in
  [`weighted_var()`](https://kylieainslie.github.io/mitey/reference/weighted_var.md)
  (`|` changed to `||` for scalar comparison)
- Fixed inverted denominator formula in
  [`weighted_var()`](https://kylieainslie.github.io/mitey/reference/weighted_var.md)
  documentation
- Added defensive guards in
  [`weighted_var()`](https://kylieainslie.github.io/mitey/reference/weighted_var.md)
  for degenerate inputs (fewer than 2 values, zero weights)
- Fixed `na.rm` handling in
  [`weighted_var()`](https://kylieainslie.github.io/mitey/reference/weighted_var.md)
  to also filter `NA` weights
- Removed unreachable dead code in
  [`calculate_si_probability_matrix()`](https://kylieainslie.github.io/mitey/reference/calculate_si_probability_matrix.md)
- Fixed ggplot2 deprecation warning (`size` changed to `linewidth` in
  [`plot_si_fit()`](https://kylieainslie.github.io/mitey/reference/plot_si_fit.md))
- Fixed
  [`ggplot_build()`](https://ggplot2.tidyverse.org/reference/ggplot_build.html)
  type check in tests for newer ggplot2 versions
- Added [`match.arg()`](https://rdrr.io/r/base/match.arg.html)
  validation for `dist` parameter in
  [`integrate_components_wrapper()`](https://kylieainslie.github.io/mitey/reference/integrate_components_wrapper.md)

### API Changes

- Reduced exported functions from 18 to 5. Internal helper functions
  (`f0`, `flower`, `fupper`, `f_norm`, `f_gam`, `conv_tri_dist`,
  `integrate_component`, `integrate_components_wrapper`, `wt_loglik`,
  `weighted_var`, `calculate_si_probability_matrix`,
  `create_day_diff_matrix`, `calculate_truncation_correction`) are now
  marked as internal. They remain accessible via `mitey:::` if needed.

### Documentation

- Fixed incomplete `@param r` documentation in
  [`conv_tri_dist()`](https://kylieainslie.github.io/mitey/reference/conv_tri_dist.md)
- Added explicit `Author` and `Maintainer` fields to DESCRIPTION for R
  4.5.2 compatibility
- Updated `.Rbuildignore` to exclude `.DS_Store`, `.Rhistory`,
  `..Rcheck`, and `.claude` directories

## mitey 0.3.0

### New Features

- **Multiple restarts for EM algorithm** - Added `n_starts` parameter to
  [`si_estim()`](https://kylieainslie.github.io/mitey/reference/si_estim.md)
  to run the EM algorithm from multiple starting points and select the
  best result based on log-likelihood. This helps avoid local optima
  when initial values are far from the true parameters (fixes
  [\#7](https://github.com/kylieainslie/mitey/issues/7))

- **Convergence diagnostics** - Added `tol` parameter to
  [`si_estim()`](https://kylieainslie.github.io/mitey/reference/si_estim.md)
  for early stopping when parameters stabilize. Returns `converged` and
  `iterations` in output

- **Log-likelihood output** -
  [`si_estim()`](https://kylieainslie.github.io/mitey/reference/si_estim.md)
  now returns `loglik` for model comparison

- **CITATION file** - Added proper citation information with Zenodo DOI
  and methodology paper reference

### Improvements

- Expanded test coverage with new tests for
  [`generate_synthetic_epidemic()`](https://kylieainslie.github.io/mitey/reference/generate_synthetic_epidemic.md),
  [`plot_si_fit()`](https://kylieainslie.github.io/mitey/reference/plot_si_fit.md),
  convergence diagnostics, and multiple restarts
- Updated R version requirement to R \>= 4.0.0
- Fixed dependency table in README (ggplot2 is required, brms is
  suggested)
- Added codecov integration for test coverage reporting

## mitey 0.2.0

CRAN release: 2025-09-02

- **CRAN release** - mitey is now available on CRAN!
- Updated function names for consistency (`rt_estim()` →
  [`wallinga_lipsitch()`](https://kylieainslie.github.io/mitey/reference/wallinga_lipsitch.md))
- Enhanced documentation and examples  
- Improved package metadata and CRAN compliance
- All features from 0.1.0 maintained with refinements

## mitey 0.1.0

- Initial release
- Implements Vink et al. (2014) method for serial interval estimation
  ([`si_estim()`](https://kylieainslie.github.io/mitey/reference/si_estim.md))
- Implements Wallinga & Lipsitch (2007) method for time-varying
  reproduction number estimation (`rt_estim()`, `rt_estim_w_boot()`)
- Includes comprehensive validation against historical datasets
- Supports both Normal and Gamma serial interval distributions
- Provides bootstrap confidence intervals and visualization functions
- Developed to support epidemiological analysis in Ainslie et al. (2025)
- Complete documentation with four vignettes demonstrating usage and
  validation

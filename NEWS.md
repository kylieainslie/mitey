# mitey (development version)

# mitey 0.3.0

## New Features

* **Multiple restarts for EM algorithm** - Added `n_starts` parameter to `si_estim()`
  to run the EM algorithm from multiple starting points and select the best result
  based on log-likelihood. This helps avoid local optima when initial values are
  far from the true parameters (fixes #7)

* **Convergence diagnostics** - Added `tol` parameter to `si_estim()` for early
  stopping when parameters stabilize. Returns `converged` and `iterations` in output

* **Log-likelihood output** - `si_estim()` now returns `loglik` for model comparison

* **CITATION file** - Added proper citation information with Zenodo DOI and
  methodology paper reference

## Improvements

* Expanded test coverage with new tests for `generate_synthetic_epidemic()`,
`plot_si_fit()`, convergence diagnostics, and multiple restarts
* Updated R version requirement to R >= 4.0.0
* Fixed dependency table in README (ggplot2 is required, brms is suggested)
* Added codecov integration for test coverage reporting

# mitey 0.2.0

* **CRAN release** - mitey is now available on CRAN!
* Updated function names for consistency (`rt_estim()` → `wallinga_lipsitch()`)
* Enhanced documentation and examples  
* Improved package metadata and CRAN compliance
* All features from 0.1.0 maintained with refinements

# mitey 0.1.0

* Initial release
* Implements Vink et al. (2014) method for serial interval estimation (`si_estim()`)
* Implements Wallinga & Lipsitch (2007) method for time-varying reproduction number estimation (`rt_estim()`, `rt_estim_w_boot()`)
* Includes comprehensive validation against historical datasets
* Supports both Normal and Gamma serial interval distributions
* Provides bootstrap confidence intervals and visualization functions
* Developed to support epidemiological analysis in Ainslie et al. (2025)
* Complete documentation with four vignettes demonstrating usage and validation

# Changelog

## mitey (development version)

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

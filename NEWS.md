# mitey 0.1.0

* Initial release
* Implements Vink et al. (2014) method for serial interval estimation (`si_estim()`)
* Implements Wallinga & Lipsitch (2007) method for time-varying reproduction number estimation (`rt_estim()`, `rt_estim_w_boot()`)
* Includes comprehensive validation against historical datasets
* Supports both Normal and Gamma serial interval distributions
* Provides bootstrap confidence intervals and visualization functions
* Developed to support epidemiological analysis in Ainslie et al. (2025)
* Complete documentation with four vignettes demonstrating usage and validation

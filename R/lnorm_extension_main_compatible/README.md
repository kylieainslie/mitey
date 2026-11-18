# Lognormal Extension — Main-Compatible R Files

This folder contains **only R source files added or modified for the Lognormal extension** of the `mitey` package. The files are adapted to the upstream `main` branch, keeping its integration window (`wind`) and dynamic number of transmission routes (`n_routes`). Earlier versions are preserved separately in the repository.

## Summary of code changes

| R file | Added or modified content | Description of code changes |
|---|---|---|
| `si_estim.R` | Distribution support | Added Lognormal as a supported distribution alongside Normal and Gamma. |
| | E-step | Initialized uniform mixture weights for each EM restart and incorporated the updated weights into posterior responsibility calculations. |
| | M-step | Replaced zero-day ICC values with 1 in the Lognormal M-step to avoid numerical errors. |
| | Convergence criterion | Introduced joint convergence checks for the SI mean, standard deviation, and mixture weights across all three distributions. |
| | Output | Added optional iteration history (`return_history = TRUE`) to record parameter estimates and mixture weights across EM iterations. Disabled by default. |
| `integrate_component.R` | Lognormal integration support | Accepts `dist = "lognormal"` and passes it to the existing integration kernels. |
| `f0.R` | Zero-day component kernel | Adds Lognormal component densities to the triangular integration kernel for zero ICC interval.|
| `flower.R` | Lower-interval integration kernel | Adds Lognormal component densities to the lower triangular integration kernel. |
| `fupper.R` | Upper-interval integration kernel | Adds Lognormal component densities to the upper triangular integration kernel. |
| `f_lnorm.R` | New continuous mixture-density function | Implements CP Half-Normal, primary-secondary (PS) Lognormal, and Fenton–Wilkinson densities for higher transmission routes. |
| `plot_si_fit.R` | Lognormal fit visualization | Uses a Half-Normal approximation for CP, a Lognormal density for PS, and Lognormal approximations based on the Fenton–Wilkinson method for PT/PQ/higher transmission routes. |



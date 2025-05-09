
<!-- README.md is generated from README.Rmd. Please edit that file -->

# The `mitey` package

<!-- badges: start -->

<!-- [![R-CMD-check](https://github.com/kylieainslie/mitey/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kylieainslie/mitey/actions/workflows/R-CMD-check.yaml) -->

<!--[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/your-package)](https://CRAN.R-project.org/package=your-package) -->

<!-- [![codecov](https://codecov.io/gh/kylieainslie/mitey/branch/main/graph/badge.svg)](https://codecov.io/gh/kylieainslie/mitey) -->

<!-- [![License](https://img.shields.io/badge/license-MIT-blue.svg)](./LICENSE) -->

<!-- ![GitHub Stars](https://img.shields.io/github/stars/kylieainslie/mitey) -->

<!-- [![DOI](https://zenodo.org/badge/DOI/<DOI>.svg)](https://doi.org/<DOI>)  -->

<!--[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/mitey)](https://CRAN.R-project.org/package=mitey) -->

<!-- badges: end -->

The `mitey` package is a lightweight package designed originally as a
companion to the analyses presented by [Ainslie et
al. 2025](http://dx.doi.org/10.2139/ssrn.5184990) on scabies
transmission. However, these methods are more widely applicable than in
the context of scabies, thus the motivation behind creating the `mitey`
package was twofold and also provides flexible, documented code for
methods to estimate epidemiological quantities of interest.

Currently, `mitey` includes methods to estimate a) the mean and standard
deviation of the serial interval distribution using a maximum likelihood
framework developed by [Vink et
al. 2014](https://doi.org/10.1093/aje/kwu209) and b) the time-varying
reproduction number using the method developed by [Walling and Lipsitch
2007](https://pmc.ncbi.nlm.nih.gov/articles/PMC1766383/).

### Estimating the serial interval

The method developed by Vink et al. uses data about the time of symptom
onset with no precise information about transmission pairs and an
assumed underlying serial interval distribution (either Gaussian or
Gamma) to estimate the mean and standard deviation of the serial
interval distribution. Briefly, the method involves calculating the
index case-to-case (ICC) interval for each person, where the person with
the earliest date of symptom onset will be considered the index case.
The rest of the individuals will have an ICC interval calculated as the
number of days between their symptom onset and the index case.

## Installation

1.  Install [R](http://cran.r-project.org)

2.  Install the development version of mitey from
    [GitHub](https://github.com/kylieainslie/mitey):

``` r
# install.packages("devtools")
devtools::install_github("kylieainslie/mitey")
```

## Vignettes

- A quick start guide showing examples of how to estimate the serial
  interval and time-varying reproduction number can be found
  [here](https://kylieainslie.github.io/mitey/articles/quick_start_guide.html).

- A script that reproduces the results from [Ainslie et
  al. 2025](http://dx.doi.org/10.2139/ssrn.5184990) can be found
  [here](https://kylieainslie.github.io/mitey/articles/reproduce_results_ainslie_et_al.html).

- Validation of the method used to estimate the mean and standard
  deviation of the serial interval proposed by [Vink et
  al. 2014](https://doi.org/10.1093/aje/kwu209) can be found
  [here](https://kylieainslie.github.io/mitey/articles/code_validation_for_Vink_method.html).

- Validation of the method used to estimate the time-varying
  reproduction number proposed by [Wallinga and Lipsitch
  2007](https://pmc.ncbi.nlm.nih.gov/articles/PMC1766383/) can be found
  [here](https://kylieainslie.github.io/mitey/articles/rt_estimation_validation.html).

## Data

Several data files are stored in the repo so that the results presented
in [Ainslie et al. 2025](http://dx.doi.org/10.2139/ssrn.5184990) are
reproducible. Data files are stored in `inst/extdata/data/`. Below is a
brief description of the different files.

- `si_data.rds`
  - Description: Data on date of symptom onset for scabies outbreaks
    described by Kaburi et al., Akunzirwe et al., Tjon-Kon-Fat et al.,
    and Ariza et al. For all outbreaks except Kaburi et al. the raw data
    was not available, thus the date of symptom onset data had to be
    reconstructed using the epidemic curve provided in the manuscript.
    The original data from Kaburi et al. is also available in the `data`
    directory (`Kaburi_et_al_data_scabies.xlsx`).
  - *Source*:
    - [Kaburi et al.](https://doi.org/10.1186/s12889-019-7085-6),
    - [Akunzirwe et al.](https://doi.org/10.21203/rs.3.rs-3205380/v1),
    - [Tjon-Kon-Fat et
      al.](https://doi.org/10.1371/journal.pntd.0009485)
    - [Ariza et al.](https://doi.org/10.1007/s10096-012-1752-1).
- `scabies_data_yearly.xlsx`
  - *Description:* Annual scabies incidence per 1000 people in the
    Netherlands from 2011-2023.
  - *Source:* [Nivel](https://www.nivel.nl/nl/zorg-en-ziekte-in-cijfers)
- `scabies_data_consultation_weekly.xslx`
  - *Description:* Weekly numbers of persons consulting for scabies (per
    100,000 people) from 2011 to 2023 in the Neltherlands as diagnosed
    by general practitioners (GPs). *Note:* Individuals in institutions
    (e.g., care homes, prisons) usually have their own health care
    provider and are generally not taken into account in GP
    registrations.
  - *Source:* Nivel

<!--You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. -->

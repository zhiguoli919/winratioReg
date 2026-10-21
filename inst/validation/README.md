# Development validation scripts

These scripts document numerical validation performed during development of `winratioReg`. They are retained in the GitHub repository for provenance and are excluded from the built package through `.Rbuildignore`.

The validation suite covers:

- `local_ipcw_validation.R`: exact numerical agreement of the package Local IPCW implementation with the finalized one- and two-covariate manuscript algorithms, including bandwidths, WR grids, BLP coefficients, and lack-of-fit statistics.
- `phcopula_validation.R`: agreement of BB1, Clayton, and Gumbel copula components, likelihood fits, dependence parameters, and log-WR grids with the original manuscript-development implementation.
- `frank_gaussian_validation.R`: independent comparison of Frank and Gaussian copula CDFs, partial derivatives, densities, independence limits, and Kendall's tau with the `copula` package, plus full-fit smoke tests.
- `degree_selection_validation.R`: agreement of automatic BIC-based Bernstein-degree selection with explicit candidate fits.
- `user_workflow_validation.R`: end-to-end checks using the synthetic example data and the main exported API.

To rerun the validation suite from the repository root after `devtools::load_all()`:

```r
source("inst/validation/local_ipcw_validation.R")
source("inst/validation/phcopula_validation.R")
source("inst/validation/frank_gaussian_validation.R")
source("inst/validation/degree_selection_validation.R")
source("inst/validation/user_workflow_validation.R")
```

The file `reference/PH_BB1_recovered_functions.R` is retained only as a numerical reference for the BB1/Clayton/Gumbel validation; it is not part of the public package API.

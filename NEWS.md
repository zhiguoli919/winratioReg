# winratioReg 1.0.0

- First stable release of `winratioReg`.
- Local IPCW implementation validated numerically against the finalized manuscript algorithms in one and two covariate dimensions.
- PH--copula implementation supports BB1, Clayton, Frank, Gaussian, and Gumbel working survival copulas. BB1/Clayton/Gumbel were validated against the original manuscript-development implementation; Frank/Gaussian were validated independently against the `copula` package.
- Supports arm-specific multiplier-bootstrap inference, best linear projections, centered tests of constancy and log-linearity, fixed or separately specified Bernstein degrees, and optional BIC selection of a common Bernstein degree.
- Includes synthetic example data, documentation, a package vignette, copula-comparison tools, and reproducibility materials.
- The validated statistical estimation code is unchanged from the check-clean 0.5.0 development version; the 1.0.0 update is release/documentation/reproducibility preparation.

# winratioReg 0.5.0

- Added `winratio_example`, a small synthetic prioritized time-to-event data set for examples and vignettes.
- Added a package introduction vignette covering Local IPCW, PH--copula fitting, Bernstein-degree handling, best linear projections, log-linearity testing, and copula comparison.
- Expanded user-facing help pages and runnable examples for the main exported functions.
- Added package-level documentation and clarified the observed-data convention (`U2 <= U1`).
- Added a user-workflow validation script for the example data and primary public API.

# winratioReg 0.4.0

* Added explicit support for separate Bernstein degrees via `degree = c(m1, m2)`.
* Added optional `degree = "auto"` selection. The automatic procedure searches a user-controlled set of common degrees (default `3:7`) and selects the minimum-total-BIC degree under the specified PH--copula family.
* Automatic degree selection is performed once in the original sample; the selected degree is held fixed during multiplier-bootstrap inference.
* `compare_copulas()` deliberately requires fixed degrees so copula-family comparison is not confounded with different marginal sieve flexibility.
* Added validation code for the automatic degree selector.

# winratioReg 0.3.0

* Added Frank and Gaussian working survival copulas to `winratio_phcopula()`.
* Extended `compare_copulas()` to BB1, Clayton, Frank, Gaussian, and Gumbel while holding the PH marginal basis and Bernstein degrees fixed.
* Added Frank and Gaussian copula component tests, independence-limit checks, and an independent validation script against the `copula` package.
* Gaussian copula evaluation uses the vectorized bivariate normal distribution implementation in `pbivnorm`.

# winratioReg 0.2.0

* Added the PH--copula estimator with Bernstein-sieve baseline hazards.
* Added BB1, Clayton, and Gumbel survival-copula working families.
* Added arm-specific Exp(1) multiplier-bootstrap pointwise inference.
* Added PH--copula best linear projections and centered multiplier-bootstrap tests of constancy and log-linearity.
* Added `compare_copulas()` with log likelihood, AIC, BIC, Kendall's tau, and dependence summaries.
* Added `print()`, `summary()`, `predict()`, `plot()`, `logLik()`, and `nobs()` methods for PH--copula fits.
* Added a validation script against the original BB1/Clayton/Gumbel implementation used during manuscript development.

# winratioReg 0.1.0

* Initial development package skeleton.
* Added Local IPCW estimation for one or two continuous baseline covariates.
* Added stability-based and rule-of-thumb bandwidth options.
* Added within-arm nonparametric bootstrap pointwise inference.
* Added best linear projections and centered-bootstrap tests of constancy and log-linearity.
* Added `print()`, `summary()`, `predict()`, and `plot()` methods for Local IPCW fits.

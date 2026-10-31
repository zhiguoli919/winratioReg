# winratioReg

## Related paper

Li, Z. and Wang, X. (2026). *Flexible Estimation of Covariate-Specific Win Ratio for Prioritized Time-to-Event Outcomes.*

Preprint available on arXiv: [arXiv:2610.34567] (https://arxiv.org/abs/2610.34567).

`winratioReg` implements flexible regression methods for **covariate-specific win ratios** with prioritized time-to-event outcomes.

Version **1.0.1** provides two complementary estimators:

- **Local IPCW**, a flexible kernel-smoothed estimator that does not impose a parametric shape on the covariate-specific win ratio;
- **PH--copula**, a model-based estimator with proportional-hazards margins, Bernstein-sieve baseline cumulative hazards, and BB1, Clayton, Frank, Gaussian, or Gumbel working survival copulas.

Both workflows provide fixed-grid best linear projections. With bootstrap inference enabled, the package also provides pointwise intervals and centered tests of constancy and log-linearity.

## Installation

Install the current release from GitHub:

```r
install.packages("remotes")
remotes::install_github("zhiguoli919/winratioReg")
```

Then load the package with

```r
library(winratioReg)
```

## Data format

Use one row per subject. The primary outcome is terminal, so observation of the secondary event ends at the primary event. Consequently the observed secondary time must satisfy `U2 <= U1`.

The package includes a small **synthetic** example data set:

```r
data(winratio_example)
head(winratio_example)
```

The key columns are `treat`, `U1`, `Delta1`, `U2`, `Delta2`, and continuous baseline covariates `x1` and `x2`.

## Local IPCW

```r
fit_local <- winratio_local(
  data = winratio_example,
  arm = "treat",
  time1 = "U1", status1 = "Delta1",
  time2 = "U2", status2 = "Delta2",
  covariates = "x1",
  tau = 5,
  bandwidth = "stability",
  B = 200,
  seed = 20261021
)

summary(fit_local)
blp(fit_local)
test_loglinearity(fit_local)
plot(fit_local)
```

For point estimation only, use `B = 0`. Local IPCW currently supports one or two continuous covariates.

## PH--copula

For a quick point-estimation fit:

```r
fit_ph <- winratio_phcopula(
  data = winratio_example,
  arm = "treat",
  time1 = "U1", status1 = "Delta1",
  time2 = "U2", status2 = "Delta2",
  formula = ~ splines::ns(x1, df = 3),
  tau = 5,
  copula = "gaussian",
  degree = 6,
  B = 0,
  multistart = FALSE
)

fit_ph
blp(fit_ph)
plot(fit_ph)
```

For bootstrap inference, `parallel = TRUE` distributes the independent multiplier-bootstrap refits across a PSOCK cluster (including on Windows):

```r
fit_ph_boot <- winratio_phcopula(
  data = winratio_example,
  arm = "treat",
  time1 = "U1", status1 = "Delta1",
  time2 = "U2", status2 = "Delta2",
  formula = ~ splines::ns(x1, df = 3),
  tau = 5,
  copula = "gaussian",
  degree = 6,
  B = 20,
  multistart = FALSE,
  parallel = TRUE,
  seed = 20261021,
  verbose = TRUE
)
```
B = 20 is used here for a quick illustration. For final inference, use a larger number of bootstrap replicates (e.g., B = 200 or more). 
The multiplier bootstrap can be parallelized with parallel = TRUE.

```r
summary(fit_ph_boot)
blp(fit_ph_boot)
test_loglinearity(fit_ph_boot)
plot(fit_ph_boot)
```

The original-sample fit is serial; `parallel = TRUE` parallelizes only the bootstrap replicates. If `ncores` is omitted, the package uses up to 16 workers while leaving at least one detected logical core free. `multistart = TRUE` remains the conservative default; `FALSE` is used above to keep the example quick.

Version 1.0.1 uses a vectorized 401-point Simpson rule by default for win/loss integration. The version 1.0.0 adaptive calculation remains available with `integration_method = "adaptive"`.

The `copula` argument accepts `"bb1"`, `"clayton"`, `"frank"`, `"gaussian"`, or `"gumbel"`. BB1, Clayton, and Gumbel are positive-dependence families in the implemented parameterizations; Frank and Gaussian also allow negative dependence.

## Bernstein degree

A fixed common degree is the default:

```r
fit_common <- winratio_phcopula(..., degree = 6)
```

Separate primary/secondary degrees can be supplied explicitly:

```r
fit_separate <- winratio_phcopula(..., degree = c(5, 6))
```

Optional BIC selection chooses one **common** degree from the original sample and then holds it fixed during bootstrap refitting:

```r
fit_auto <- winratio_phcopula(
  ...,
  degree = "auto",
  degree_candidates = 3:7
)

fit_auto$degree
fit_auto$degree_selection$table
```

Automatic selection is optional; `degree = 6` remains the default.

## Comparing working copulas

Copula comparison deliberately holds the PH basis and Bernstein degrees fixed:

```r
cmp <- compare_copulas(
  data = winratio_example,
  arm = "treat",
  time1 = "U1", status1 = "Delta1",
  time2 = "U2", status2 = "Delta2",
  formula = ~ splines::ns(x1, df = 3),
  tau = 5,
  copulas = c("bb1", "clayton", "frank", "gaussian", "gumbel"),
  degree = 6
)

cmp
```

AIC/BIC are working-model comparison tools rather than proof of copula correctness. Fitting several plausible families is useful as a sensitivity analysis.

## Vignette

After installation, run:

```r
vignette("winratioReg-introduction", package = "winratioReg")
```

for a compact walkthrough of both estimators and the recommended copula-comparison workflow.

## Reproducibility

The repository contains manuscript-development scripts under `reproducibility/`. The simulation scripts are separated from the package estimation engine. Application scripts for CALGB 9633 are also provided, but the participant-level C9633 data are **not distributed** with the package or repository; authorized access to those data is required.

## Citation

To obtain the package citation in R, use:

```r
citation("winratioReg")
```

## Repository

https://github.com/zhiguoli919/winratioReg

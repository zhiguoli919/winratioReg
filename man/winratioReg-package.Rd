\name{winratioReg-package}
\alias{winratioReg}
\alias{winratioReg-package}
\docType{package}
\title{Flexible Regression Methods for Covariate-Specific Win Ratios}
\description{
The \pkg{winratioReg} package implements flexible fixed-horizon estimation and
inference for covariate-specific win ratios with prioritized time-to-event
outcomes. It provides a local inverse probability of censoring weighted (IPCW)
estimator and a model-based proportional hazards (PH)--copula estimator with
Bernstein-sieve baseline cumulative hazards.
}
\details{
The Local IPCW procedure estimates the win ratio without imposing a parametric
shape on its covariate dependence. The PH--copula procedure supports BB1,
Clayton, Frank, Gaussian, and Gumbel working survival copulas. Both workflows
provide best linear projections; bootstrap inference additionally supports
pointwise intervals and tests of constancy and log-linearity.

The observed-data convention uses one row per subject. The primary event is
terminal, so the observed secondary-event time must not exceed the observed
primary follow-up time. See \code{winratio_example} for a synthetic example.
}
\seealso{\code{winratio_local}, \code{winratio_phcopula},
\code{compare_copulas}, \code{blp}, \code{test_loglinearity}}
\keyword{package}

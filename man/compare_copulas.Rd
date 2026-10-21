\name{compare_copulas}
\alias{compare_copulas}
\alias{print.winratio_copula_comparison}
\title{Compare candidate PH--copula families}
\usage{
compare_copulas(
  data, arm, time1, status1, time2, status2, formula, tau,
  treatment = 1, control = 0,
  copulas = c("bb1", "clayton", "frank", "gaussian", "gumbel"), degree = 6L,
  sieve_upper = NULL, grid = NULL, grid_length = NULL,
  grid_probs = c(0.05, 0.95), projection_weights = NULL,
  multistart = TRUE, maxit = 2000L, integration_reltol = 1e-7,
  integration_subdivisions = 300L, verbose = FALSE
)
\method{print}{winratio_copula_comparison}(x, ...)
}
\arguments{
  \item{data, arm, time1, status1, time2, status2, formula, tau, treatment, control, degree, sieve_upper, grid, grid_length, grid_probs, projection_weights, multistart, maxit, integration_reltol, integration_subdivisions, verbose}{See \code{winratio_phcopula}. The \code{degree} argument must be explicitly specified; \code{"auto"} is not allowed here.}
  \item{copulas}{Candidate families selected from BB1, Clayton, Frank, Gaussian, and Gumbel.}
  \item{x}{A copula-comparison object.}
  \item{...}{Additional printing arguments.}
}
\value{
A data frame containing total log likelihood, AIC, BIC, convergence
information, and arm-specific dependence summaries. For the Gaussian family,
the dependence estimate is reported under \code{control_rho} and
\code{treatment_rho}. The fitted objects are retained in the \code{fits}
attribute.
}
\description{
Fits the same PH marginal basis and Bernstein degrees under several working
survival-copula families for likelihood-based comparison.
}
\details{
The Bernstein degrees are deliberately held fixed across candidate copulas.
Automatic degree selection is not performed separately within
\code{compare_copulas()}, because changing both the copula family and marginal
sieve flexibility would confound the comparison. AIC and BIC are useful
working-model comparison summaries, but they do not establish that any copula
family is correctly specified; sensitivity analysis across plausible families
is encouraged.
}
\examples{
\donttest{
data(winratio_example)
cmp <- compare_copulas(
  data = winratio_example,
  arm = "treat",
  time1 = "U1", status1 = "Delta1",
  time2 = "U2", status2 = "Delta2",
  formula = ~ x1,
  tau = 4,
  copulas = c("frank", "gaussian"),
  degree = 2,
  grid_length = 7,
  maxit = 800
)
cmp
}
}

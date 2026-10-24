\name{winratio_phcopula}
\alias{winratio_phcopula}
\title{PH--copula estimation of a covariate-specific win ratio}
\usage{
winratio_phcopula(
  data, arm, time1, status1, time2, status2, formula, tau,
  treatment = 1, control = 0,
  copula = "bb1", degree = 6L, degree_candidates = 3:7,
  sieve_upper = NULL, grid = NULL, grid_length = NULL,
  grid_probs = c(0.05, 0.95), B = 200L, conf.level = 0.95,
  projection_weights = NULL, seed = NULL, multistart = TRUE,
  maxit = 2000L, bootstrap_maxit = 1200L,
  integration_method = c("simpson", "adaptive"), integration_grid = 401L,
  integration_reltol = 1e-7, integration_subdivisions = 300L,
  parallel = FALSE, ncores = NULL, verbose = FALSE
)
}
\arguments{
  \item{data}{Data frame with one row per subject.}
  \item{arm}{Character string naming the treatment-arm column.}
  \item{time1, status1}{Character strings naming the primary terminal-event follow-up time and 0/1 event indicator.}
  \item{time2, status2}{Character strings naming the secondary-event time and 0/1 event indicator. The observation scheme requires \code{time2 <= time1}.}
  \item{formula}{One-sided formula defining the marginal PH covariate basis, such as \code{~ x1}, \code{~ x1 + x2}, or \code{~ splines::ns(x1, df = 3)}.}
  \item{tau}{Positive fixed analysis horizon.}
  \item{treatment, control}{Values of \code{arm} identifying treatment and control.}
  \item{copula}{Working survival-copula family: \code{"bb1"}, \code{"clayton"}, \code{"frank"}, \code{"gaussian"}, or \code{"gumbel"}.}
  \item{degree}{Bernstein polynomial degree. A single positive integer uses the same degree for both margins; a length-two vector gives primary and secondary degrees separately. The value \code{"auto"} selects one common degree for both margins by BIC under the specified copula family.}
  \item{degree_candidates}{Positive integer candidate degrees considered when \code{degree = "auto"}; default \code{3:7}. Ignored for fixed degrees.}
  \item{sieve_upper}{Upper endpoint for the Bernstein sieve. Defaults to the largest supplied follow-up time or \code{tau}.}
  \item{grid}{Optional raw-covariate evaluation grid.}
  \item{grid_length}{Number of grid points per raw covariate when \code{grid} is omitted.}
  \item{grid_probs}{Two probabilities defining the automatic grid range.}
  \item{B}{Number of arm-specific Exp(1) multiplier-bootstrap replicates; zero gives point estimation only.}
  \item{conf.level}{Confidence level.}
  \item{projection_weights}{Optional nonnegative weights over the evaluation grid.}
  \item{seed}{Optional multiplier-bootstrap seed.}
  \item{multistart}{Use multiple dependence-parameter starting values for the original fit.}
  \item{maxit, bootstrap_maxit}{Maximum optimizer iterations for original and bootstrap fits.}
  \item{integration_method}{Win/loss integration method. The default \code{"simpson"} uses a vectorized fixed Simpson grid; \code{"adaptive"} retains the adaptive integration used in version 1.0.0.}
  \item{integration_grid}{Odd number of time-grid points for Simpson integration; default 401. Ignored for adaptive integration.}
  \item{integration_reltol}{Relative tolerance for adaptive numerical integration of win/loss probabilities.}
  \item{integration_subdivisions}{Maximum subdivisions for adaptive integration.}
  \item{parallel}{Logical; if \code{TRUE}, multiplier-bootstrap replicates are run on a PSOCK cluster. Original arm fits remain serial.}
  \item{ncores}{Number of bootstrap workers when \code{parallel = TRUE}. If \code{NULL}, up to 16 workers are used while leaving at least one detected logical core free.}
  \item{verbose}{Print fitting and bootstrap progress.}
}
\value{
An object of class \code{winratio_phcopula}. The selected or fixed degrees are
stored in \code{degree}; automatic-selection information is stored in
\code{degree_selection}; arm-specific model fits are in \code{fit0} and
\code{fit1}; and grid-specific estimates are in \code{estimates}.
}
\description{
Fits arm-specific proportional hazards marginal models with Bernstein-sieve
baseline cumulative hazards and a working copula dependence model. The fitted
joint survival distributions are used to estimate the fixed-horizon
covariate-specific win ratio.
}
\details{
The copula is applied on the survival scale. Frank dependence uses its usual
scalar theta; Gaussian dependence uses the latent-normal correlation rho.
Pointwise uncertainty uses an arm-specific Exp(1) weighted multiplier bootstrap.
The same bootstrap provides inference for best linear projections and centered
tests of constancy and log-linearity. Bootstrap replicates can be evaluated in
parallel with \code{parallel = TRUE}; multiplier weights are generated on the
master process before parallel fitting so a fixed seed gives the same weights
in serial and parallel modes.

By default, win/loss probabilities are evaluated with a vectorized 401-point
Simpson rule. Set \code{integration_method = "adaptive"} to use the adaptive
integration strategy from version 1.0.0.

When \code{degree = "auto"}, one common Bernstein degree for the two margins is
selected by minimizing total BIC across the two arm-specific likelihoods. The
selected degree is chosen once from the original sample and held fixed in all
bootstrap refits. To compare copula families, hold the degree fixed and use
\code{compare_copulas()}.
}
\examples{
data(winratio_example)

# Quick point-estimation example under a Gaussian working copula.
fit_ph <- winratio_phcopula(
  data = winratio_example,
  arm = "treat",
  time1 = "U1", status1 = "Delta1",
  time2 = "U2", status2 = "Delta2",
  formula = ~ x1,
  tau = 4,
  copula = "gaussian",
  degree = 2,
  grid_length = 7,
  B = 0,
  maxit = 800
)
fit_ph
blp(fit_ph)
predict(fit_ph, data.frame(x1 = c(-0.5, 0, 0.5)))

\donttest{
# BIC selection of a common Bernstein degree.
fit_auto <- winratio_phcopula(
  data = winratio_example,
  arm = "treat",
  time1 = "U1", status1 = "Delta1",
  time2 = "U2", status2 = "Delta2",
  formula = ~ x1,
  tau = 4,
  copula = "gaussian",
  degree = "auto",
  degree_candidates = 2:4,
  grid_length = 7,
  B = 0
)
fit_auto$degree_selection$table
}
}

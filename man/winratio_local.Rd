\name{winratio_local}
\alias{winratio_local}
\title{Local IPCW estimation of a covariate-specific win ratio}
\usage{
winratio_local(
  data, arm, time1, status1, time2, status2, covariates, tau,
  treatment = 1, control = 0, grid = NULL, grid_length = NULL,
  grid_probs = c(0.05, 0.95), bandwidth = c("stability", "rot"),
  h_mult = c(0.7, 1, 1.3), c_ref = 1.25,
  bandwidth_tol = 1.1, B = 200L, conf.level = 0.95,
  projection_weights = NULL, seed = NULL, verbose = FALSE
)
}
\arguments{
\item{data}{A data frame with one row per subject.}
\item{arm}{Character string naming the treatment-arm column.}
\item{time1, status1}{Character strings naming the observed primary terminal-event follow-up time and 0/1 event indicator.}
\item{time2, status2}{Character strings naming the observed secondary-event time and 0/1 event indicator. The observation scheme requires \code{time2 <= time1}.}
\item{covariates}{Character vector naming one or two continuous baseline covariates.}
\item{tau}{Positive fixed analysis horizon.}
\item{treatment, control}{Values of \code{arm} identifying treatment and control.}
\item{grid}{Optional evaluation grid. A numeric vector is accepted with one covariate; use a two-column matrix or data frame with two covariates.}
\item{grid_length}{Number of grid points per covariate when the grid is generated automatically. Defaults to 21 in one dimension and 5 per dimension in two dimensions.}
\item{grid_probs}{Two probabilities defining the automatic grid range.}
\item{bandwidth}{\code{"stability"}, \code{"rot"}, or a positive numeric bandwidth with one value per covariate.}
\item{h_mult}{Candidate multipliers for stability bandwidth selection.}
\item{c_ref}{Reference-bandwidth inflation factor.}
\item{bandwidth_tol}{Tolerance above the minimum stability score.}
\item{B}{Number of within-arm ordinary bootstrap resamples; use zero for point estimation only.}
\item{conf.level}{Confidence level for pointwise intervals and bootstrap projection/test inference.}
\item{projection_weights}{Optional fixed nonnegative weights over the evaluation grid.}
\item{seed}{Optional bootstrap seed.}
\item{verbose}{Logical progress indicator.}
}
\value{
An object of class \code{winratio_local}. The component \code{estimates}
contains the grid-specific log win ratio and win ratio; \code{bandwidth}
stores the bandwidth used; \code{projection} stores constant and linear best
linear projections; and \code{bootstrap} contains bootstrap inference when
\code{B > 0}.
}
\description{
Estimates a fixed-horizon covariate-specific win ratio by local IPCW with exact
event-time win/loss integration.
}
\details{
The method uses kernel smoothing within treatment arm and an arm-specific local
Kaplan--Meier estimator of the common censoring survival function. The primary
terminal event has priority over the secondary event. The selected bandwidth is
computed from the original sample and held fixed throughout bootstrap
resampling.
}
\examples{
data(winratio_example)

fit_local <- winratio_local(
  data = winratio_example,
  arm = "treat",
  time1 = "U1", status1 = "Delta1",
  time2 = "U2", status2 = "Delta2",
  covariates = "x1",
  tau = 4,
  grid_length = 9,
  bandwidth = "rot",
  B = 0
)
fit_local
blp(fit_local)
predict(fit_local, data.frame(x1 = c(-0.5, 0, 0.5)))

\donttest{
# Bootstrap inference (increase B for substantive analyses).
fit_local_boot <- winratio_local(
  data = winratio_example,
  arm = "treat",
  time1 = "U1", status1 = "Delta1",
  time2 = "U2", status2 = "Delta2",
  covariates = "x1",
  tau = 4,
  grid_length = 9,
  bandwidth = "rot",
  B = 50,
  seed = 1
)
test_loglinearity(fit_local_boot)
}
}

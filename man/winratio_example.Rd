\name{winratio_example}
\alias{winratio_example}
\docType{data}
\title{Synthetic prioritized time-to-event example data}
\description{
A synthetic two-arm data set for illustrating \pkg{winratioReg}. It is not
derived from CALGB 9633 or any other clinical trial.
}
\usage{data(winratio_example)}
\format{
A data frame with 180 rows and 8 variables:
\describe{
  \item{\code{id}}{Synthetic subject identifier.}
  \item{\code{treat}}{Treatment indicator: 0 = control, 1 = treatment.}
  \item{\code{U1}}{Observed follow-up time for the primary terminal event.}
  \item{\code{Delta1}}{Primary-event indicator.}
  \item{\code{U2}}{Observed secondary-event time; by construction \code{U2 <= U1}.}
  \item{\code{Delta2}}{Secondary-event indicator, observed only before the primary terminal event/censoring.}
  \item{\code{x1}}{Continuous baseline covariate.}
  \item{\code{x2}}{Second continuous baseline covariate.}
}
}
\details{
The data were generated solely for software examples using proportional-hazards
margins, dependent latent event times, and right censoring. They are intended
for demonstrating syntax and package output, not for reproducing manuscript
simulation results.
}
\examples{
data(winratio_example)
head(winratio_example)
stopifnot(all(winratio_example$U2 <= winratio_example$U1 + 1e-12))
}
\keyword{datasets}

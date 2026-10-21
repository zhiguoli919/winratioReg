\name{blp}
\alias{blp}
\alias{blp.winratio_local}
\title{Best linear projection of an estimated log win ratio}
\usage{
blp(object, ...)
\method{blp}{winratio_local}(object, type = c("linear", "constant"),
  interval = c("normal", "basic"), conf.level = object$conf.level, ...)
}
\arguments{
\item{object}{A fitted object. Local IPCW and PH--copula fits are supported.}
\item{type}{\code{"linear"} for a best linear projection or \code{"constant"} for a constant projection.}
\item{interval}{Bootstrap interval type.}
\item{conf.level}{Confidence level.}
\item{...}{Additional arguments passed to the class-specific method.}
}
\value{A data frame containing projection coefficients and bootstrap uncertainty summaries. Standard errors and intervals are \code{NA} when the original fit used \code{B = 0}.}
\description{Extracts the fixed-grid best linear projection of the estimated log win ratio.}
\seealso{\code{winratio_local}, \code{winratio_phcopula}}

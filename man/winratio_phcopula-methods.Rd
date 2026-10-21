\name{winratio_phcopula-methods}
\alias{print.winratio_phcopula}
\alias{summary.winratio_phcopula}
\alias{predict.winratio_phcopula}
\alias{plot.winratio_phcopula}
\alias{blp.winratio_phcopula}
\alias{test_loglinearity.winratio_phcopula}
\alias{test_constancy.winratio_phcopula}
\alias{logLik.winratio_phcopula}
\alias{nobs.winratio_phcopula}
\title{Methods for PH--copula win-ratio fits}
\usage{
\method{print}{winratio_phcopula}(x, ...)
\method{summary}{winratio_phcopula}(object, ...)
\method{predict}{winratio_phcopula}(object, newdata, type = c("response", "link"), ...)
\method{plot}{winratio_phcopula}(x, interval = c("basic", "normal", "none"), ...)
\method{blp}{winratio_phcopula}(object, type = c("linear", "constant"),
  interval = c("normal", "basic"), conf.level = object$conf.level, ...)
\method{test_loglinearity}{winratio_phcopula}(object, ...)
\method{test_constancy}{winratio_phcopula}(object, ...)
\method{logLik}{winratio_phcopula}(object, ...)
\method{nobs}{winratio_phcopula}(object, ...)
}
\arguments{
  \item{x, object}{A fitted \code{winratio_phcopula} object.}
  \item{newdata}{Data frame or matrix of raw covariate values.}
  \item{type}{For prediction, \code{"response"} returns WR and \code{"link"} returns log WR. For \code{blp}, projection type.}
  \item{interval}{Pointwise or projection interval type.}
  \item{conf.level}{Confidence level for projection intervals.}
  \item{...}{Additional arguments.}
}
\description{Printing, prediction, plotting, best linear projection, testing, and likelihood methods for PH--copula fits.}

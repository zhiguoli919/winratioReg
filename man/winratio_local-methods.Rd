\name{winratio_local-methods}
\alias{print.winratio_local}
\alias{summary.winratio_local}
\alias{predict.winratio_local}
\alias{plot.winratio_local}
\title{Methods for Local IPCW win-ratio fits}
\usage{
\method{print}{winratio_local}(x, ...)
\method{summary}{winratio_local}(object, ...)
\method{predict}{winratio_local}(object, newdata, type = c("response", "link"), ...)
\method{plot}{winratio_local}(x, interval = c("basic", "normal", "none"), ...)
}
\arguments{
\item{x, object}{A \code{winratio_local} object.}
\item{newdata}{New covariate values for prediction.}
\item{type}{\code{"response"} returns WR; \code{"link"} returns log WR.}
\item{interval}{Pointwise interval to draw for one-dimensional fits.}
\item{...}{Additional arguments passed to plotting methods where appropriate.}
}
\description{Print, summarize, predict from, and plot a Local IPCW fit.}

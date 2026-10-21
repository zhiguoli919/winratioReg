\name{test_loglinearity}
\alias{test_loglinearity}
\alias{test_loglinearity.winratio_local}
\title{Bootstrap test of log-linearity}
\usage{
test_loglinearity(object, ...)
\method{test_loglinearity}{winratio_local}(object, ...)
}
\arguments{\item{object}{A fitted object obtained with \code{B > 0}. Local IPCW and PH--copula fits are supported.}\item{...}{Additional arguments passed to the class-specific method.}}
\value{An object of class \code{htest}.}
\description{Returns the centered-bootstrap lack-of-fit test of a linear log win-ratio specification. Local IPCW uses the ordinary within-arm bootstrap; PH--copula uses the arm-specific Exp(1) multiplier bootstrap.}
\seealso{\code{winratio_local}, \code{winratio_phcopula}}

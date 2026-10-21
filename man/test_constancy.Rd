\name{test_constancy}
\alias{test_constancy}
\alias{test_constancy.winratio_local}
\title{Bootstrap test of constancy}
\usage{
test_constancy(object, ...)
\method{test_constancy}{winratio_local}(object, ...)
}
\arguments{\item{object}{A fitted object obtained with \code{B > 0}. Local IPCW and PH--copula fits are supported.}\item{...}{Additional arguments passed to the class-specific method.}}
\value{An object of class \code{htest}.}
\description{Returns the centered-bootstrap lack-of-fit test of a constant log win-ratio specification.}
\seealso{\code{winratio_local}, \code{winratio_phcopula}}

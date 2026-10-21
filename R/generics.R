#' Best linear projection
#'
#' Extract a low-dimensional projection of an estimated log win-ratio surface.
#'
#' @param object A fitted object.
#' @param ... Additional arguments passed to a method.
#' @export
blp <- function(object, ...) {
  UseMethod("blp")
}

#' Test log-linearity
#'
#' Test whether an estimated covariate-specific log win ratio is adequately
#' represented by a linear function of the supplied covariates.
#'
#' @param object A fitted object.
#' @param ... Additional arguments passed to a method.
#' @export
test_loglinearity <- function(object, ...) {
  UseMethod("test_loglinearity")
}

#' Test constancy
#'
#' Test whether an estimated covariate-specific log win ratio is constant over
#' the evaluation grid.
#'
#' @param object A fitted object.
#' @param ... Additional arguments passed to a method.
#' @export
test_constancy <- function(object, ...) {
  UseMethod("test_constancy")
}

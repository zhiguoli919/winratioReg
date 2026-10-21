.local_trapz_weights <- function(x) {
  J <- length(x)
  if (J < 2L) stop("At least two grid points are required.", call. = FALSE)
  w <- numeric(J)
  w[1L] <- (x[2L] - x[1L]) / 2
  w[J] <- (x[J] - x[J - 1L]) / 2
  if (J > 2L) {
    for (j in 2L:(J - 1L)) w[j] <- (x[j + 1L] - x[j - 1L]) / 2
  }
  w / sum(w)
}

.local_projection_weights <- function(grid, supplied = NULL) {
  if (!is.null(supplied)) {
    w <- as.numeric(supplied)
    if (length(w) != nrow(grid) || any(!is.finite(w)) || any(w < 0) || sum(w) <= 0) {
      stop("projection_weights must be nonnegative finite values with one weight per grid point.", call. = FALSE)
    }
    return(w / sum(w))
  }
  if (ncol(grid) == 1L) .local_trapz_weights(grid[, 1L]) else rep(1 / nrow(grid), nrow(grid))
}

.local_project <- function(y, grid, w, type = c("linear", "constant")) {
  type <- match.arg(type)
  if (any(!is.finite(y))) return(NULL)
  if (type == "constant") {
    coef <- sum(w * y)
    fitted <- rep(coef, length(y))
    return(list(coefficients = coef, fitted = fitted, resid = y - fitted))
  }

  X <- cbind(1, grid)
  XtWX <- crossprod(X, w * X)
  if (qr(XtWX)$rank < ncol(X)) return(NULL)
  beta <- solve(XtWX, crossprod(X, w * y))
  fitted <- as.vector(X %*% beta)
  list(coefficients = as.numeric(beta), fitted = fitted, resid = y - fitted)
}

.local_lack_of_fit <- function(y, grid, w, type) {
  pr <- .local_project(y, grid, w, type)
  if (is.null(pr)) return(NULL)
  list(D = sum(w * pr$resid^2), fitted = pr$fitted, resid = pr$resid,
       coefficients = pr$coefficients)
}

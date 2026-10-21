#' @export
print.winratio_local <- function(x, ...) {
  cat("Local IPCW covariate-specific win-ratio fit\n")
  cat("  Horizon tau:", x$tau, "\n")
  cat("  Covariates:", paste(x$covariates, collapse = ", "), "\n")
  cat("  N:", x$data_info$n, "(treatment", x$data_info$n_treatment,
      ", control", x$data_info$n_control, ")\n")
  cat("  Grid points:", nrow(x$grid), "\n")
  cat("  Bandwidth:", paste(signif(x$bandwidth, 5), collapse = ", "), "\n")
  if (x$B > 0L) cat("  Bootstrap:", x$bootstrap$n_success, "of", x$B, "successful\n")
  invisible(x)
}

#' @export
summary.winratio_local <- function(object, ...) {
  print(object)
  cat("\nBest linear projection of log WR:\n")
  print(blp(object, type = "linear"), row.names = FALSE)
  if (object$B > 0L) {
    cat("\nLog-linearity test:\n")
    print(test_loglinearity(object))
    cat("\nConstancy test:\n")
    print(test_constancy(object))
  }
  invisible(object)
}

#' @export
blp.winratio_local <- function(object, type = c("linear", "constant"),
                               interval = c("normal", "basic"),
                               conf.level = object$conf.level, ...) {
  type <- match.arg(type)
  interval <- match.arg(interval)
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      conf.level <= 0 || conf.level >= 1) stop("Invalid conf.level.", call. = FALSE)

  pr <- object$projection[[type]]
  est <- as.numeric(pr$coefficients)
  terms <- if (type == "constant") "(Constant)" else c("(Intercept)", object$covariates)
  se <- lower <- upper <- rep(NA_real_, length(est))

  if (!is.null(object$bootstrap)) {
    bmat <- if (type == "constant") matrix(object$bootstrap$projection$constant, ncol = 1L) else object$bootstrap$projection$linear
    for (j in seq_along(est)) {
      vals <- bmat[, j]
      vals <- vals[is.finite(vals)]
      if (length(vals) >= 2L) {
        se[j] <- stats::sd(vals)
        alpha <- 1 - conf.level
        if (interval == "normal") {
          z <- stats::qnorm(1 - alpha / 2)
          lower[j] <- est[j] - z * se[j]
          upper[j] <- est[j] + z * se[j]
        } else {
          diff <- vals - est[j]
          qq <- stats::quantile(diff, probs = c(alpha / 2, 1 - alpha / 2), names = FALSE)
          lower[j] <- est[j] - qq[2L]
          upper[j] <- est[j] - qq[1L]
        }
      }
    }
  }
  data.frame(term = terms, estimate = est, SE = se, lower = lower, upper = upper,
             row.names = NULL, check.names = FALSE)
}

#' @export
test_loglinearity.winratio_local <- function(object, ...) {
  if (is.null(object$bootstrap)) stop("Bootstrap inference is unavailable; refit with B > 0.", call. = FALSE)
  z <- object$bootstrap$tests$loglinearity
  structure(list(
    statistic = c(T = z$statistic),
    parameter = c(B_success = z$n_boot_success),
    p.value = z$p.value,
    method = "Centered bootstrap test of log-linearity for Local IPCW log win ratio",
    data.name = paste(object$covariates, collapse = ", ")
  ), class = "htest")
}

#' @export
test_constancy.winratio_local <- function(object, ...) {
  if (is.null(object$bootstrap)) stop("Bootstrap inference is unavailable; refit with B > 0.", call. = FALSE)
  z <- object$bootstrap$tests$constancy
  structure(list(
    statistic = c(T = z$statistic),
    parameter = c(B_success = z$n_boot_success),
    p.value = z$p.value,
    method = "Centered bootstrap test of constancy for Local IPCW log win ratio",
    data.name = paste(object$covariates, collapse = ", ")
  ), class = "htest")
}

#' @export
predict.winratio_local <- function(object, newdata, type = c("response", "link"), ...) {
  type <- match.arg(type)
  if (!is.data.frame(newdata) && !is.matrix(newdata)) {
    stop("newdata must be a data frame or matrix.", call. = FALSE)
  }
  if (is.data.frame(newdata)) {
    missing_vars <- setdiff(object$covariates, names(newdata))
    if (length(missing_vars) > 0L) stop("newdata is missing: ", paste(missing_vars, collapse = ", "), call. = FALSE)
    G <- as.matrix(newdata[, object$covariates, drop = FALSE])
  } else {
    G <- as.matrix(newdata)
    if (ncol(G) != object$dimension) stop("newdata has the wrong number of columns.", call. = FALSE)
  }
  storage.mode(G) <- "double"
  if (any(!is.finite(G))) stop("newdata must contain finite numeric covariates.", call. = FALSE)
  wr <- .local_estimate_grid(object$analysis_data$dat1, object$analysis_data$dat0,
                             G, object$bandwidth, object$tau)
  if (type == "link") log(wr) else wr
}

#' @export
plot.winratio_local <- function(x, interval = c("basic", "normal", "none"), ...) {
  interval <- match.arg(interval)
  if (x$dimension == 1L) {
    xx <- x$grid[, 1L]
    yy <- x$estimates$WR
    ylim <- range(c(yy, 1), finite = TRUE)
    if (interval == "basic" && any(is.finite(x$estimates$Basic_lower))) {
      ylim <- range(c(ylim, x$estimates$Basic_lower, x$estimates$Basic_upper), finite = TRUE)
    }
    if (interval == "normal" && any(is.finite(x$estimates$Normal_lower))) {
      ylim <- range(c(ylim, x$estimates$Normal_lower, x$estimates$Normal_upper), finite = TRUE)
    }
    graphics::plot(xx, yy, type = "l", lwd = 2, ylim = ylim,
                   xlab = x$covariates[1L], ylab = "Win ratio", ...)
    if (interval == "basic" && any(is.finite(x$estimates$Basic_lower))) {
      graphics::lines(xx, x$estimates$Basic_lower, lty = 2)
      graphics::lines(xx, x$estimates$Basic_upper, lty = 2)
    }
    if (interval == "normal" && any(is.finite(x$estimates$Normal_lower))) {
      graphics::lines(xx, x$estimates$Normal_lower, lty = 2)
      graphics::lines(xx, x$estimates$Normal_upper, lty = 2)
    }
    graphics::abline(h = 1, lty = 3)
    return(invisible(x))
  }

  xs <- sort(unique(x$grid[, 1L]))
  ys <- sort(unique(x$grid[, 2L]))
  if (length(xs) * length(ys) != nrow(x$grid)) {
    stop("Two-dimensional plotting currently requires a complete rectangular grid.", call. = FALSE)
  }
  z <- matrix(NA_real_, nrow = length(xs), ncol = length(ys))
  for (i in seq_len(nrow(x$grid))) {
    ix <- match(x$grid[i, 1L], xs)
    iy <- match(x$grid[i, 2L], ys)
    z[ix, iy] <- x$estimates$WR[i]
  }
  graphics::filled.contour(xs, ys, z,
                           xlab = x$covariates[1L], ylab = x$covariates[2L], ...)
  invisible(x)
}

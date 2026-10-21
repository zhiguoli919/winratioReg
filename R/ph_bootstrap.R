.ph_bootstrap <- function(dat1, dat0, fit1, fit0, X_grid, raw_grid, tau,
                          weights, B, conf.level, seed, n_total,
                          m1, m2, maxit, rel.tol, subdivisions, verbose) {
  if (!is.null(seed)) set.seed(seed)
  G <- nrow(X_grid)
  boot_logWR <- matrix(NA_real_, nrow = B, ncol = G)
  success <- logical(B)
  every <- max(1L, floor(B / 10L))

  logWR_hat <- .ph_logwr_grid(fit1, fit0, X_grid, tau,
                              rel.tol = rel.tol, subdivisions = subdivisions)

  for (b in seq_len(B)) {
    w1 <- stats::rexp(nrow(dat1))
    w0 <- stats::rexp(nrow(dat0))
    w1 <- w1 / mean(w1)
    w0 <- w0 / mean(w0)

    f1b <- try(.ph_fit_arm(
      dat1, family = fit1$family, m1 = m1, m2 = m2,
      weights = w1, start = fit1$par_unconstr,
      multistart = FALSE, maxit = maxit
    ), silent = TRUE)
    f0b <- try(.ph_fit_arm(
      dat0, family = fit0$family, m1 = m1, m2 = m2,
      weights = w0, start = fit0$par_unconstr,
      multistart = FALSE, maxit = maxit
    ), silent = TRUE)

    if (!inherits(f1b, "try-error") && !inherits(f0b, "try-error") &&
        f1b$convergence == 0L && f0b$convergence == 0L) {
      vals <- try(.ph_logwr_grid(f1b, f0b, X_grid, tau,
                                 rel.tol = rel.tol, subdivisions = subdivisions),
                  silent = TRUE)
      if (!inherits(vals, "try-error") && all(is.finite(vals))) {
        boot_logWR[b, ] <- vals
        success[b] <- TRUE
      }
    }

    if (isTRUE(verbose) && (b %% every == 0L || b == B)) {
      message("PH--copula bootstrap ", b, "/", B,
              "; successful curves: ", sum(success))
    }
  }

  n_success <- sum(success)
  required <- if (B >= 20L) max(20L, ceiling(0.80 * B)) else max(2L, ceiling(0.80 * B))
  if (n_success < 2L) stop("Fewer than two PH--copula bootstrap curves were successful.", call. = FALSE)
  if (n_success < required) {
    warning("Only ", n_success, " of ", B,
            " PH--copula bootstrap curves were successful; test p-values are set to NA.",
            call. = FALSE)
  }

  alpha <- 1 - conf.level
  se <- apply(boot_logWR, 2L, stats::sd, na.rm = TRUE)
  zcrit <- stats::qnorm(1 - alpha / 2)
  normal_lower <- logWR_hat - zcrit * se
  normal_upper <- logWR_hat + zcrit * se

  boot_diff <- sweep(boot_logWR, 2L, logWR_hat, FUN = "-")
  qlo <- apply(boot_diff, 2L, stats::quantile, probs = alpha / 2,
               na.rm = TRUE, names = FALSE)
  qhi <- apply(boot_diff, 2L, stats::quantile, probs = 1 - alpha / 2,
               na.rm = TRUE, names = FALSE)
  basic_lower <- logWR_hat - qhi
  basic_upper <- logWR_hat - qlo

  obs_const <- .local_lack_of_fit(logWR_hat, raw_grid, weights, "constant")
  obs_linear <- .local_lack_of_fit(logWR_hat, raw_grid, weights, "linear")
  scale_test <- n_total
  T_const <- rep(NA_real_, B)
  T_linear <- rep(NA_real_, B)
  const_coef <- rep(NA_real_, B)
  linear_coef <- matrix(NA_real_, nrow = B, ncol = ncol(raw_grid) + 1L)

  for (b in which(success)) {
    log_b <- boot_logWR[b, ]
    diff_b <- log_b - logWR_hat
    tc <- .local_lack_of_fit(diff_b, raw_grid, weights, "constant")
    tl <- .local_lack_of_fit(diff_b, raw_grid, weights, "linear")
    if (!is.null(tc)) T_const[b] <- scale_test * tc$D
    if (!is.null(tl)) T_linear[b] <- scale_test * tl$D
    pc <- .local_project(log_b, raw_grid, weights, "constant")
    pl <- .local_project(log_b, raw_grid, weights, "linear")
    if (!is.null(pc)) const_coef[b] <- pc$coefficients
    if (!is.null(pl)) linear_coef[b, ] <- pl$coefficients
  }

  make_test <- function(obs, boot) {
    vals <- boot[is.finite(boot)]
    Tobs <- scale_test * obs$D
    if (n_success < required || length(vals) < required) {
      return(list(statistic = Tobs, critical = NA_real_, p.value = NA_real_,
                  n_boot_success = length(vals), scale = scale_test))
    }
    list(
      statistic = Tobs,
      critical = as.numeric(stats::quantile(vals, probs = conf.level, names = FALSE)),
      p.value = mean(vals >= Tobs),
      n_boot_success = length(vals),
      scale = scale_test
    )
  }

  list(
    logWR = boot_logWR,
    success = success,
    n_success = n_success,
    se = se,
    normal_lower = normal_lower,
    normal_upper = normal_upper,
    basic_lower = basic_lower,
    basic_upper = basic_upper,
    projection = list(constant = const_coef, linear = linear_coef),
    tests = list(
      constancy = make_test(obs_const, T_const),
      loglinearity = make_test(obs_linear, T_linear)
    )
  )
}

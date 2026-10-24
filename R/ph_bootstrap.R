.ph_bootstrap_one <- function(dat1, dat0, fit1, fit0, X_grid, tau,
                              m1, m2, maxit, w1, w0,
                              integration_method, integration_grid,
                              rel.tol, subdivisions) {
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

  if (inherits(f1b, "try-error") || inherits(f0b, "try-error") ||
      is.null(f1b) || is.null(f0b) ||
      f1b$convergence != 0L || f0b$convergence != 0L) {
    return(list(success = FALSE, logWR = rep(NA_real_, nrow(X_grid))))
  }

  vals <- try(.ph_logwr_grid(
    f1b, f0b, X_grid, tau,
    method = integration_method, n_grid = integration_grid,
    rel.tol = rel.tol, subdivisions = subdivisions
  ), silent = TRUE)

  if (inherits(vals, "try-error") || length(vals) != nrow(X_grid) ||
      any(!is.finite(vals))) {
    return(list(success = FALSE, logWR = rep(NA_real_, nrow(X_grid))))
  }
  list(success = TRUE, logWR = as.numeric(vals))
}

.ph_bootstrap <- function(dat1, dat0, fit1, fit0, X_grid, raw_grid, tau,
                          weights, B, conf.level, seed, n_total,
                          m1, m2, maxit, logWR_hat,
                          integration_method, integration_grid,
                          rel.tol, subdivisions, verbose,
                          parallel = FALSE, ncores = NULL) {
  if (!is.null(seed)) set.seed(seed)
  G <- nrow(X_grid)
  boot_logWR <- matrix(NA_real_, nrow = B, ncol = G)
  success <- logical(B)
  every <- max(1L, floor(B / 10L))

  ## Generate multiplier weights on the master process in the same order as
  ## the original serial implementation. This keeps seeded serial and parallel
  ## runs reproducible relative to one another.
  boot_weights <- vector("list", B)
  for (b in seq_len(B)) {
    w1 <- stats::rexp(nrow(dat1))
    w0 <- stats::rexp(nrow(dat0))
    boot_weights[[b]] <- list(w1 = w1 / mean(w1), w0 = w0 / mean(w0))
  }

  run_serial <- function() {
    ans <- vector("list", B)
    for (b in seq_len(B)) {
      ans[[b]] <- .ph_bootstrap_one(
        dat1, dat0, fit1, fit0, X_grid, tau,
        m1, m2, maxit,
        boot_weights[[b]]$w1, boot_weights[[b]]$w0,
        integration_method, integration_grid, rel.tol, subdivisions
      )
      if (isTRUE(verbose) && (b %% every == 0L || b == B)) {
        n_ok <- sum(vapply(ans[seq_len(b)], function(z) isTRUE(z$success), logical(1)))
        message("PH--copula bootstrap ", b, "/", B,
                "; successful curves: ", n_ok)
      }
    }
    ans
  }

  use_parallel <- isTRUE(parallel) && B > 1L
  if (use_parallel) {
    if (is.null(ncores)) {
      available <- parallel::detectCores(logical = TRUE)
      if (!is.finite(available) || available < 2L) available <- 2L
      ncores <- max(1L, min(B, 16L, available - 1L))
    }
    ncores <- max(1L, min(B, as.integer(ncores)))
    use_parallel <- ncores > 1L
  }

  if (use_parallel) {
    if (isTRUE(verbose)) {
      message("PH--copula bootstrap: running ", B,
              " replicates on ", ncores, " parallel workers...")
    }
    cl <- parallel::makePSOCKcluster(ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, loadNamespace("winratioReg"))
    parallel::clusterExport(
      cl,
      varlist = c(
        "dat1", "dat0", "fit1", "fit0", "X_grid", "tau",
        "m1", "m2", "maxit", "boot_weights", "integration_method",
        "integration_grid", "rel.tol", "subdivisions"
      ),
      envir = environment()
    )
    ans <- parallel::parLapply(cl, seq_len(B), function(b) {
      utils::getFromNamespace(".ph_bootstrap_one", "winratioReg")(
        dat1, dat0, fit1, fit0, X_grid, tau,
        m1, m2, maxit,
        boot_weights[[b]]$w1, boot_weights[[b]]$w0,
        integration_method, integration_grid, rel.tol, subdivisions
      )
    })
    parallel::stopCluster(cl)
    on.exit(NULL, add = FALSE)
    if (isTRUE(verbose)) {
      message("PH--copula bootstrap complete; successful curves: ",
              sum(vapply(ans, function(z) isTRUE(z$success), logical(1))),
              "/", B)
    }
  } else {
    ans <- run_serial()
  }

  for (b in seq_len(B)) {
    if (isTRUE(ans[[b]]$success)) {
      boot_logWR[b, ] <- ans[[b]]$logWR
      success[b] <- TRUE
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
    ),
    parallel = use_parallel,
    ncores = if (use_parallel) ncores else 1L
  )
}

.local_bootstrap <- function(dat1, dat0, grid, h, tau, logWR_hat, weights,
                             B, conf.level, seed, n_total, verbose) {
  if (!is.null(seed)) set.seed(seed)
  G <- nrow(grid)
  boot_logWR <- matrix(NA_real_, nrow = B, ncol = G)
  success <- logical(B)
  n1 <- nrow(dat1)
  n0 <- nrow(dat0)
  every <- max(1L, floor(B / 10L))

  for (b in seq_len(B)) {
    d1 <- dat1[sample.int(n1, n1, replace = TRUE), , drop = FALSE]
    d0 <- dat0[sample.int(n0, n0, replace = TRUE), , drop = FALSE]
    wrb <- try(.local_estimate_grid(d1, d0, grid, h, tau), silent = TRUE)
    if (!inherits(wrb, "try-error") && all(is.finite(wrb)) && all(wrb > 0)) {
      boot_logWR[b, ] <- log(wrb)
      success[b] <- TRUE
    }
    if (isTRUE(verbose) && (b %% every == 0L || b == B)) {
      message("Bootstrap ", b, "/", B, "; successful curves: ", sum(success))
    }
  }

  n_success <- sum(success)
  required <- if (B >= 20L) max(20L, ceiling(0.80 * B)) else max(2L, ceiling(0.80 * B))
  if (n_success < 2L) stop("Fewer than two bootstrap curves were successful.", call. = FALSE)
  if (n_success < required) {
    warning("Only ", n_success, " of ", B, " bootstrap curves were successful; test p-values are set to NA.", call. = FALSE)
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

  obs_const <- .local_lack_of_fit(logWR_hat, grid, weights, "constant")
  obs_linear <- .local_lack_of_fit(logWR_hat, grid, weights, "linear")
  scale_test <- n_total * prod(h)
  T_const <- rep(NA_real_, B)
  T_linear <- rep(NA_real_, B)
  const_coef <- rep(NA_real_, B)
  linear_coef <- matrix(NA_real_, nrow = B, ncol = ncol(grid) + 1L)

  for (b in which(success)) {
    log_b <- boot_logWR[b, ]
    diff_b <- log_b - logWR_hat
    tc <- .local_lack_of_fit(diff_b, grid, weights, "constant")
    tl <- .local_lack_of_fit(diff_b, grid, weights, "linear")
    if (!is.null(tc)) T_const[b] <- scale_test * tc$D
    if (!is.null(tl)) T_linear[b] <- scale_test * tl$D
    pc <- .local_project(log_b, grid, weights, "constant")
    pl <- .local_project(log_b, grid, weights, "linear")
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

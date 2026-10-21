.ph_bernstein_basis <- function(t, m, lower = 0, upper) {
  if (!is.numeric(m) || length(m) != 1L || m < 0 || m != as.integer(m)) {
    stop("Bernstein degree must be a nonnegative integer.", call. = FALSE)
  }
  if (!is.finite(lower) || !is.finite(upper) || upper <= lower) {
    stop("Invalid Bernstein interval.", call. = FALSE)
  }
  t <- as.numeric(t)
  s <- (t - lower) / (upper - lower)
  s <- pmin(pmax(s, 0), 1)
  B <- vapply(0:as.integer(m), function(k) {
    choose(m, k) * s^k * (1 - s)^(m - k)
  }, numeric(length(s)))
  if (m == 0L) B <- matrix(B, ncol = 1L)
  B
}

.ph_phi_from_eta <- function(eta_inc) {
  c(0, cumsum(exp(as.numeric(eta_inc))))
}

.ph_baseline_eval <- function(t, phi, lower = 0, upper) {
  m <- length(phi) - 1L
  if (m < 1L) stop("Bernstein degree must be at least 1.", call. = FALSE)
  Bm <- .ph_bernstein_basis(t, m = m, lower = lower, upper = upper)
  Bm1 <- .ph_bernstein_basis(t, m = m - 1L, lower = lower, upper = upper)
  dphi <- diff(phi)
  Lambda <- as.vector(Bm %*% phi)
  lambda <- as.vector((m / (upper - lower)) * (Bm1 %*% dphi))
  list(
    Lambda = pmax(Lambda, 0),
    lambda = pmax(lambda, .Machine$double.xmin)
  )
}

.ph_baseline_eval_observed <- function(dat, phi, margin = 1L) {
  m <- length(phi) - 1L
  upper <- attr(dat, "tmax")
  if (margin == 1L) {
    Bm <- attr(dat, "B1_m")
    Bm1 <- attr(dat, "B1_m1")
  } else {
    Bm <- attr(dat, "B2_m")
    Bm1 <- attr(dat, "B2_m1")
  }
  dphi <- diff(phi)
  Lambda <- as.vector(Bm %*% phi)
  lambda <- as.vector((m / upper) * (Bm1 %*% dphi))
  list(
    Lambda = pmax(Lambda, 0),
    lambda = pmax(lambda, .Machine$double.xmin)
  )
}

.ph_bb1_components <- function(u, v, theta, delta) {
  eps <- 1e-12
  if (!is.finite(theta) || !is.finite(delta) || theta <= 0 || delta < 1) {
    stop("Invalid BB1 parameters.", call. = FALSE)
  }
  u <- pmin(pmax(as.numeric(u), eps), 1 - eps)
  v <- pmin(pmax(as.numeric(v), eps), 1 - eps)
  logu <- log(u)
  logv <- log(v)
  A <- pmax(expm1(-theta * logu), 1e-300)
  B <- pmax(expm1(-theta * logv), 1e-300)
  logA <- log(A)
  logB <- log(B)
  x1 <- delta * logA
  x2 <- delta * logB
  mx <- pmax(x1, x2)
  logQ <- mx + log(exp(x1 - mx) + exp(x2 - mx))
  logR <- logQ / delta
  R <- exp(pmin(logR, 700))
  log1pR <- log1p(R)
  logC <- -(1 / theta) * log1pR
  C <- exp(pmax(logC, -745))

  logCu <- (-1 / theta - 1) * log1pR +
    (1 / delta - 1) * logQ +
    (delta - 1) * logA +
    (-theta - 1) * logu
  logCv <- (-1 / theta - 1) * log1pR +
    (1 / delta - 1) * logQ +
    (delta - 1) * logB +
    (-theta - 1) * logv
  Cu <- exp(pmin(pmax(logCu, -745), 0))
  Cv <- exp(pmin(pmax(logCv, -745), 0))

  bracket <- pmax(theta * (delta - 1) * (1 + R) + (1 + theta) * R, 1e-300)
  logDensity <- (-theta - 1) * logu +
    (-theta - 1) * logv +
    (delta - 1) * logA +
    (delta - 1) * logB +
    (1 / delta - 2) * logQ +
    (-1 / theta - 2) * log1pR +
    log(bracket)
  density <- exp(pmin(pmax(logDensity, -745), 700))

  list(
    C = pmin(pmax(C, eps), 1),
    Cu = pmin(pmax(Cu, eps), 1),
    Cv = pmin(pmax(Cv, eps), 1),
    density = pmax(density, eps)
  )
}

.ph_clayton_components <- function(u, v, theta) {
  eps <- 1e-12
  if (!is.finite(theta) || theta <= 0) stop("Clayton theta must be positive.", call. = FALSE)
  u <- pmin(pmax(as.numeric(u), eps), 1 - eps)
  v <- pmin(pmax(as.numeric(v), eps), 1 - eps)
  A <- pmax(u^(-theta) + v^(-theta) - 1, eps)
  C <- A^(-1 / theta)
  Cu <- u^(-theta - 1) * A^(-1 / theta - 1)
  Cv <- v^(-theta - 1) * A^(-1 / theta - 1)
  density <- (1 + theta) * u^(-theta - 1) * v^(-theta - 1) *
    A^(-1 / theta - 2)
  list(
    C = pmax(C, eps),
    Cu = pmin(pmax(Cu, eps), 1),
    Cv = pmin(pmax(Cv, eps), 1),
    density = pmax(density, eps)
  )
}

.ph_gumbel_components <- function(u, v, theta) {
  eps <- 1e-12
  if (!is.finite(theta) || theta <= 1) stop("Gumbel theta must exceed 1.", call. = FALSE)
  u <- pmin(pmax(as.numeric(u), eps), 1 - eps)
  v <- pmin(pmax(as.numeric(v), eps), 1 - eps)
  a <- -log(u)
  b <- -log(v)
  A <- pmax(a^theta + b^theta, eps)
  r <- A^(1 / theta)
  C <- exp(-r)
  common <- A^(1 / theta - 1)
  Cu <- C * common * a^(theta - 1) / u
  Cv <- C * common * b^(theta - 1) / v
  density <- C / (u * v) * a^(theta - 1) * b^(theta - 1) *
    A^(2 / theta - 2) * (1 + (theta - 1) * A^(-1 / theta))
  list(
    C = pmax(C, eps),
    Cu = pmin(pmax(Cu, eps), 1),
    Cv = pmin(pmax(Cv, eps), 1),
    density = pmax(density, eps)
  )
}


.ph_frank_components <- function(u, v, theta) {
  eps <- 1e-12
  if (!is.finite(theta)) stop("Frank theta must be finite.", call. = FALSE)
  u <- pmin(pmax(as.numeric(u), eps), 1 - eps)
  v <- pmin(pmax(as.numeric(v), eps), 1 - eps)

  ## The Frank copula approaches independence as theta -> 0. Treat a tiny
  ## neighborhood of zero by its continuous limit to avoid cancellation.
  if (abs(theta) < 1e-7) {
    return(list(
      C = u * v,
      Cu = v,
      Cv = u,
      density = rep(1, length(u))
    ))
  }

  A <- expm1(-theta * u)
  B <- expm1(-theta * v)
  D <- expm1(-theta)
  H <- D + A * B

  inside <- 1 + (A * B) / D
  inside <- pmax(inside, 1e-300)
  C <- -log(inside) / theta

  eu <- exp(-theta * u)
  ev <- exp(-theta * v)
  Cu <- eu * B / H
  Cv <- ev * A / H
  density <- (-theta * D) * exp(-theta * (u + v)) / (H^2)

  list(
    C = pmin(pmax(C, eps), 1),
    Cu = pmin(pmax(Cu, eps), 1),
    Cv = pmin(pmax(Cv, eps), 1),
    density = pmax(density, eps)
  )
}

.ph_gaussian_components <- function(u, v, rho) {
  eps <- 1e-12
  if (!is.finite(rho) || abs(rho) >= 1) {
    stop("Gaussian rho must lie strictly between -1 and 1.", call. = FALSE)
  }
  u <- pmin(pmax(as.numeric(u), eps), 1 - eps)
  v <- pmin(pmax(as.numeric(v), eps), 1 - eps)
  x <- stats::qnorm(u)
  y <- stats::qnorm(v)
  s2 <- 1 - rho^2
  s <- sqrt(s2)

  C <- pbivnorm::pbivnorm(x, y, rho = rho)
  Cu <- stats::pnorm((y - rho * x) / s)
  Cv <- stats::pnorm((x - rho * y) / s)
  log_density <- -log(s) +
    (2 * rho * x * y - rho^2 * (x^2 + y^2)) / (2 * s2)
  density <- exp(pmin(pmax(log_density, -745), 700))

  list(
    C = pmin(pmax(C, eps), 1),
    Cu = pmin(pmax(Cu, eps), 1),
    Cv = pmin(pmax(Cv, eps), 1),
    density = pmax(density, eps)
  )
}

.ph_frank_debye1 <- function(theta) {
  if (abs(theta) < 1e-6) {
    ## D_1(theta) = 1 - theta/4 + theta^2/36 + O(theta^4)
    return(1 - theta / 4 + theta^2 / 36)
  }
  integrand <- function(t) {
    out <- t / expm1(t)
    out[abs(t) < 1e-8] <- 1
    out
  }
  if (theta > 0) {
    val <- stats::integrate(integrand, 0, theta, rel.tol = 1e-9)$value
  } else {
    val <- -stats::integrate(integrand, theta, 0, rel.tol = 1e-9)$value
  }
  val / theta
}

.ph_copula_components <- function(u, v, dependence, family) {
  family <- match.arg(family, c("bb1", "clayton", "frank", "gaussian", "gumbel"))
  if (family == "bb1") {
    .ph_bb1_components(u, v, dependence$theta, dependence$delta)
  } else if (family == "clayton") {
    .ph_clayton_components(u, v, dependence$theta)
  } else if (family == "frank") {
    .ph_frank_components(u, v, dependence$theta)
  } else if (family == "gaussian") {
    .ph_gaussian_components(u, v, dependence$theta)
  } else {
    .ph_gumbel_components(u, v, dependence$theta)
  }
}

.ph_dependence_summary <- function(family, dependence) {
  if (family == "bb1") {
    theta <- dependence$theta
    delta <- dependence$delta
    return(list(
      kendall_tau = 1 - 2 / (delta * (theta + 2)),
      lower_tail = 2^(-1 / (theta * delta)),
      upper_tail = 2 - 2^(1 / delta)
    ))
  }
  if (family == "clayton") {
    theta <- dependence$theta
    return(list(
      kendall_tau = theta / (theta + 2),
      lower_tail = 2^(-1 / theta),
      upper_tail = 0
    ))
  }
  if (family == "frank") {
    theta <- dependence$theta
    tau <- if (abs(theta) < 1e-7) 0 else {
      1 - 4 / theta + 4 * .ph_frank_debye1(theta) / theta
    }
    return(list(kendall_tau = tau, lower_tail = 0, upper_tail = 0))
  }
  if (family == "gaussian") {
    rho <- dependence$theta
    return(list(
      kendall_tau = 2 * asin(rho) / pi,
      lower_tail = 0,
      upper_tail = 0
    ))
  }
  theta <- dependence$theta
  list(
    kendall_tau = 1 - 1 / theta,
    lower_tail = 0,
    upper_tail = 2 - 2^(1 / theta)
  )
}

.ph_dep_unconstrained <- function(family, dependence) {
  if (family == "bb1") {
    c(log(dependence$theta), log(dependence$delta - 1))
  } else if (family == "clayton") {
    log(dependence$theta)
  } else if (family == "frank") {
    dependence$theta
  } else if (family == "gaussian") {
    atanh(dependence$theta)
  } else {
    log(dependence$theta - 1)
  }
}

.ph_dep_from_unconstrained <- function(family, z) {
  if (family == "bb1") {
    list(theta = exp(z[1L]), delta = 1 + exp(z[2L]))
  } else if (family == "clayton") {
    list(theta = exp(z[1L]))
  } else if (family == "frank") {
    list(theta = z[1L])
  } else if (family == "gaussian") {
    list(theta = tanh(z[1L]))
  } else {
    list(theta = 1 + exp(z[1L]))
  }
}

.ph_unpack <- function(par, p, m1, m2, family) {
  off <- 0L
  beta1 <- par[off + seq_len(p)]
  off <- off + p
  beta2 <- par[off + seq_len(p)]
  off <- off + p
  qdep <- if (family == "bb1") 2L else 1L
  dep_z <- par[off + seq_len(qdep)]
  off <- off + qdep
  eta1 <- par[off + seq_len(m1)]
  off <- off + m1
  eta2 <- par[off + seq_len(m2)]
  list(
    beta1 = beta1,
    beta2 = beta2,
    dependence = .ph_dep_from_unconstrained(family, dep_z),
    phi1 = .ph_phi_from_eta(eta1),
    phi2 = .ph_phi_from_eta(eta2)
  )
}

.ph_dep_bounds <- function(family) {
  if (family == "bb1") {
    return(list(lower = c(log(1e-4), log(1e-4)), upper = c(log(20), log(9))))
  }
  if (family == "clayton") {
    return(list(lower = log(0.005), upper = log(20)))
  }
  if (family == "frank") {
    return(list(lower = -30, upper = 30))
  }
  if (family == "gaussian") {
    return(list(lower = atanh(-0.95), upper = atanh(0.95)))
  }
  list(lower = log(0.001), upper = log(9))
}

.ph_start_dependence <- function(family) {
  if (family == "bb1") {
    return(list(
      list(theta = 0.02, delta = 2.5),
      list(theta = 0.10, delta = 2.5),
      list(theta = 0.50, delta = 2.0),
      list(theta = 1.50, delta = 1.4)
    ))
  }
  if (family == "clayton") {
    return(lapply(c(0.25, 0.75, 1.5, 3), function(z) list(theta = z)))
  }
  if (family == "frank") {
    return(lapply(c(-5, -1, 1, 5), function(z) list(theta = z)))
  }
  if (family == "gaussian") {
    return(lapply(c(-0.6, -0.2, 0.2, 0.6), function(z) list(theta = z)))
  }
  lapply(c(1.15, 1.5, 2, 3), function(z) list(theta = z))
}

.ph_cox_beta_start <- function(dat, margin) {
  X <- as.matrix(attr(dat, "X_mat"))
  p <- ncol(X)
  fit <- try({
    if (margin == 1L) {
      survival::coxph(survival::Surv(dat$U1, dat$Delta1) ~ X, ties = "efron")
    } else {
      survival::coxph(survival::Surv(dat$U2, dat$Delta2) ~ X, ties = "efron")
    }
  }, silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(0, p))
  b <- as.numeric(stats::coef(fit))
  if (length(b) != p || any(!is.finite(b))) rep(0, p) else b
}

.ph_make_start <- function(dat, m1, m2, family, dependence) {
  X <- as.matrix(attr(dat, "X_mat"))
  p <- ncol(X)
  upper <- attr(dat, "tmax")
  lam1 <- max(sum(dat$Delta1) / max(sum(dat$U1), 1e-8), 1e-3)
  lam2 <- max(sum(dat$Delta2) / max(sum(dat$U2), 1e-8), 1e-3)
  inc1 <- rep(lam1 * upper / m1, m1)
  inc2 <- rep(lam2 * upper / m2, m2)
  if (family %in% c("bb1", "frank", "gaussian")) {
    beta1 <- .ph_cox_beta_start(dat, 1L)
    beta2 <- .ph_cox_beta_start(dat, 2L)
  } else {
    beta1 <- rep(0, p)
    beta2 <- rep(0, p)
  }
  c(beta1, beta2, .ph_dep_unconstrained(family, dependence), log(inc1), log(inc2))
}

.ph_individual_loglik <- function(par, dat, m1, m2, family) {
  X <- as.matrix(attr(dat, "X_mat"))
  p <- ncol(X)
  pars <- .ph_unpack(par, p, m1, m2, family)
  lp1 <- as.vector(X %*% pars$beta1)
  lp2 <- as.vector(X %*% pars$beta2)
  if (any(!is.finite(lp1)) || any(!is.finite(lp2)) ||
      any(abs(lp1) > 40) || any(abs(lp2) > 40)) {
    return(rep(-Inf, nrow(dat)))
  }
  b1 <- .ph_baseline_eval_observed(dat, pars$phi1, margin = 1L)
  b2 <- .ph_baseline_eval_observed(dat, pars$phi2, margin = 2L)
  logS1 <- -exp(lp1) * b1$Lambda
  logS2 <- -exp(lp2) * b2$Lambda
  S1 <- exp(pmax(logS1, -700))
  S2 <- exp(pmax(logS2, -700))
  logf1 <- log(b1$lambda) + lp1 + logS1
  logf2 <- log(b2$lambda) + lp2 + logS2
  cc <- try(.ph_copula_components(S1, S2, pars$dependence, family), silent = TRUE)
  if (inherits(cc, "try-error")) return(rep(-Inf, nrow(dat)))

  d1 <- dat$Delta1
  d2 <- dat$Delta2
  both_events <- d1 == 1L & d2 == 1L
  only1 <- d1 == 1L & d2 == 0L
  only2 <- d1 == 0L & d2 == 1L
  both_cens <- d1 == 0L & d2 == 0L
  ll <- numeric(nrow(dat))
  ll[both_events] <- logf1[both_events] + logf2[both_events] + log(cc$density[both_events])
  ll[only1] <- logf1[only1] + log(cc$Cu[only1])
  ll[only2] <- logf2[only2] + log(cc$Cv[only2])
  ll[both_cens] <- log(cc$C[both_cens])
  ll
}

.ph_fit_arm <- function(dat, family, m1, m2, weights = NULL, start = NULL,
                        multistart = TRUE, maxit = 2000L) {
  family <- match.arg(family, c("bb1", "clayton", "frank", "gaussian", "gumbel"))
  n <- nrow(dat)
  p <- ncol(attr(dat, "X_mat"))
  if (is.null(weights)) weights <- rep(1, n)
  if (length(weights) != n || any(!is.finite(weights)) || any(weights <= 0)) {
    stop("Likelihood weights must be positive finite values, one per observation.", call. = FALSE)
  }
  weights <- weights / mean(weights)

  objective <- function(par) {
    ll <- .ph_individual_loglik(par, dat, m1, m2, family)
    if (length(ll) != n || any(!is.finite(ll))) return(1e100)
    val <- sum(weights * ll)
    if (!is.finite(val)) 1e100 else -val
  }

  if (!is.null(start)) {
    starts <- list(as.numeric(start))
  } else {
    deps <- .ph_start_dependence(family)
    starts <- lapply(deps, function(dep) .ph_make_start(dat, m1, m2, family, dep))
    if (!isTRUE(multistart)) starts <- starts[1L]
  }

  dep_bounds <- .ph_dep_bounds(family)
  qdep <- if (family == "bb1") 2L else 1L
  lower <- c(rep(-20, 2 * p), dep_bounds$lower, rep(-20, m1 + m2))
  upper <- c(rep(20, 2 * p), dep_bounds$upper, rep(5, m1 + m2))
  if (length(lower) != 2 * p + qdep + m1 + m2) stop("Internal parameter-bound error.", call. = FALSE)

  fits <- vector("list", length(starts))
  for (k in seq_along(starts)) {
    z <- try(stats::optim(
      par = starts[[k]], fn = objective, method = "L-BFGS-B",
      lower = lower, upper = upper,
      control = list(maxit = as.integer(maxit), factr = 1e7, pgtol = 1e-6)
    ), silent = TRUE)
    if (!inherits(z, "try-error") && is.finite(z$value)) fits[[k]] <- z
  }
  good <- which(vapply(fits, function(z) !is.null(z) && is.finite(z$value), logical(1)))
  if (length(good) == 0L) stop("PH--copula optimization failed for family '", family, "'.", call. = FALSE)
  converged <- good[vapply(fits[good], function(z) z$convergence == 0L, logical(1))]
  candidate <- if (length(converged) > 0L) converged else good
  values <- vapply(fits[candidate], function(z) z$value, numeric(1))
  best <- fits[[candidate[which.min(values)]]]
  pars <- .ph_unpack(best$par, p, m1, m2, family)
  dep_sum <- .ph_dependence_summary(family, pars$dependence)
  npar <- 2 * p + qdep + m1 + m2

  structure(list(
    par_unconstr = best$par,
    convergence = best$convergence,
    message = best$message,
    negloglik = best$value,
    loglik = -best$value,
    AIC = 2 * npar + 2 * best$value,
    BIC = log(n) * npar + 2 * best$value,
    npar = npar,
    nobs = n,
    beta1 = pars$beta1,
    beta2 = pars$beta2,
    dependence = pars$dependence,
    theta = pars$dependence$theta,
    rho = if (family == "gaussian") pars$dependence$theta else NULL,
    delta = if (family == "bb1") pars$dependence$delta else NULL,
    kendall_tau = dep_sum$kendall_tau,
    lower_tail = dep_sum$lower_tail,
    upper_tail = dep_sum$upper_tail,
    phi1 = pars$phi1,
    phi2 = pars$phi2,
    p = p,
    m1 = m1,
    m2 = m2,
    lower = 0,
    upper = attr(dat, "tmax"),
    family = family
  ), class = "phcopula_arm_fit")
}

.ph_margin_from_fit <- function(fit, margin, t, x) {
  x <- as.numeric(x)
  if (length(x) != fit$p) stop("x has wrong dimension.", call. = FALSE)
  if (margin == 1L) {
    beta <- fit$beta1
    phi <- fit$phi1
  } else {
    beta <- fit$beta2
    phi <- fit$phi2
  }
  b <- .ph_baseline_eval(t, phi, lower = fit$lower, upper = fit$upper)
  lp <- sum(x * beta)
  S <- exp(-exp(lp) * b$Lambda)
  f <- b$lambda * exp(lp) * S
  list(S = S, f = f)
}

.ph_wr_from_fits <- function(fit1, fit0, x, tau, rel.tol = 1e-7, subdivisions = 300L) {
  if (!identical(fit1$family, fit0$family)) stop("Treatment arms use different copula families.", call. = FALSE)
  family <- fit1$family
  x <- as.numeric(x)

  W_primary <- stats::integrate(function(t) {
    .ph_margin_from_fit(fit1, 1L, t, x)$S * .ph_margin_from_fit(fit0, 1L, t, x)$f
  }, 0, tau, rel.tol = rel.tol, subdivisions = subdivisions)$value
  L_primary <- stats::integrate(function(t) {
    .ph_margin_from_fit(fit0, 1L, t, x)$S * .ph_margin_from_fit(fit1, 1L, t, x)$f
  }, 0, tau, rel.tol = rel.tol, subdivisions = subdivisions)$value

  S11_tau <- .ph_margin_from_fit(fit1, 1L, tau, x)$S
  S10_tau <- .ph_margin_from_fit(fit0, 1L, tau, x)$S

  W_secondary <- stats::integrate(function(t) {
    m21 <- .ph_margin_from_fit(fit1, 2L, t, x)
    m20 <- .ph_margin_from_fit(fit0, 2L, t, x)
    c1 <- .ph_copula_components(rep(S11_tau, length(t)), m21$S, fit1$dependence, family)
    c0 <- .ph_copula_components(rep(S10_tau, length(t)), m20$S, fit0$dependence, family)
    c1$C * c0$Cv * m20$f
  }, 0, tau, rel.tol = rel.tol, subdivisions = subdivisions)$value

  L_secondary <- stats::integrate(function(t) {
    m21 <- .ph_margin_from_fit(fit1, 2L, t, x)
    m20 <- .ph_margin_from_fit(fit0, 2L, t, x)
    c1 <- .ph_copula_components(rep(S11_tau, length(t)), m21$S, fit1$dependence, family)
    c0 <- .ph_copula_components(rep(S10_tau, length(t)), m20$S, fit0$dependence, family)
    c0$C * c1$Cv * m21$f
  }, 0, tau, rel.tol = rel.tol, subdivisions = subdivisions)$value

  W <- W_primary + W_secondary
  L <- L_primary + L_secondary
  WR <- W / L
  list(
    W_primary = W_primary, W_secondary = W_secondary,
    L_primary = L_primary, L_secondary = L_secondary,
    W = W, L = L, WR = WR, logWR = log(WR)
  )
}

.ph_logwr_grid <- function(fit1, fit0, X_grid, tau, rel.tol = 1e-7, subdivisions = 300L) {
  vapply(seq_len(nrow(X_grid)), function(j) {
    .ph_wr_from_fits(fit1, fit0, X_grid[j, ], tau,
                     rel.tol = rel.tol, subdivisions = subdivisions)$logWR
  }, numeric(1))
}

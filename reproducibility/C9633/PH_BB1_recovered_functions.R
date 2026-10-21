baseline_eval <-
function(
    t,
    phi,
    lower = 0,
    upper) {

  m <-
    length(phi) -
    1L

  if (m < 1L) {
    stop(
      "Bernstein degree must be at least 1."
    )
  }

  Bm <- bernstein_basis(
    t,
    m = m,
    lower = lower,
    upper = upper
  )

  Lambda <- as.vector(
    Bm %*% phi
  )

  dphi <-
    diff(phi)

  Bm1 <- bernstein_basis(
    t,
    m = m - 1L,
    lower = lower,
    upper = upper
  )

  lambda <- as.vector(
    (
      m /
        (upper - lower)
    ) *
      (
        Bm1 %*%
          dphi
      )
  )

  list(
    Lambda =
      pmax(
        Lambda,
        0
      ),

    lambda =
      pmax(
        lambda,
        .Machine$double.xmin
      )
  )
}
bb1_components <-
function(
    u,
    v,
    theta,
    delta) {

  eps <- 1e-12

  if (
    !is.finite(theta) ||
    !is.finite(delta) ||
    theta <= 0 ||
    delta < 1
  ) {

    stop(
      "Invalid BB1 parameters."
    )
  }

  u <- pmin(
    pmax(
      as.numeric(u),
      eps
    ),
    1 - eps
  )

  v <- pmin(
    pmax(
      as.numeric(v),
      eps
    ),
    1 - eps
  )

  logu <- log(u)
  logv <- log(v)

  ## u^{-theta} - 1
  ## expm1 is numerically stable when theta is small.
  A <- expm1(
    -theta * logu
  )

  B <- expm1(
    -theta * logv
  )

  A <- pmax(
    A,
    1e-300
  )

  B <- pmax(
    B,
    1e-300
  )

  logA <- log(A)
  logB <- log(B)

  ## log Q, where
  ## Q = A^delta + B^delta

  x1 <- delta * logA
  x2 <- delta * logB

  mx <- pmax(
    x1,
    x2
  )

  logQ <-
    mx +
    log(
      exp(x1 - mx) +
        exp(x2 - mx)
    )

  logR <-
    logQ /
    delta

  R <- exp(
    pmin(
      logR,
      700
    )
  )

  log1pR <-
    log1p(R)

  ##########################################################
  ## C
  ##########################################################

  logC <-
    -(1 / theta) *
    log1pR

  C <- exp(
    pmax(
      logC,
      -745
    )
  )

  ##########################################################
  ## Partial derivative with respect to u
  ##########################################################

  logCu <-
    (-1 / theta - 1) *
    log1pR +
    (1 / delta - 1) *
    logQ +
    (delta - 1) *
    logA +
    (-theta - 1) *
    logu

  Cu <- exp(
    pmin(
      pmax(
        logCu,
        -745
      ),
      0
    )
  )

  ##########################################################
  ## Partial derivative with respect to v
  ##########################################################

  logCv <-
    (-1 / theta - 1) *
    log1pR +
    (1 / delta - 1) *
    logQ +
    (delta - 1) *
    logB +
    (-theta - 1) *
    logv

  Cv <- exp(
    pmin(
      pmax(
        logCv,
        -745
      ),
      0
    )
  )

  ##########################################################
  ## Copula density
  ##########################################################

  bracket <-
    theta *
    (delta - 1) *
    (1 + R) +
    (1 + theta) *
    R

  bracket <- pmax(
    bracket,
    1e-300
  )

  logDensity <-
    (-theta - 1) *
    logu +
    (-theta - 1) *
    logv +
    (delta - 1) *
    logA +
    (delta - 1) *
    logB +
    (1 / delta - 2) *
    logQ +
    (-1 / theta - 2) *
    log1pR +
    log(bracket)

  density <- exp(
    pmin(
      pmax(
        logDensity,
        -745
      ),
      700
    )
  )

  list(
    C =
      pmin(
        pmax(C, eps),
        1
      ),

    Cu =
      pmin(
        pmax(Cu, eps),
        1
      ),

    Cv =
      pmin(
        pmax(Cv, eps),
        1
      ),

    density =
      pmax(
        density,
        eps
      )
  )
}
bb1_kendall_tau <-
function(
    theta,
    delta) {

  1 -
    2 /
    (
      delta *
        (theta + 2)
    )
}
bb1_lower_tail <-
function(
    theta,
    delta) {

  2^(
    -1 /
      (
        theta *
          delta
      )
  )
}
bb1_upper_tail <-
function(
    delta) {

  2 -
    2^(1 / delta)
}
bernstein_basis <-
function(
    t,
    m,
    lower = 0,
    upper) {

  t <- as.numeric(t)

  s <-
    (t - lower) /
    (upper - lower)

  s <- pmin(
    pmax(s, 0),
    1
  )

  B <- vapply(
    0:m,
    function(k) {

      choose(m, k) *
        s^k *
        (1 - s)^(m - k)
    },
    numeric(length(s))
  )

  matrix(
    B,
    nrow = length(s),
    ncol = m + 1L
  )
}
check_linearity <-
function(data,
                            timevar,
                            eventvar,
                            covariate,
                            arm,
                            endpoint_name,
                            covariate_name,
                            spline_df = 3) {

  d <- data[data$treat == arm, ]

  keep <- complete.cases(
    d[, c(timevar, eventvar, covariate)]
  )

  d <- d[keep, ]

  form_linear <- as.formula(
    paste0(
      "Surv(", timevar, ", ", eventvar, ") ~ ",
      covariate
    )
  )

  form_spline <- as.formula(
    paste0(
      "Surv(", timevar, ", ", eventvar, ") ~ ",
      "ns(", covariate, ", df = ", spline_df, ")"
    )
  )

  fit_linear <- coxph(
    form_linear,
    data = d,
    ties = "efron",
    x = TRUE
  )

  fit_spline <- coxph(
    form_spline,
    data = d,
    ties = "efron",
    x = TRUE
  )

  ll_linear <- as.numeric(logLik(fit_linear))
  ll_spline <- as.numeric(logLik(fit_spline))

  df_linear <- attr(logLik(fit_linear), "df")
  df_spline <- attr(logLik(fit_spline), "df")

  LR <- 2 * (ll_spline - ll_linear)
  df_diff <- df_spline - df_linear

  p_value <- pchisq(
    LR,
    df = df_diff,
    lower.tail = FALSE
  )

  cat("\n")
  cat("============================================================\n")
  cat("Endpoint:   ", endpoint_name, "\n")
  cat("Covariate:  ", covariate_name, "\n")
  cat("Treatment:  ", arm, "\n")
  cat("N:          ", nrow(d), "\n")
  cat("Events:     ", sum(d[[eventvar]]), "\n")
  cat("============================================================\n")

  cat("\nLinear Cox model:\n")
  print(summary(fit_linear))

  cat("\nLikelihood-ratio test for nonlinearity:\n")
  cat("LR statistic =", round(LR, 4), "\n")
  cat("df           =", df_diff, "\n")
  cat("p-value      =", signif(p_value, 5), "\n")

  cat("\nAIC:\n")
  cat("Linear model =", round(AIC(fit_linear), 2), "\n")
  cat("Spline model =", round(AIC(fit_spline), 2), "\n")

  invisible(
    list(
      fit_linear = fit_linear,
      fit_spline = fit_spline,
      LR = LR,
      df = df_diff,
      p_value = p_value,
      n = nrow(d),
      events = sum(d[[eventvar]]),
      AIC_linear = AIC(fit_linear),
      AIC_spline = AIC(fit_spline)
    )
  )
}
check_ph <-
function(data, timevar, eventvar, covariate,
                     arm, endpoint_name, covariate_name) {

  d <- data[data$treat == arm, ]

  keep <- complete.cases(
    d[, c(timevar, eventvar, covariate)]
  )

  d <- d[keep, ]

  form <- as.formula(
    paste0(
      "Surv(", timevar, ", ", eventvar, ") ~ ",
      covariate
    )
  )

  fit <- coxph(
    form,
    data = d,
    x = TRUE
  )

  zph <- cox.zph(fit)

  cat("\n")
  cat("============================================================\n")
  cat("Endpoint:   ", endpoint_name, "\n")
  cat("Covariate:  ", covariate_name, "\n")
  cat("Treatment:  ", arm, "\n")
  cat("N:          ", nrow(d), "\n")
  cat("Events:     ", sum(d[[eventvar]]), "\n")
  cat("============================================================\n\n")

  print(summary(fit))

  cat("\nSchoenfeld residual test:\n")
  print(zph)

  invisible(
    list(
      fit = fit,
      zph = zph
    )
  )
}
clayton_components <-
function(u, v, theta) {

  eps <- 1e-12

  u <- pmin(pmax(as.numeric(u), eps), 1 - eps)
  v <- pmin(pmax(as.numeric(v), eps), 1 - eps)

  A <-
    u^(-theta) +
    v^(-theta) -
    1

  A <- pmax(A, eps)

  C <- A^(-1 / theta)

  Cu <-
    u^(-theta - 1) *
    A^(-1 / theta - 1)

  Cv <-
    v^(-theta - 1) *
    A^(-1 / theta - 1)

  density <-
    (1 + theta) *
    u^(-theta - 1) *
    v^(-theta - 1) *
    A^(-1 / theta - 2)

  list(
    C = pmax(C, eps),
    Cu = pmax(Cu, eps),
    Cv = pmax(Cv, eps),
    density = pmax(density, eps)
  )
}
copula_components <-
function(
    u,
    v,
    theta,
    family) {

  if (family == "clayton") {

    clayton_components(
      u,
      v,
      theta
    )

  } else if (family == "gumbel") {

    gumbel_components(
      u,
      v,
      theta
    )

  } else {

    stop("Unknown copula family.")
  }
}
cox_beta_start <-
function(
    dat,
    margin) {

  X <- as.matrix(
    attr(
      dat,
      "X_mat"
    )
  )

  p <- ncol(X)

  if (margin == 1L) {

    fit <- try(
      coxph(
        Surv(
          dat$U1,
          dat$Delta1
        ) ~ X,
        ties = "efron"
      ),
      silent = TRUE
    )

  } else {

    fit <- try(
      coxph(
        Surv(
          dat$U2,
          dat$Delta2
        ) ~ X,
        ties = "efron"
      ),
      silent = TRUE
    )
  }

  if (
    inherits(
      fit,
      "try-error"
    )
  ) {

    return(
      rep(
        0,
        p
      )
    )
  }

  b <-
    as.numeric(
      coef(fit)
    )

  if (
    length(b) != p ||
    any(
      !is.finite(b)
    )
  ) {

    return(
      rep(
        0,
        p
      )
    )
  }

  b
}
eval_ns_spec <-
function(
    z,
    spec) {

  B <- ns(
    z,
    knots =
      spec$knots,
    Boundary.knots =
      spec$Boundary.knots,
    intercept = FALSE
  )

  Bref <- ns(
    spec$reference,
    knots =
      spec$knots,
    Boundary.knots =
      spec$Boundary.knots,
    intercept = FALSE
  )

  sweep(
    as.matrix(B),
    2,
    as.numeric(Bref),
    FUN = "-"
  )
}
fit_bb1_ph_sieve <-
function(
    dat,
    m1 = 6,
    m2 = 6,
    weights = NULL,
    start = NULL,
    multistart = TRUE,
    maxit = 2000) {

  X <- as.matrix(attr(dat, "X_mat"))
  p <- ncol(X)

  objective <- function(par) {

    val <- loglik_bb1_terminal(
      par = par,
      dat = dat,
      m1 = m1,
      m2 = m2,
      weights = weights
    )

    if (!is.finite(val)) {
      return(1e100)
    }

    -val
  }

  lower_bounds <- c(
    rep(-20, 2 * p),
    log(1e-4),       # theta > 0
    log(1e-4),       # delta = 1 + exp(z_delta)
    rep(-20, m1 + m2)
  )

  upper_bounds <- c(
    rep(20, 2 * p),
    log(20),
    log(9),
    rep(5, m1 + m2)
  )

  starts <- list()

  if (!is.null(start)) {

    starts[[1L]] <- start

  } else if (multistart) {

    start_pairs <- list(
      c(theta = 0.02, delta = 2.5),
      c(theta = 0.10, delta = 2.5),
      c(theta = 0.50, delta = 2.0),
      c(theta = 1.50, delta = 1.4)
    )

    for (ss in start_pairs) {

      starts[[length(starts) + 1L]] <- make_bb1_start(
        dat = dat,
        m1 = m1,
        m2 = m2,
        theta0 = as.numeric(ss["theta"]),
        delta0 = as.numeric(ss["delta"])
      )
    }

  } else {

    starts[[1L]] <- make_bb1_start(
      dat = dat,
      m1 = m1,
      m2 = m2,
      theta0 = 0.05,
      delta0 = 2.5
    )
  }

  fits <- vector("list", length(starts))

  for (k in seq_along(starts)) {

    fit_k <- try(
      optim(
        par = starts[[k]],
        fn = objective,
        method = "L-BFGS-B",
        lower = lower_bounds,
        upper = upper_bounds,
        control = list(
          maxit = maxit,
          factr = 1e7,
          pgtol = 1e-6
        )
      ),
      silent = TRUE
    )

    if (!inherits(fit_k, "try-error") &&
        is.finite(fit_k$value)) {

      fits[[k]] <- fit_k
    }
  }

  good <- which(
    vapply(
      fits,
      function(x) {
        !is.null(x) && is.finite(x$value)
      },
      logical(1)
    )
  )

  if (length(good) == 0L) {
    stop("PH--BB1 optimization failed.")
  }

  converged <- good[
    vapply(
      fits[good],
      function(x) x$convergence == 0L,
      logical(1)
    )
  ]

  if (length(converged) > 0L) {
    candidate <- converged
  } else {
    candidate <- good
  }

  values <- vapply(
    fits[candidate],
    function(x) x$value,
    numeric(1)
  )

  best_index <- candidate[which.min(values)]

  best <- fits[[best_index]]

  pars <- unpack_bb1_par(
    best$par,
    p = p,
    m1 = m1,
    m2 = m2
  )

  npar <- 2 * p + 2 + m1 + m2

  structure(
    list(
      par_unconstr = best$par,
      convergence = best$convergence,
      negloglik = best$value,
      loglik = -best$value,
      AIC = 2 * npar + 2 * best$value,
      npar = npar,

      beta1 = pars$beta1,
      beta2 = pars$beta2,

      theta = pars$theta,
      delta = pars$delta,

      kendall_tau = bb1_kendall_tau(
        pars$theta,
        pars$delta
      ),

      lower_tail = bb1_lower_tail(
        pars$theta,
        pars$delta
      ),

      upper_tail = bb1_upper_tail(
        pars$delta
      ),

      phi1 = pars$phi1,
      phi2 = pars$phi2,

      p = p,
      m1 = m1,
      m2 = m2,
      lower = 0,
      upper = attr(dat, "tmax")
    ),
    class = "phbb1_sieve_fit"
  )
}
fit_family_pair <-
function(family) {

  cat(
    "\nFitting",
    family,
    "copula...\n"
  )

  fit0 <- fit_phcopula_sieve(
    dat0,
    family = family,
    m1 = m1,
    m2 = m2
  )

  fit1 <- fit_phcopula_sieve(
    dat1,
    family = family,
    m1 = m1,
    m2 = m2
  )

  list(
    fit0 = fit0,
    fit1 = fit1,
    total_AIC =
      fit0$AIC +
      fit1$AIC,
    total_loglik =
      fit0$loglik +
      fit1$loglik
  )
}
fit_phcopula_sieve <-
function(
    dat,
    family,
    m1 = 6,
    m2 = 6,
    weights = NULL,
    start = NULL,
    multistart = TRUE,
    maxit = 2000) {

  X <- as.matrix(
    attr(dat, "X_mat")
  )

  p <- ncol(X)

  if (family == "clayton") {

    theta_starts <-
      c(0.25, 0.75, 1.5, 3)

  } else {

    theta_starts <-
      c(1.15, 1.5, 2, 3)
  }

  objective <- function(par) {

    val <- loglik_terminal(
      par = par,
      dat = dat,
      family = family,
      m1 = m1,
      m2 = m2,
      weights = weights
    )

    if (!is.finite(val)) {
      return(1e100)
    }

    -val
  }

  ## Parameter bounds
  lower_beta <- rep(-20, 2 * p)
  upper_beta <- rep(20, 2 * p)

  if (family == "clayton") {

    lower_zeta <- log(0.005)
    upper_zeta <- log(20)

  } else {

    lower_zeta <- log(0.001)
    upper_zeta <- log(9)
  }

  lower_eta <-
    rep(-20, m1 + m2)

  upper_eta <-
    rep(5, m1 + m2)

  lower_bounds <- c(
    lower_beta,
    lower_zeta,
    lower_eta
  )

  upper_bounds <- c(
    upper_beta,
    upper_zeta,
    upper_eta
  )

  starts <- list()

  if (!is.null(start)) {

    starts[[1L]] <- start

  } else if (multistart) {

    for (th in theta_starts) {

      starts[[length(starts) + 1L]] <-
        make_start(
          dat,
          m1 = m1,
          m2 = m2,
          family = family,
          theta0 = th
        )
    }

  } else {

    starts[[1L]] <-
      make_start(
        dat,
        m1 = m1,
        m2 = m2,
        family = family,
        theta0 = theta_starts[2]
      )
  }

  fits <- vector(
    "list",
    length(starts)
  )

  for (k in seq_along(starts)) {

    fit_k <- try(
      optim(
        par = starts[[k]],
        fn = objective,
        method = "L-BFGS-B",
        lower = lower_bounds,
        upper = upper_bounds,
        control = list(
          maxit = maxit
        )
      ),
      silent = TRUE
    )

    if (!inherits(
      fit_k,
      "try-error"
    )) {

      fits[[k]] <- fit_k
    }
  }

  good <- which(
    vapply(
      fits,
      function(x) {
        !is.null(x) &&
        is.finite(x$value)
      },
      logical(1)
    )
  )

  if (length(good) == 0L) {
    stop(
      paste(
        "PH--copula optimization failed for",
        family
      )
    )
  }

  converged <- good[
    vapply(
      fits[good],
      function(x) x$convergence == 0L,
      logical(1)
    )
  ]

  if (length(converged) > 0L) {

    candidate <- converged

  } else {

    candidate <- good
  }

  values <- vapply(
    fits[candidate],
    function(x) x$value,
    numeric(1)
  )

  best_index <-
    candidate[
      which.min(values)
    ]

  best <- fits[[best_index]]

  pars <- unpack_par(
    best$par,
    p = p,
    m1 = m1,
    m2 = m2,
    family = family
  )

  npar <-
    2 * p +
    1 +
    m1 +
    m2

  structure(
    list(
      par_unconstr = best$par,
      convergence = best$convergence,
      negloglik = best$value,
      loglik = -best$value,
      AIC =
        2 * npar +
        2 * best$value,
      npar = npar,

      beta1 = pars$beta1,
      beta2 = pars$beta2,

      theta = pars$theta,

      kendall_tau =
        kendall_tau_from_theta(
          pars$theta,
          family
        ),

      phi1 = pars$phi1,
      phi2 = pars$phi2,

      m1 = m1,
      m2 = m2,
      p = p,

      lower = 0,
      upper = attr(dat, "tmax"),

      family = family
    ),
    class = "phcopula_sieve_fit"
  )
}
get_family_summary <-
function(
    object,
    family) {

  if (inherits(
    object,
    "try-error"
  )) {

    return(
      data.frame(
        Family = family,
        Control_theta = NA,
        Control_Kendall_tau = NA,
        Treatment_theta = NA,
        Treatment_Kendall_tau = NA,
        Total_loglik = NA,
        Total_AIC = Inf
      )
    )
  }

  data.frame(
    Family = family,

    Control_theta =
      object$fit0$theta,

    Control_Kendall_tau =
      object$fit0$kendall_tau,

    Treatment_theta =
      object$fit1$theta,

    Treatment_Kendall_tau =
      object$fit1$kendall_tau,

    Total_loglik =
      object$total_loglik,

    Total_AIC =
      object$total_AIC
  )
}
get_p <-
function(obj) {
  as.numeric(obj$zph$table[1, "p"])
}
gumbel_components <-
function(u, v, theta) {

  eps <- 1e-12

  u <- pmin(pmax(as.numeric(u), eps), 1 - eps)
  v <- pmin(pmax(as.numeric(v), eps), 1 - eps)

  a <- -log(u)
  b <- -log(v)

  A <- a^theta + b^theta
  A <- pmax(A, eps)

  r <- A^(1 / theta)

  C <- exp(-r)

  common <-
    A^(1 / theta - 1)

  Cu <-
    C *
    common *
    a^(theta - 1) /
    u

  Cv <-
    C *
    common *
    b^(theta - 1) /
    v

  density <-
    C / (u * v) *
    a^(theta - 1) *
    b^(theta - 1) *
    A^(2 / theta - 2) *
    (
      1 +
      (theta - 1) *
      A^(-1 / theta)
    )

  list(
    C = pmax(C, eps),
    Cu = pmax(Cu, eps),
    Cv = pmax(Cv, eps),
    density = pmax(density, eps)
  )
}
kendall_tau_from_theta <-
function(
    theta,
    family) {

  if (family == "clayton") {

    theta / (theta + 2)

  } else {

    1 - 1 / theta
  }
}
loglik_bb1_terminal <-
function(
    par,
    dat,
    m1 = 6,
    m2 = 6,
    weights = NULL) {

  X <- as.matrix(
    attr(
      dat,
      "X_mat"
    )
  )

  p <- ncol(X)

  upper <-
    attr(
      dat,
      "tmax"
    )

  if (
    is.null(weights)
  ) {

    weights <-
      rep(
        1,
        nrow(dat)
      )
  }

  pars <-
    unpack_bb1_par(
      par,
      p = p,
      m1 = m1,
      m2 = m2
    )

  if (
    pars$theta < 1e-4 ||
    pars$theta > 20 ||
    pars$delta < 1.0001 ||
    pars$delta > 10
  ) {

    return(-1e100)
  }

  lp1 <-
    as.vector(
      X %*%
        pars$beta1
    )

  lp2 <-
    as.vector(
      X %*%
        pars$beta2
    )

  if (
    any(
      !is.finite(lp1)
    ) ||
    any(
      !is.finite(lp2)
    ) ||
    any(
      abs(lp1) > 40
    ) ||
    any(
      abs(lp2) > 40
    )
  ) {

    return(-1e100)
  }

  b1 <- baseline_eval(
    dat$U1,
    pars$phi1,
    lower = 0,
    upper = upper
  )

  b2 <- baseline_eval(
    dat$U2,
    pars$phi2,
    lower = 0,
    upper = upper
  )

  logS1 <-
    -exp(lp1) *
    b1$Lambda

  logS2 <-
    -exp(lp2) *
    b2$Lambda

  S1 <- exp(
    pmax(
      logS1,
      -700
    )
  )

  S2 <- exp(
    pmax(
      logS2,
      -700
    )
  )

  logf1 <-
    log(
      b1$lambda
    ) +
    lp1 +
    logS1

  logf2 <-
    log(
      b2$lambda
    ) +
    lp2 +
    logS2

  cc <- try(
    bb1_components(
      S1,
      S2,
      pars$theta,
      pars$delta
    ),
    silent = TRUE
  )

  if (
    inherits(
      cc,
      "try-error"
    )
  ) {

    return(-1e100)
  }

  d1 <- dat$Delta1
  d2 <- dat$Delta2

  ll <-
    numeric(
      nrow(dat)
    )

  both_events <-
    d1 == 1L &
    d2 == 1L

  only1 <-
    d1 == 1L &
    d2 == 0L

  only2 <-
    d1 == 0L &
    d2 == 1L

  both_cens <-
    d1 == 0L &
    d2 == 0L

  ll[both_events] <-
    logf1[both_events] +
    logf2[both_events] +
    log(
      cc$density[
        both_events
      ]
    )

  ll[only1] <-
    logf1[only1] +
    log(
      cc$Cu[
        only1
      ]
    )

  ll[only2] <-
    logf2[only2] +
    log(
      cc$Cv[
        only2
      ]
    )

  ll[both_cens] <-
    log(
      cc$C[
        both_cens
      ]
    )

  if (
    any(
      !is.finite(ll)
    )
  ) {

    return(-1e100)
  }

  sum(
    weights *
      ll
  )
}
loglik_terminal <-
function(
    par,
    dat,
    family,
    m1 = 6,
    m2 = 6,
    lower = 0,
    upper = NULL,
    weights = NULL) {

  X <- as.matrix(
    attr(dat, "X_mat")
  )

  if (is.null(X)) {
    stop("X_mat attribute is missing.")
  }

  p <- ncol(X)

  if (is.null(upper)) {
    upper <- attr(dat, "tmax")
  }

  if (is.null(weights)) {
    weights <- rep(1, nrow(dat))
  }

  pars <- unpack_par(
    par,
    p = p,
    m1 = m1,
    m2 = m2,
    family = family
  )

  ## Keep optimization away from numerically extreme dependence.
  if (family == "clayton") {
    if (pars$theta < 0.005 ||
        pars$theta > 20) {
      return(-1e100)
    }
  }

  if (family == "gumbel") {
    if (pars$theta < 1.001 ||
        pars$theta > 10) {
      return(-1e100)
    }
  }

  lp1 <- as.vector(
    X %*% pars$beta1
  )

  lp2 <- as.vector(
    X %*% pars$beta2
  )

  if (any(!is.finite(lp1)) ||
      any(!is.finite(lp2)) ||
      any(abs(lp1) > 40) ||
      any(abs(lp2) > 40)) {

    return(-1e100)
  }

  b1 <- baseline_eval(
    dat$U1,
    pars$phi1,
    lower = lower,
    upper = upper
  )

  b2 <- baseline_eval(
    dat$U2,
    pars$phi2,
    lower = lower,
    upper = upper
  )

  logS1 <-
    -exp(lp1) *
    b1$Lambda

  logS2 <-
    -exp(lp2) *
    b2$Lambda

  S1 <- exp(
    pmax(logS1, -700)
  )

  S2 <- exp(
    pmax(logS2, -700)
  )

  logf1 <-
    log(b1$lambda) +
    lp1 +
    logS1

  logf2 <-
    log(b2$lambda) +
    lp2 +
    logS2

  cc <- copula_components(
    S1,
    S2,
    pars$theta,
    family
  )

  d1 <- dat$Delta1
  d2 <- dat$Delta2

  ll <- numeric(nrow(dat))

  both_events <-
    d1 == 1L &
    d2 == 1L

  only1 <-
    d1 == 1L &
    d2 == 0L

  only2 <-
    d1 == 0L &
    d2 == 1L

  both_cens <-
    d1 == 0L &
    d2 == 0L

  ll[both_events] <-
    logf1[both_events] +
    logf2[both_events] +
    log(cc$density[both_events])

  ll[only1] <-
    logf1[only1] +
    log(cc$Cu[only1])

  ll[only2] <-
    logf2[only2] +
    log(cc$Cv[only2])

  ll[both_cens] <-
    log(cc$C[both_cens])

  if (any(!is.finite(ll))) {
    return(-1e100)
  }

  sum(
    weights * ll
  )
}
logwr_grid <-
function(
    fit1,
    fit0,
    X_grid,
    tau) {

  vapply(
    seq_len(
      nrow(X_grid)
    ),
    function(j) {

      wr_from_fits(
        fit1 = fit1,
        fit0 = fit0,
        x = X_grid[j, ],
        tau = tau
      )$logWR
    },
    numeric(1)
  )
}
logwr_grid_bb1 <-
function(
    fit1,
    fit0,
    X_grid,
    tau) {

  vapply(
    seq_len(
      nrow(X_grid)
    ),
    function(j) {

      wr_from_bb1_fits(
        fit1 =
          fit1,

        fit0 =
          fit0,

        x =
          X_grid[
            j,
            ,
            drop = TRUE
          ],

        tau =
          tau
      )$logWR
    },
    numeric(1)
  )
}
make_bb1_start <-
function(
    dat,
    m1,
    m2,
    theta0,
    delta0) {

  X <- as.matrix(
    attr(
      dat,
      "X_mat"
    )
  )

  p <- ncol(X)

  upper <-
    attr(
      dat,
      "tmax"
    )

  beta1_start <-
    cox_beta_start(
      dat,
      margin = 1L
    )

  beta2_start <-
    cox_beta_start(
      dat,
      margin = 2L
    )

  lam1 <-
    sum(
      dat$Delta1
    ) /
    max(
      sum(dat$U1),
      1e-8
    )

  lam2 <-
    sum(
      dat$Delta2
    ) /
    max(
      sum(dat$U2),
      1e-8
    )

  lam1 <- max(
    lam1,
    1e-3
  )

  lam2 <- max(
    lam2,
    1e-3
  )

  inc1 <-
    rep(
      lam1 *
        upper /
        m1,
      m1
    )

  inc2 <-
    rep(
      lam2 *
        upper /
        m2,
      m2
    )

  c(
    beta1_start,
    beta2_start,

    log(theta0),
    log(delta0 - 1),

    log(inc1),
    log(inc2)
  )
}
make_ns_spec <-
function(z, df = 3) {

  B <- ns(
    z,
    df = df,
    intercept = FALSE
  )

  list(
    knots =
      attr(B, "knots"),

    Boundary.knots =
      attr(B, "Boundary.knots"),

    reference =
      median(
        z,
        na.rm = TRUE
      )
  )
}
make_start <-
function(
    dat,
    m1,
    m2,
    family,
    theta0) {

  X <- as.matrix(
    attr(dat, "X_mat")
  )

  p <- ncol(X)

  upper <- attr(dat, "tmax")

  lam1 <-
    sum(dat$Delta1) /
    max(sum(dat$U1), 1e-8)

  lam2 <-
    sum(dat$Delta2) /
    max(sum(dat$U2), 1e-8)

  lam1 <- max(lam1, 1e-3)
  lam2 <- max(lam2, 1e-3)

  inc1 <-
    rep(
      lam1 * upper / m1,
      m1
    )

  inc2 <-
    rep(
      lam2 * upper / m2,
      m2
    )

  c(
    rep(0, p),
    rep(0, p),
    zeta_from_theta(
      theta0,
      family
    ),
    log(inc1),
    log(inc2)
  )
}
margin_from_bb1_fit <-
function(
    fit,
    margin,
    t,
    x) {

  x <-
    as.numeric(x)

  if (
    length(x) !=
    fit$p
  ) {

    stop(
      "x has wrong dimension."
    )
  }

  if (
    margin == 1L
  ) {

    beta <-
      fit$beta1

    phi <-
      fit$phi1

  } else {

    beta <-
      fit$beta2

    phi <-
      fit$phi2
  }

  b <- baseline_eval(
    t,
    phi,
    lower =
      fit$lower,
    upper =
      fit$upper
  )

  lp <-
    sum(
      x *
      beta
    )

  S <-
    exp(
      -exp(lp) *
        b$Lambda
    )

  f <-
    b$lambda *
    exp(lp) *
    S

  list(
    S = S,
    f = f
  )
}
margin_from_fit <-
function(
    fit,
    margin = 1L,
    t,
    x) {

  x <- as.numeric(x)

  if (length(x) != fit$p) {
    stop("x has wrong dimension.")
  }

  if (margin == 1L) {

    beta <- fit$beta1
    phi <- fit$phi1

  } else {

    beta <- fit$beta2
    phi <- fit$phi2
  }

  b <- baseline_eval(
    t,
    phi,
    lower = fit$lower,
    upper = fit$upper
  )

  lp <- sum(
    x * beta
  )

  S <-
    exp(
      -exp(lp) *
      b$Lambda
    )

  f <-
    b$lambda *
    exp(lp) *
    S

  list(
    S = S,
    f = f
  )
}
phi_from_eta <-
function(
    eta_inc) {

  c(
    0,
    cumsum(
      exp(
        as.numeric(
          eta_inc
        )
      )
    )
  )
}
theta_from_zeta <-
function(
    zeta,
    family) {

  if (family == "clayton") {

    exp(zeta)

  } else {

    1 + exp(zeta)
  }
}
unpack_bb1_par <-
function(
    par,
    p,
    m1,
    m2) {

  off <- 0L

  beta1 <-
    par[
      off +
        seq_len(p)
    ]

  off <- off + p

  beta2 <-
    par[
      off +
        seq_len(p)
    ]

  off <- off + p

  z_theta <-
    par[off + 1L]

  off <- off + 1L

  z_delta <-
    par[off + 1L]

  off <- off + 1L

  eta1 <-
    par[
      off +
        seq_len(m1)
    ]

  off <- off + m1

  eta2 <-
    par[
      off +
        seq_len(m2)
    ]

  theta <-
    exp(z_theta)

  delta <-
    1 +
    exp(z_delta)

  list(
    beta1 = beta1,
    beta2 = beta2,

    theta = theta,
    delta = delta,

    eta1 = eta1,
    eta2 = eta2,

    phi1 =
      phi_from_eta(
        eta1
      ),

    phi2 =
      phi_from_eta(
        eta2
      )
  )
}
unpack_par <-
function(
    par,
    p,
    m1,
    m2,
    family) {

  off <- 0L

  beta1 <-
    par[off + seq_len(p)]

  off <- off + p

  beta2 <-
    par[off + seq_len(p)]

  off <- off + p

  zeta <- par[off + 1L]

  off <- off + 1L

  eta1 <-
    par[off + seq_len(m1)]

  off <- off + m1

  eta2 <-
    par[off + seq_len(m2)]

  list(
    beta1 = beta1,
    beta2 = beta2,
    theta = theta_from_zeta(
      zeta,
      family
    ),
    eta1 = eta1,
    eta2 = eta2,
    phi1 = phi_from_eta(eta1),
    phi2 = phi_from_eta(eta2)
  )
}
weighted_bootstrap <-
function(
    dat1,
    dat0,
    fit1,
    fit0,
    X_grid,
    tau,
    B,
    m1,
    m2,
    seed = NULL,
    maxit = 1200) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  G <- nrow(X_grid)

  boot_logWR <-
    matrix(
      NA_real_,
      nrow = B,
      ncol = G
    )

  logWR_original <-
    logwr_grid(
      fit1,
      fit0,
      X_grid,
      tau
    )

  success <- logical(B)

  for (b in seq_len(B)) {

    if (
      b == 1L ||
      b %% 25L == 0L ||
      b == B
    ) {

      cat(
        "Bootstrap completed",
        b,
        "of",
        B,
        "\n"
      )

      flush.console()
    }

    ## Exp(1) multiplier weights
    w1 <- rexp(nrow(dat1))
    w0 <- rexp(nrow(dat0))

    ## Normalization does not alter the arm-specific maximizer.
    w1 <- w1 / mean(w1)
    w0 <- w0 / mean(w0)

    f1b <- try(
      fit_phcopula_sieve(
        dat1,
        family = fit1$family,
        m1 = m1,
        m2 = m2,
        weights = w1,
        start = fit1$par_unconstr,
        multistart = FALSE,
        maxit = maxit
      ),
      silent = TRUE
    )

    f0b <- try(
      fit_phcopula_sieve(
        dat0,
        family = fit0$family,
        m1 = m1,
        m2 = m2,
        weights = w0,
        start = fit0$par_unconstr,
        multistart = FALSE,
        maxit = maxit
      ),
      silent = TRUE
    )

    if (
      inherits(f1b, "try-error") ||
      inherits(f0b, "try-error")
    ) {

      next
    }

    if (
      f1b$convergence != 0L ||
      f0b$convergence != 0L
    ) {

      next
    }

    vals <- try(
      logwr_grid(
        f1b,
        f0b,
        X_grid,
        tau
      ),
      silent = TRUE
    )

    if (
      inherits(vals, "try-error") ||
      any(!is.finite(vals))
    ) {

      next
    }

    boot_logWR[b, ] <- vals

    success[b] <- TRUE
  }

  n_success <- sum(success)

  cat(
    "\nSuccessful bootstrap fits:",
    n_success,
    "of",
    B,
    "\n"
  )

  se <- apply(
    boot_logWR,
    2,
    sd,
    na.rm = TRUE
  )

  ##########################################################
  ## Normal-approximation CI
  ##########################################################

  zcrit <- qnorm(0.975)

  normal_lower <-
    logWR_original -
    zcrit * se

  normal_upper <-
    logWR_original +
    zcrit * se

  ##########################################################
  ## Centered/basic multiplier-bootstrap CI
  ##########################################################

  boot_diff <- sweep(
    boot_logWR,
    2,
    logWR_original,
    FUN = "-"
  )

  basic_lower <-
    rep(
      NA_real_,
      G
    )

  basic_upper <-
    rep(
      NA_real_,
      G
    )

  for (j in seq_len(G)) {

    vals <-
      boot_diff[, j]

    vals <-
      vals[
        is.finite(vals)
      ]

    if (length(vals) < 20L) {
      next
    }

    qq <- quantile(
      vals,
      probs = c(
        0.025,
        0.975
      ),
      type = 7,
      names = FALSE
    )

    basic_lower[j] <-
      logWR_original[j] -
      qq[2]

    basic_upper[j] <-
      logWR_original[j] -
      qq[1]
  }

  list(
    logWR_hat = logWR_original,
    boot_logWR = boot_logWR,
    success = success,
    n_success = n_success,
    se = se,

    normal_lower = normal_lower,
    normal_upper = normal_upper,

    basic_lower = basic_lower,
    basic_upper = basic_upper
  )
}
weighted_bootstrap_bb1 <-
function(
    dat1,
    dat0,
    fit1,
    fit0,
    X_grid,
    tau,
    B,
    m1,
    m2,
    seed = NULL,
    maxit = 1000) {

  if (
    !is.null(seed)
  ) {

    set.seed(seed)
  }

  G <-
    nrow(X_grid)

  original <-
    logwr_grid_bb1(
      fit1,
      fit0,
      X_grid,
      tau
    )

  boot_logWR <-
    matrix(
      NA_real_,
      nrow = B,
      ncol = G
    )

  success <-
    logical(B)

  for (
    b in
    seq_len(B)
  ) {

    ########################################################
    ## Independent Exp(1) multiplier weights by arm
    ########################################################

    w1 <-
      rexp(
        nrow(dat1)
      )

    w0 <-
      rexp(
        nrow(dat0)
      )

    w1 <-
      w1 /
      mean(w1)

    w0 <-
      w0 /
      mean(w0)

    ########################################################
    ## Warm-start treatment fit
    ########################################################

    f1b <- try(
      fit_bb1_ph_sieve(
        dat1,
        m1 = m1,
        m2 = m2,

        weights =
          w1,

        start =
          fit1$par_unconstr,

        multistart =
          FALSE,

        maxit =
          maxit
      ),
      silent = TRUE
    )

    ########################################################
    ## Warm-start control fit
    ########################################################

    f0b <- try(
      fit_bb1_ph_sieve(
        dat0,
        m1 = m1,
        m2 = m2,

        weights =
          w0,

        start =
          fit0$par_unconstr,

        multistart =
          FALSE,

        maxit =
          maxit
      ),
      silent = TRUE
    )

    ########################################################
    ## Compute bootstrap WR curve if both fits converged
    ########################################################

    if (
      !inherits(
        f1b,
        "try-error"
      ) &&
      !inherits(
        f0b,
        "try-error"
      ) &&
      f1b$convergence == 0L &&
      f0b$convergence == 0L
    ) {

      vals <- try(
        logwr_grid_bb1(
          f1b,
          f0b,
          X_grid,
          tau
        ),
        silent = TRUE
      )

      if (
        !inherits(
          vals,
          "try-error"
        ) &&
        all(
          is.finite(vals)
        )
      ) {

        boot_logWR[
          b,
        ] <- vals

        success[b] <-
          TRUE
      }
    }

    cat(
      "Completed bootstrap",
      b,
      "of",
      B,
      "- successful:",
      sum(success),
      "\n"
    )

    flush.console()
  }

  ##########################################################
  ## Bootstrap SE
  ##########################################################

  se <- apply(
    boot_logWR,
    2,
    function(z) {

      z <-
        z[
          is.finite(z)
        ]

      if (
        length(z) <
        2L
      ) {

        NA_real_

      } else {

        sd(z)
      }
    }
  )

  ##########################################################
  ## Normal bootstrap-SE confidence intervals
  ##########################################################

  zcrit <-
    qnorm(0.975)

  normal_lower <-
    original -
    zcrit *
    se

  normal_upper <-
    original +
    zcrit *
    se

  ##########################################################
  ## Basic weighted-bootstrap confidence intervals
  ##########################################################

  boot_diff <- sweep(
    boot_logWR,
    2,
    original,
    FUN = "-"
  )

  basic_lower <-
    rep(
      NA_real_,
      G
    )

  basic_upper <-
    rep(
      NA_real_,
      G
    )

  for (
    j in
    seq_len(G)
  ) {

    vals <-
      boot_diff[
        ,
        j
      ]

    vals <-
      vals[
        is.finite(vals)
      ]

    if (
      length(vals) <
      20L
    ) {

      next
    }

    qq <- quantile(
      vals,
      probs =
        c(
          0.025,
          0.975
        ),
      type = 7,
      names = FALSE
    )

    basic_lower[j] <-
      original[j] -
      qq[2]

    basic_upper[j] <-
      original[j] -
      qq[1]
  }

  list(
    logWR_hat =
      original,

    boot_logWR =
      boot_logWR,

    success =
      success,

    n_success =
      sum(success),

    se =
      se,

    normal_lower =
      normal_lower,

    normal_upper =
      normal_upper,

    basic_lower =
      basic_lower,

    basic_upper =
      basic_upper
  )
}
wr_from_bb1_fits <-
function(
    fit1,
    fit0,
    x,
    tau,
    rel.tol = 1e-7,
    subdivisions = 300L) {

  x <-
    as.numeric(x)

  ##########################################################
  ## Primary comparison
  ##########################################################

  W_primary <-
    integrate(
      function(t) {

        S11 <-
          margin_from_bb1_fit(
            fit1,
            1L,
            t,
            x
          )$S

        f10 <-
          margin_from_bb1_fit(
            fit0,
            1L,
            t,
            x
          )$f

        S11 *
          f10
      },

      lower = 0,
      upper = tau,

      rel.tol =
        rel.tol,

      subdivisions =
        subdivisions
    )$value

  L_primary <-
    integrate(
      function(t) {

        S10 <-
          margin_from_bb1_fit(
            fit0,
            1L,
            t,
            x
          )$S

        f11 <-
          margin_from_bb1_fit(
            fit1,
            1L,
            t,
            x
          )$f

        S10 *
          f11
      },

      lower = 0,
      upper = tau,

      rel.tol =
        rel.tol,

      subdivisions =
        subdivisions
    )$value

  ##########################################################
  ## Secondary comparison
  ##########################################################

  S11_tau <-
    margin_from_bb1_fit(
      fit1,
      1L,
      tau,
      x
    )$S

  S10_tau <-
    margin_from_bb1_fit(
      fit0,
      1L,
      tau,
      x
    )$S

  W_secondary <-
    integrate(
      function(t) {

        m21 <-
          margin_from_bb1_fit(
            fit1,
            2L,
            t,
            x
          )

        m20 <-
          margin_from_bb1_fit(
            fit0,
            2L,
            t,
            x
          )

        c1 <-
          bb1_components(
            rep(
              S11_tau,
              length(t)
            ),
            m21$S,
            fit1$theta,
            fit1$delta
          )

        c0 <-
          bb1_components(
            rep(
              S10_tau,
              length(t)
            ),
            m20$S,
            fit0$theta,
            fit0$delta
          )

        J1 <-
          c1$C

        minus_dJ0 <-
          c0$Cv *
          m20$f

        J1 *
          minus_dJ0
      },

      lower = 0,
      upper = tau,

      rel.tol =
        rel.tol,

      subdivisions =
        subdivisions
    )$value

  L_secondary <-
    integrate(
      function(t) {

        m21 <-
          margin_from_bb1_fit(
            fit1,
            2L,
            t,
            x
          )

        m20 <-
          margin_from_bb1_fit(
            fit0,
            2L,
            t,
            x
          )

        c1 <-
          bb1_components(
            rep(
              S11_tau,
              length(t)
            ),
            m21$S,
            fit1$theta,
            fit1$delta
          )

        c0 <-
          bb1_components(
            rep(
              S10_tau,
              length(t)
            ),
            m20$S,
            fit0$theta,
            fit0$delta
          )

        J0 <-
          c0$C

        minus_dJ1 <-
          c1$Cv *
          m21$f

        J0 *
          minus_dJ1
      },

      lower = 0,
      upper = tau,

      rel.tol =
        rel.tol,

      subdivisions =
        subdivisions
    )$value

  W <-
    W_primary +
    W_secondary

  L <-
    L_primary +
    L_secondary

  WR <-
    W / L

  list(
    W_primary =
      W_primary,

    W_secondary =
      W_secondary,

    L_primary =
      L_primary,

    L_secondary =
      L_secondary,

    W = W,

    L = L,

    WR = WR,

    logWR =
      log(WR)
  )
}
wr_from_fits <-
function(
    fit1,
    fit0,
    x,
    tau,
    rel.tol = 1e-7,
    subdivisions = 300L) {

  x <- as.numeric(x)

  ##########################################################
  ## Primary-event win probability
  ##########################################################

  f_W_primary <- function(t) {

    S11 <-
      margin_from_fit(
        fit1,
        margin = 1L,
        t = t,
        x = x
      )$S

    f10 <-
      margin_from_fit(
        fit0,
        margin = 1L,
        t = t,
        x = x
      )$f

    S11 * f10
  }

  f_L_primary <- function(t) {

    S10 <-
      margin_from_fit(
        fit0,
        margin = 1L,
        t = t,
        x = x
      )$S

    f11 <-
      margin_from_fit(
        fit1,
        margin = 1L,
        t = t,
        x = x
      )$f

    S10 * f11
  }

  W_primary <- integrate(
    f_W_primary,
    lower = 0,
    upper = tau,
    rel.tol = rel.tol,
    subdivisions = subdivisions
  )$value

  L_primary <- integrate(
    f_L_primary,
    lower = 0,
    upper = tau,
    rel.tol = rel.tol,
    subdivisions = subdivisions
  )$value

  ##########################################################
  ## Secondary-event comparison among pairs surviving
  ## the primary event beyond tau
  ##########################################################

  S11_tau <-
    margin_from_fit(
      fit1,
      margin = 1L,
      t = tau,
      x = x
    )$S

  S10_tau <-
    margin_from_fit(
      fit0,
      margin = 1L,
      t = tau,
      x = x
    )$S

  f_W_secondary <- function(t) {

    m21 <-
      margin_from_fit(
        fit1,
        margin = 2L,
        t = t,
        x = x
      )

    m20 <-
      margin_from_fit(
        fit0,
        margin = 2L,
        t = t,
        x = x
      )

    c1 <- copula_components(
      rep(
        S11_tau,
        length(t)
      ),
      m21$S,
      fit1$theta,
      fit1$family
    )

    c0 <- copula_components(
      rep(
        S10_tau,
        length(t)
      ),
      m20$S,
      fit0$theta,
      fit0$family
    )

    J1 <- c1$C

    minus_dJ0 <-
      c0$Cv *
      m20$f

    J1 *
      minus_dJ0
  }

  f_L_secondary <- function(t) {

    m21 <-
      margin_from_fit(
        fit1,
        margin = 2L,
        t = t,
        x = x
      )

    m20 <-
      margin_from_fit(
        fit0,
        margin = 2L,
        t = t,
        x = x
      )

    c1 <- copula_components(
      rep(
        S11_tau,
        length(t)
      ),
      m21$S,
      fit1$theta,
      fit1$family
    )

    c0 <- copula_components(
      rep(
        S10_tau,
        length(t)
      ),
      m20$S,
      fit0$theta,
      fit0$family
    )

    J0 <- c0$C

    minus_dJ1 <-
      c1$Cv *
      m21$f

    J0 *
      minus_dJ1
  }

  W_secondary <- integrate(
    f_W_secondary,
    lower = 0,
    upper = tau,
    rel.tol = rel.tol,
    subdivisions = subdivisions
  )$value

  L_secondary <- integrate(
    f_L_secondary,
    lower = 0,
    upper = tau,
    rel.tol = rel.tol,
    subdivisions = subdivisions
  )$value

  W <-
    W_primary +
    W_secondary

  L <-
    L_primary +
    L_secondary

  WR <- W / L

  list(
    W_primary = W_primary,
    W_secondary = W_secondary,
    L_primary = L_primary,
    L_secondary = L_secondary,
    W = W,
    L = L,
    WR = WR,
    logWR = log(WR)
  )
}
zeta_from_theta <-
function(
    theta,
    family) {

  if (family == "clayton") {

    log(theta)

  } else {

    log(theta - 1)
  }
}

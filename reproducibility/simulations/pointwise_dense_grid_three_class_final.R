############################################################
## DENSE-GRID POINTWISE SIMULATION FOR FINAL FIGURE
##
## Three classes, three estimators:
##   1. Mao--Wang
##   2. Local IPCW
##   3. PH--copula
##
## PURPOSE
## -------
## This run is FOR FIGURE CONSTRUCTION ONLY.  It does not replace
## the already-finalized five-point inferential simulations.
##
## The dense evaluation grid uses the 31 points
##     x = -0.75, -0.70, ..., 0.70, 0.75
## plus x = -0.375 and 0.375, so that ALL FIVE points from the
## finalized inferential simulations are nested exactly.  Thus the
## final dense grid contains 33 points.
##
## IMPORTANT CONSISTENCY FEATURE
## -----------------------------
## The Local IPCW bandwidth is selected on the ORIGINAL five-point
## grid used in the finalized simulations:
##     -0.75, -0.375, 0, 0.375, 0.75.
## The selected bandwidth is then held fixed while the Local IPCW
## curve is evaluated on the 33-point dense grid.  Thus making the
## plotting grid denser does not itself change the bandwidth rule.
##
## There is NO bootstrap in this program.
##
## Final design retained from the frozen pointwise simulations:
##   R = 500
##   n0 = n1 = 250
##   tau = 5
##   tmax = 8
##   common censoring: min{Exp(0.20), 8}
##   m1 = m2 = 6 for the PH--copula sieve fits
##
## Class-specific Monte Carlo seeds are exactly those used in the
## finalized runs:
##   Class 1: 20261002
##   Class 2: 20261003
##   Class 3: 20261010
##
## Outputs:
##   pointwise_dense_grid_three_class_final.rds
##   pointwise_dense_grid_summary.csv
##   pointwise_dense_grid_original5.csv
##   pointwise_dense_grid_figure.pdf
##   pointwise_dense_grid_figure.png
## plus class-specific checkpoint files.
############################################################

suppressPackageStartupMessages({
  library(copula)
  library(WR)
  library(parallel)
})

############################################################
## SETTINGS
############################################################

R_MC <- 500L
N0 <- 250L
N1 <- 250L
TAU <- 5
TMAX <- 8
CENSOR_LAMBDA <- 0.20
M1 <- 6L
M2 <- 6L
THETA_CLAYTON <- 0.8571428571428571

X_SELECT <- seq(-0.75, 0.75, length.out = 5)
X_DENSE <- sort(unique(c(
  seq(-0.75, 0.75, by = 0.05),
  X_SELECT
)))

H_MULT <- c(0.7, 1.0, 1.3)
C_REF <- 1.25
BW_TOL <- 1.10

CLASS_SEEDS <- c(
  Class1 = 20261002L,
  Class2 = 20261003L,
  Class3 = 20261010L
)

AVAILABLE_CORES <- parallel::detectCores(logical = TRUE)
if (is.na(AVAILABLE_CORES) || AVAILABLE_CORES < 1L) AVAILABLE_CORES <- 1L
N_CORES <- max(1L, min(16L, AVAILABLE_CORES - 1L))
BATCH_SIZE <- 50L
RESUME <- TRUE

OUTPUT_PREFIX <- "pointwise_dense_grid"

`%||%` <- function(a, b) if (!is.null(a)) a else b
rtruncnorm_centered <- function(n, a = -1, b = 1) {
  out <- numeric(n)
  k <- 0L
  while (k < n) {
    z <- rnorm(n - k)
    keep <- (z >= a & z <= b)
    m <- sum(keep)
    if (m > 0L) {
      out[(k + 1L):(k + m)] <- z[keep]
      k <- k + m
    }
  }
  out
}
as_X_matrix <- function(x) {
  X <- as.matrix(x)
  if (is.null(dim(X))) X <- matrix(X, ncol = 1L)
  storage.mode(X) <- "double"
  X
}
resample_arm_dat <- function(dat, idx) {
  out <- dat[idx, , drop = FALSE]
  X <- attr(dat, "X_mat")
  if (!is.null(X)) attr(out, "X_mat") <- X[idx, , drop = FALSE]
  attr(out, "tmax") <- attr(dat, "tmax")
  out
}
simulate_class1_gumbel <- function(
    n0 = 250,
    n1 = 250,
    beta_x = 0.30,
    beta_A = log(1.5),
    gamma = 0.50,
    lambda1 = 0.20,
    lambda2 = 2.00,
    gumbel_theta = 2.0,
    censor_lambda = 0.20,
    censor_beta_x = 0.00,
    censor_beta_A = 0.00,
    tmax = 8,
    x_lower = -1,
    x_upper = 1,
    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(n0 > 0, n1 > 0, gumbel_theta > 1, tmax > 0)
  n <- n0 + n1
  x <- rtruncnorm_centered(n, x_lower, x_upper)
  A <- c(rep(0L, n0), rep(1L, n1))
  ord <- sample.int(n)
  x <- x[ord]
  A <- A[ord]
  cop <- copula::gumbelCopula(param = gumbel_theta, dim = 2)
  U <- copula::rCopula(n, cop)
  benefit <- beta_x * x + beta_A * A + gamma * A * x
  rate1 <- lambda1 * exp(-benefit)
  rate2 <- lambda2 * exp(-benefit)
  T1 <- -log(U[, 1L]) / rate1
  T2 <- -log(U[, 2L]) / rate2
  censor_rate <- censor_lambda * exp(censor_beta_x * x + censor_beta_A * A)
  C_random <- rexp(n, rate = censor_rate)
  C <- pmin(C_random, tmax)
  U1 <- pmin(T1, C)
  Delta1 <- as.integer(T1 <= C)
  U2 <- pmin(T2, T1, C)
  Delta2 <- as.integer(T2 <= pmin(T1, C))
  dat <- data.frame(
    id = seq_len(n),
    x = x,
    A = A,
    U1 = U1,
    Delta1 = Delta1,
    U2 = U2,
    Delta2 = Delta2,
    T1_true = T1,
    T2_true = T2,
    C = C
  )
  attr(dat, "tmax") <- tmax
  attr(dat, "truth") <- list(
    beta_x = beta_x,
    beta_A = beta_A,
    gamma = gamma,
    lambda1 = lambda1,
    lambda2 = lambda2,
    gumbel_theta = gumbel_theta,
    censor_lambda = censor_lambda,
    censor_beta_x = censor_beta_x,
    censor_beta_A = censor_beta_A
  )
  split_arm <- function(a) {
    da <- dat[dat$A == a, , drop = FALSE]
    attr(da, "X_mat") <- matrix(da$x, ncol = 1L,
                                 dimnames = list(NULL, "x"))
    attr(da, "tmax") <- tmax
    da
  }
  list(
    dat = dat,
    dat0 = split_arm(0L),
    dat1 = split_arm(1L),
    true_logWR = function(x0) beta_A + gamma * x0,
    settings = attr(dat, "truth")
  )
}
bernstein_basis <- function(t, m, lower = 0, upper) {
  stopifnot(m >= 0, upper > lower)
  t <- as.numeric(t)
  s <- (t - lower) / (upper - lower)
  s <- pmin(pmax(s, 0), 1)
  B <- vapply(0:m, function(k) {
    choose(m, k) * s^k * (1 - s)^(m - k)
  }, numeric(length(s)))
  if (m == 0L) B <- matrix(B, ncol = 1L)
  B
}
phi_from_eta <- function(eta_inc) {
  c(0, cumsum(exp(as.numeric(eta_inc))))
}
baseline_eval <- function(t, phi, lower = 0, upper) {
  m <- length(phi) - 1L
  if (m < 1L) stop("Bernstein degree must be at least 1.")
  Bm <- bernstein_basis(t, m = m, lower = lower, upper = upper)
  Lambda <- as.vector(Bm %*% phi)
  dphi <- diff(phi)
  Bm1 <- bernstein_basis(t, m = m - 1L, lower = lower, upper = upper)
  lambda <- as.vector((m / (upper - lower)) * (Bm1 %*% dphi))
  list(
    Lambda = pmax(Lambda, 0),
    lambda = pmax(lambda, .Machine$double.xmin)
  )
}
gumbel_components <- function(u, v, theta) {
  if (theta <= 1) stop("Gumbel theta must be > 1 for the interior model.")
  eps <- 1e-12
  u <- pmin(pmax(as.numeric(u), eps), 1 - eps)
  v <- pmin(pmax(as.numeric(v), eps), 1 - eps)
  a <- -log(u)
  b <- -log(v)
  A <- a^theta + b^theta
  A <- pmax(A, eps)
  r <- A^(1 / theta)
  C <- exp(-r)
  common <- A^(1 / theta - 1)
  Cu <- C * common * a^(theta - 1) / u
  Cv <- C * common * b^(theta - 1) / v
  dens <- C / (u * v) *
    a^(theta - 1) * b^(theta - 1) *
    A^(2 / theta - 2) *
    (1 + (theta - 1) * A^(-1 / theta))
  list(
    C = pmax(C, eps),
    Cu = pmax(Cu, eps),
    Cv = pmax(Cv, eps),
    density = pmax(dens, eps)
  )
}
theta_from_zeta <- function(zeta) 1 + exp(zeta)
zeta_from_theta <- function(theta) {
  if (theta <= 1) stop("theta must exceed 1.")
  log(theta - 1)
}
unpack_gumbel_par <- function(par, p, m1, m2) {
  off <- 0L
  beta1 <- par[off + seq_len(p)]
  off <- off + p
  beta2 <- par[off + seq_len(p)]
  off <- off + p
  zeta <- par[off + 1L]
  off <- off + 1L
  eta1 <- par[off + seq_len(m1)]
  off <- off + m1
  eta2 <- par[off + seq_len(m2)]
  list(
    beta1 = beta1,
    beta2 = beta2,
    theta = theta_from_zeta(zeta),
    eta1 = eta1,
    eta2 = eta2,
    phi1 = phi_from_eta(eta1),
    phi2 = phi_from_eta(eta2)
  )
}
make_gumbel_start <- function(dat, m1, m2, theta0 = 2) {
  X <- as.matrix(attr(dat, "X_mat"))
  p <- ncol(X)
  upper <- attr(dat, "tmax") %||% max(dat$U1, dat$U2)
  lam1 <- sum(dat$Delta1) / max(sum(dat$U1), 1e-8)
  lam2 <- sum(dat$Delta2) / max(sum(dat$U2), 1e-8)
  lam1 <- max(lam1, 1e-3)
  lam2 <- max(lam2, 1e-3)
  inc1 <- rep(lam1 * upper / m1, m1)
  inc2 <- rep(lam2 * upper / m2, m2)
  c(
    rep(0, p),
    rep(0, p),
    zeta_from_theta(theta0),
    log(inc1),
    log(inc2)
  )
}
loglik_gumbel_terminal <- function(par, dat, m1 = 6, m2 = 6,
                                   lower = 0, upper = NULL) {
  X <- as.matrix(attr(dat, "X_mat"))
  if (is.null(X)) stop("dat must have an X_mat attribute.")
  p <- ncol(X)
  if (is.null(upper)) upper <- attr(dat, "tmax") %||% max(dat$U1, dat$U2)
  pars <- unpack_gumbel_par(par, p = p, m1 = m1, m2 = m2)
  lp1 <- as.vector(X %*% pars$beta1)
  lp2 <- as.vector(X %*% pars$beta2)
  b1 <- baseline_eval(dat$U1, pars$phi1, lower = lower, upper = upper)
  b2 <- baseline_eval(dat$U2, pars$phi2, lower = lower, upper = upper)
  logS1 <- -exp(lp1) * b1$Lambda
  logS2 <- -exp(lp2) * b2$Lambda
  S1 <- exp(pmax(logS1, -745))
  S2 <- exp(pmax(logS2, -745))
  logf1 <- log(b1$lambda) + lp1 + logS1
  logf2 <- log(b2$lambda) + lp2 + logS2
  cc <- gumbel_components(S1, S2, pars$theta)
  d1 <- dat$Delta1
  d2 <- dat$Delta2
  ll <- numeric(nrow(dat))
  both_events <- (d1 == 1L & d2 == 1L)
  only1 <- (d1 == 1L & d2 == 0L)
  only2 <- (d1 == 0L & d2 == 1L)
  both_cens <- (d1 == 0L & d2 == 0L)
  ll[both_events] <- logf1[both_events] + logf2[both_events] +
    log(cc$density[both_events])
  ll[only1] <- logf1[only1] + log(cc$Cu[only1])
  ll[only2] <- logf2[only2] + log(cc$Cv[only2])
  ll[both_cens] <- log(cc$C[both_cens])
  if (any(!is.finite(ll))) return(-1e100)
  sum(ll)
}
fit_gumbel_ph_sieve <- function(dat,
                                m1 = 6,
                                m2 = 6,
                                lower = 0,
                                upper = NULL,
                                start = NULL,
                                theta_starts = c(1.4, 2.0, 3.5),
                                multistart = TRUE,
                                maxit = 1500,
                                reltol = 1e-8) {
  X <- as.matrix(attr(dat, "X_mat"))
  if (is.null(X)) stop("dat must have an X_mat attribute.")
  p <- ncol(X)
  if (is.null(upper)) upper <- attr(dat, "tmax") %||% max(dat$U1, dat$U2)
  objective <- function(par) {
    val <- loglik_gumbel_terminal(par, dat, m1 = m1, m2 = m2,
                                  lower = lower, upper = upper)
    if (!is.finite(val)) return(1e100)
    -val
  }
  starts <- list()
  if (!is.null(start)) {
    starts[[1L]] <- start
  } else if (multistart) {
    for (th in theta_starts) starts[[length(starts) + 1L]] <-
      make_gumbel_start(dat, m1, m2, theta0 = th)
  } else {
    starts[[1L]] <- make_gumbel_start(dat, m1, m2, theta0 = 2)
  }
  best <- NULL
  for (st in starts) {
    fit1 <- try(
      optim(st, objective, method = "BFGS",
            control = list(maxit = maxit, reltol = reltol)),
      silent = TRUE
    )
    if (!inherits(fit1, "try-error") && is.finite(fit1$value)) {
      if (is.null(best) || fit1$value < best$value) best <- fit1
    }
  }
  if (is.null(best) || best$convergence != 0L) {
    st <- if (!is.null(best)) best$par else starts[[1L]]
    fit2 <- try(
      optim(st, objective, method = "Nelder-Mead",
            control = list(maxit = 2L * maxit, reltol = reltol)),
      silent = TRUE
    )
    if (!inherits(fit2, "try-error") && is.finite(fit2$value) &&
        (is.null(best) || fit2$value < best$value)) best <- fit2
  }
  if (is.null(best)) stop("PH--Gumbel sieve optimization failed.")
  pars <- unpack_gumbel_par(best$par, p = p, m1 = m1, m2 = m2)
  structure(list(
    par_unconstr = best$par,
    convergence = best$convergence,
    negloglik = best$value,
    beta1 = pars$beta1,
    beta2 = pars$beta2,
    theta = pars$theta,
    phi1 = pars$phi1,
    phi2 = pars$phi2,
    m1 = m1,
    m2 = m2,
    p = p,
    lower = lower,
    upper = upper
  ), class = "phgumbel_sieve_fit")
}
margin_from_fit <- function(fit, margin = 1L, t, x) {
  x <- as.numeric(x)
  if (length(x) != fit$p) stop("x has wrong dimension.")
  if (margin == 1L) {
    beta <- fit$beta1
    phi <- fit$phi1
  } else {
    beta <- fit$beta2
    phi <- fit$phi2
  }
  b <- baseline_eval(t, phi, lower = fit$lower, upper = fit$upper)
  lp <- sum(x * beta)
  S <- exp(-exp(lp) * b$Lambda)
  f <- b$lambda * exp(lp) * S
  list(S = S, f = f, Lambda0 = b$Lambda, lambda0 = b$lambda)
}
joint_surv_from_fit <- function(fit, s, t, x) {
  S1 <- margin_from_fit(fit, 1L, s, x)$S
  S2 <- margin_from_fit(fit, 2L, t, x)$S
  gumbel_components(S1, S2, fit$theta)$C
}
wr_from_gumbel_fits <- function(fit1, fit0, x, tau,
                                rel.tol = 1e-7,
                                subdivisions = 300L) {
  x <- as.numeric(x)
  f_Wp <- function(t) {
    S11 <- margin_from_fit(fit1, 1L, t, x)$S
    f10 <- margin_from_fit(fit0, 1L, t, x)$f
    S11 * f10
  }
  f_Lp <- function(t) {
    S10 <- margin_from_fit(fit0, 1L, t, x)$S
    f11 <- margin_from_fit(fit1, 1L, t, x)$f
    S10 * f11
  }
  W_prim <- integrate(f_Wp, 0, tau, rel.tol = rel.tol,
                      subdivisions = subdivisions, stop.on.error = TRUE)$value
  L_prim <- integrate(f_Lp, 0, tau, rel.tol = rel.tol,
                      subdivisions = subdivisions, stop.on.error = TRUE)$value
  S11_tau <- margin_from_fit(fit1, 1L, tau, x)$S
  S10_tau <- margin_from_fit(fit0, 1L, tau, x)$S
  f_Ws <- function(t) {
    m21 <- margin_from_fit(fit1, 2L, t, x)
    m20 <- margin_from_fit(fit0, 2L, t, x)
    c1 <- gumbel_components(rep(S11_tau, length(t)), m21$S, fit1$theta)
    c0 <- gumbel_components(rep(S10_tau, length(t)), m20$S, fit0$theta)
    J1 <- c1$C
    minus_dJ0 <- c0$Cv * m20$f
    J1 * minus_dJ0
  }
  f_Ls <- function(t) {
    m21 <- margin_from_fit(fit1, 2L, t, x)
    m20 <- margin_from_fit(fit0, 2L, t, x)
    c1 <- gumbel_components(rep(S11_tau, length(t)), m21$S, fit1$theta)
    c0 <- gumbel_components(rep(S10_tau, length(t)), m20$S, fit0$theta)
    J0 <- c0$C
    minus_dJ1 <- c1$Cv * m21$f
    J0 * minus_dJ1
  }
  W_sec <- integrate(f_Ws, 0, tau, rel.tol = rel.tol,
                     subdivisions = subdivisions, stop.on.error = TRUE)$value
  L_sec <- integrate(f_Ls, 0, tau, rel.tol = rel.tol,
                     subdivisions = subdivisions, stop.on.error = TRUE)$value
  W <- W_prim + W_sec
  L <- L_prim + L_sec
  WR <- W / L
  list(
    W_prim = W_prim,
    L_prim = L_prim,
    W_sec = W_sec,
    L_sec = L_sec,
    W = W,
    L = L,
    WR = WR,
    logWR = log(WR)
  )
}
logwr_grid_from_fits <- function(fit1, fit0, x_grid, tau, ...) {
  vapply(x_grid, function(x0) {
    wr_from_gumbel_fits(fit1, fit0, x = x0, tau = tau, ...)$logWR
  }, numeric(1))
}
class1_true_margin <- function(t, x, arm, margin,
                               beta_x, beta_A, gamma,
                               lambda1, lambda2) {
  benefit <- beta_x * x + beta_A * arm + gamma * arm * x
  rate <- if (margin == 1L) lambda1 else lambda2
  rate <- rate * exp(-benefit)
  S <- exp(-rate * t)
  f <- rate * S
  list(S = S, f = f)
}
class1_true_wr_numeric <- function(x, tau,
                                   beta_x = 0.30,
                                   beta_A = log(1.5),
                                   gamma = 0.50,
                                   lambda1 = 0.20,
                                   lambda2 = 2.00,
                                   gumbel_theta = 2.0,
                                   rel.tol = 1e-10) {
  mp <- function(arm, margin, t) {
    class1_true_margin(t, x, arm, margin,
                       beta_x, beta_A, gamma, lambda1, lambda2)
  }
  Wp <- integrate(function(t) mp(1, 1L, t)$S * mp(0, 1L, t)$f,
                  0, tau, rel.tol = rel.tol)$value
  Lp <- integrate(function(t) mp(0, 1L, t)$S * mp(1, 1L, t)$f,
                  0, tau, rel.tol = rel.tol)$value
  S11_tau <- mp(1, 1L, tau)$S
  S10_tau <- mp(0, 1L, tau)$S
  Ws <- integrate(function(t) {
    m21 <- mp(1, 2L, t)
    m20 <- mp(0, 2L, t)
    c1 <- gumbel_components(rep(S11_tau, length(t)), m21$S, gumbel_theta)
    c0 <- gumbel_components(rep(S10_tau, length(t)), m20$S, gumbel_theta)
    c1$C * c0$Cv * m20$f
  }, 0, tau, rel.tol = rel.tol)$value
  Ls <- integrate(function(t) {
    m21 <- mp(1, 2L, t)
    m20 <- mp(0, 2L, t)
    c1 <- gumbel_components(rep(S11_tau, length(t)), m21$S, gumbel_theta)
    c0 <- gumbel_components(rep(S10_tau, length(t)), m20$S, gumbel_theta)
    c0$C * c1$Cv * m21$f
  }, 0, tau, rel.tol = rel.tol)$value
  WR <- (Wp + Ws) / (Lp + Ls)
  list(
    WR = WR,
    logWR = log(WR),
    analytic_logWR = beta_A + gamma * x,
    discrepancy = log(WR) - (beta_A + gamma * x)
  )
}
suppressPackageStartupMessages(library(WR))



clayton_components <- function(u, v, theta) {
  if (!is.finite(theta) || theta <= 0) stop("Clayton theta must be positive.")

  eps <- 1e-12
  u <- pmin(pmax(as.numeric(u), eps), 1 - eps)
  v <- pmin(pmax(as.numeric(v), eps), 1 - eps)

  A <- u^(-theta) + v^(-theta) - 1
  A <- pmax(A, eps)

  C <- A^(-1 / theta)
  Cu <- u^(-theta - 1) * A^(-1 / theta - 1)
  Cv <- v^(-theta - 1) * A^(-1 / theta - 1)
  density <- (1 + theta) * u^(-theta - 1) * v^(-theta - 1) *
    A^(-1 / theta - 2)

  list(
    C = pmax(C, eps),
    Cu = pmax(Cu, eps),
    Cv = pmax(Cv, eps),
    density = pmax(density, eps)
  )
}

class2_rates <- function(
    x,
    a,
    lambda1 = 0.20,
    lambda2 = 2.00,
    beta1_x = 0.30,
    beta1_A = log(1.5),
    gamma1 = 0.50,
    beta2_x = 0.10,
    beta2_A = log(1.2),
    gamma2 = -0.25) {

  benefit1 <- beta1_x * x + beta1_A * a + gamma1 * a * x
  benefit2 <- beta2_x * x + beta2_A * a + gamma2 * a * x

  c(
    rate1 = lambda1 * exp(-benefit1),
    rate2 = lambda2 * exp(-benefit2)
  )
}

class2_true_wr <- function(x, tau = 5, theta = 0.8571428571428571) {
  r0 <- class2_rates(x, 0)
  r1 <- class2_rates(x, 1)

  r10 <- unname(r0["rate1"])
  r20 <- unname(r0["rate2"])
  r11 <- unname(r1["rate1"])
  r21 <- unname(r1["rate2"])

  S <- function(t, rate) exp(-rate * t)
  f <- function(t, rate) rate * exp(-rate * t)

  W_primary <- integrate(
    function(t) S(t, r11) * f(t, r10),
    lower = 0,
    upper = tau,
    rel.tol = 1e-10
  )$value

  L_primary <- integrate(
    function(t) S(t, r10) * f(t, r11),
    lower = 0,
    upper = tau,
    rel.tol = 1e-10
  )$value

  S11_tau <- S(tau, r11)
  S10_tau <- S(tau, r10)

  W_secondary <- integrate(
    function(t) {
      S21 <- S(t, r21)
      S20 <- S(t, r20)
      c1 <- clayton_components(rep(S11_tau, length(t)), S21, theta)
      c0 <- clayton_components(rep(S10_tau, length(t)), S20, theta)
      c1$C * c0$Cv * f(t, r20)
    },
    lower = 0,
    upper = tau,
    rel.tol = 1e-10
  )$value

  L_secondary <- integrate(
    function(t) {
      S21 <- S(t, r21)
      S20 <- S(t, r20)
      c1 <- clayton_components(rep(S11_tau, length(t)), S21, theta)
      c0 <- clayton_components(rep(S10_tau, length(t)), S20, theta)
      c0$C * c1$Cv * f(t, r21)
    },
    lower = 0,
    upper = tau,
    rel.tol = 1e-10
  )$value

  W <- W_primary + W_secondary
  L <- L_primary + L_secondary

  c(
    W_primary = W_primary,
    W_secondary = W_secondary,
    W = W,
    L_primary = L_primary,
    L_secondary = L_secondary,
    L = L,
    WR = W / L,
    logWR = log(W / L)
  )
}

rtruncnorm_centered <- function(n, a = -1, b = 1) {
  out <- numeric(n)
  k <- 0L

  while (k < n) {
    z <- rnorm(n - k)
    keep <- z >= a & z <= b
    m <- sum(keep)

    if (m > 0L) {
      out[(k + 1L):(k + m)] <- z[keep]
      k <- k + m
    }
  }

  out
}

simulate_class2_clayton <- function(
    n0 = 250,
    n1 = 250,
    theta = 0.8571428571428571,
    censor_lambda = 0.20,
    tmax = 8,
    seed = 20261003) {

  set.seed(seed)

  make_arm <- function(n, a) {
    x <- rtruncnorm_centered(n, -1, 1)

    cop <- copula::claytonCopula(param = theta, dim = 2)
    U <- copula::rCopula(n, cop)

    benefit1 <- 0.30 * x + log(1.5) * a + 0.50 * a * x
    benefit2 <- 0.10 * x + log(1.2) * a - 0.25 * a * x

    rate1 <- 0.20 * exp(-benefit1)
    rate2 <- 2.00 * exp(-benefit2)

    T1 <- -log(U[, 1L]) / rate1
    T2 <- -log(U[, 2L]) / rate2

    C_random <- rexp(n, rate = censor_lambda)
    C <- pmin(C_random, tmax)

    U1 <- pmin(T1, C)
    Delta1 <- as.integer(T1 <= C)

    U2 <- pmin(T2, T1, C)
    Delta2 <- as.integer(T2 <= pmin(T1, C))

    dat <- data.frame(
      x = x,
      A = a,
      U1 = U1,
      Delta1 = Delta1,
      U2 = U2,
      Delta2 = Delta2,
      T1_true = T1,
      T2_true = T2,
      C = C
    )

    attr(dat, "X_mat") <- matrix(x, ncol = 1L, dimnames = list(NULL, "x"))
    attr(dat, "tmax") <- tmax
    dat
  }

  list(
    dat0 = make_arm(n0, 0L),
    dat1 = make_arm(n1, 1L)
  )
}

bernstein_basis <- function(t, m, lower = 0, upper) {
  stopifnot(m >= 0, upper > lower)

  t <- as.numeric(t)
  s <- (t - lower) / (upper - lower)
  s <- pmin(pmax(s, 0), 1)

  B <- vapply(
    0:m,
    function(k) choose(m, k) * s^k * (1 - s)^(m - k),
    numeric(length(s))
  )

  if (m == 0L) B <- matrix(B, ncol = 1L)
  B
}

phi_from_eta <- function(eta_inc) {
  c(0, cumsum(exp(as.numeric(eta_inc))))
}

baseline_eval <- function(t, phi, lower = 0, upper) {
  m <- length(phi) - 1L
  if (m < 1L) stop("Bernstein degree must be at least 1.")

  Bm <- bernstein_basis(t, m = m, lower = lower, upper = upper)
  Lambda <- as.vector(Bm %*% phi)

  dphi <- diff(phi)
  Bm1 <- bernstein_basis(t, m = m - 1L, lower = lower, upper = upper)
  lambda <- as.vector((m / (upper - lower)) * (Bm1 %*% dphi))

  list(
    Lambda = pmax(Lambda, 0),
    lambda = pmax(lambda, .Machine$double.xmin)
  )
}

theta_from_zeta_clayton <- function(zeta) exp(zeta)

zeta_from_theta_clayton <- function(theta) {
  if (!is.finite(theta) || theta <= 0) stop("theta must be positive.")
  log(theta)
}

unpack_clayton_par <- function(par, p, m1, m2) {
  off <- 0L

  beta1 <- par[off + seq_len(p)]
  off <- off + p

  beta2 <- par[off + seq_len(p)]
  off <- off + p

  zeta <- par[off + 1L]
  off <- off + 1L

  eta1 <- par[off + seq_len(m1)]
  off <- off + m1

  eta2 <- par[off + seq_len(m2)]

  list(
    beta1 = beta1,
    beta2 = beta2,
    theta = theta_from_zeta_clayton(zeta),
    eta1 = eta1,
    eta2 = eta2,
    phi1 = phi_from_eta(eta1),
    phi2 = phi_from_eta(eta2)
  )
}

make_clayton_start <- function(dat, m1, m2, theta0 = 0.8571428571428571) {
  X <- as.matrix(attr(dat, "X_mat"))
  p <- ncol(X)
  upper <- attr(dat, "tmax")

  lam1 <- sum(dat$Delta1) / max(sum(dat$U1), 1e-8)
  lam2 <- sum(dat$Delta2) / max(sum(dat$U2), 1e-8)
  lam1 <- max(lam1, 1e-3)
  lam2 <- max(lam2, 1e-3)

  inc1 <- rep(lam1 * upper / m1, m1)
  inc2 <- rep(lam2 * upper / m2, m2)

  c(
    rep(0, p),
    rep(0, p),
    zeta_from_theta_clayton(theta0),
    log(inc1),
    log(inc2)
  )
}

loglik_clayton_terminal <- function(
    par,
    dat,
    m1 = 6,
    m2 = 6,
    lower = 0,
    upper = NULL) {

  X <- as.matrix(attr(dat, "X_mat"))
  if (is.null(X)) stop("dat must have an X_mat attribute.")

  p <- ncol(X)
  if (is.null(upper)) upper <- attr(dat, "tmax")

  pars <- unpack_clayton_par(par, p = p, m1 = m1, m2 = m2)

  lp1 <- as.vector(X %*% pars$beta1)
  lp2 <- as.vector(X %*% pars$beta2)

  b1 <- baseline_eval(dat$U1, pars$phi1, lower = lower, upper = upper)
  b2 <- baseline_eval(dat$U2, pars$phi2, lower = lower, upper = upper)

  logS1 <- -exp(lp1) * b1$Lambda
  logS2 <- -exp(lp2) * b2$Lambda

  S1 <- exp(pmax(logS1, -745))
  S2 <- exp(pmax(logS2, -745))

  logf1 <- log(b1$lambda) + lp1 + logS1
  logf2 <- log(b2$lambda) + lp2 + logS2

  cc <- clayton_components(S1, S2, pars$theta)

  d1 <- dat$Delta1
  d2 <- dat$Delta2

  both_events <- d1 == 1L & d2 == 1L
  only1 <- d1 == 1L & d2 == 0L
  only2 <- d1 == 0L & d2 == 1L
  both_cens <- d1 == 0L & d2 == 0L

  ll <- numeric(nrow(dat))

  ll[both_events] <- logf1[both_events] + logf2[both_events] +
    log(cc$density[both_events])

  ll[only1] <- logf1[only1] + log(cc$Cu[only1])
  ll[only2] <- logf2[only2] + log(cc$Cv[only2])
  ll[both_cens] <- log(cc$C[both_cens])

  if (any(!is.finite(ll))) return(-1e100)
  sum(ll)
}

fit_clayton_ph_sieve <- function(
    dat,
    m1 = 6,
    m2 = 6,
    lower = 0,
    upper = NULL,
    start = NULL,
    theta_starts = c(0.4, 0.8571428571428571, 1.5),
    multistart = TRUE,
    maxit = 1500,
    reltol = 1e-8) {

  X <- as.matrix(attr(dat, "X_mat"))
  if (is.null(X)) stop("dat must have an X_mat attribute.")

  p <- ncol(X)
  if (is.null(upper)) upper <- attr(dat, "tmax")

  objective <- function(par) {
    val <- loglik_clayton_terminal(
      par,
      dat,
      m1 = m1,
      m2 = m2,
      lower = lower,
      upper = upper
    )

    if (!is.finite(val)) return(1e100)
    -val
  }

  starts <- list()

  if (!is.null(start)) {
    starts[[1L]] <- start
  } else if (multistart) {
    for (th in theta_starts) {
      starts[[length(starts) + 1L]] <- make_clayton_start(dat, m1, m2, theta0 = th)
    }
  } else {
    starts[[1L]] <- make_clayton_start(dat, m1, m2)
  }

  best <- NULL

  for (st in starts) {
    z <- try(
      optim(
        st,
        objective,
        method = "BFGS",
        control = list(maxit = maxit, reltol = reltol)
      ),
      silent = TRUE
    )

    if (!inherits(z, "try-error") && is.finite(z$value)) {
      if (is.null(best) || z$value < best$value) best <- z
    }
  }

  if (is.null(best) || best$convergence != 0L) {
    st <- if (is.null(best)) starts[[1L]] else best$par

    z <- try(
      optim(
        st,
        objective,
        method = "Nelder-Mead",
        control = list(maxit = 2L * maxit, reltol = reltol)
      ),
      silent = TRUE
    )

    if (!inherits(z, "try-error") && is.finite(z$value) &&
        (is.null(best) || z$value < best$value)) {
      best <- z
    }
  }

  if (is.null(best)) stop("PH--Clayton sieve optimization failed.")

  pars <- unpack_clayton_par(best$par, p = p, m1 = m1, m2 = m2)

  list(
    par_unconstr = best$par,
    convergence = best$convergence,
    negloglik = best$value,
    beta1 = pars$beta1,
    beta2 = pars$beta2,
    theta = pars$theta,
    phi1 = pars$phi1,
    phi2 = pars$phi2,
    m1 = m1,
    m2 = m2,
    p = p,
    lower = lower,
    upper = upper
  )
}

margin_from_clayton_fit <- function(fit, margin = 1L, t, x) {
  x <- as.numeric(x)
  if (length(x) != fit$p) stop("x has wrong dimension.")

  if (margin == 1L) {
    beta <- fit$beta1
    phi <- fit$phi1
  } else {
    beta <- fit$beta2
    phi <- fit$phi2
  }

  b <- baseline_eval(t, phi, lower = fit$lower, upper = fit$upper)
  lp <- sum(x * beta)
  S <- exp(-exp(lp) * b$Lambda)
  f <- b$lambda * exp(lp) * S

  list(S = S, f = f)
}

wr_from_clayton_fits <- function(
    fit1,
    fit0,
    x,
    tau = 5,
    rel.tol = 1e-7,
    subdivisions = 300L) {

  Wp <- integrate(
    function(t) {
      S11 <- margin_from_clayton_fit(fit1, 1L, t, x)$S
      f10 <- margin_from_clayton_fit(fit0, 1L, t, x)$f
      S11 * f10
    },
    0,
    tau,
    rel.tol = rel.tol,
    subdivisions = subdivisions,
    stop.on.error = TRUE
  )$value

  Lp <- integrate(
    function(t) {
      S10 <- margin_from_clayton_fit(fit0, 1L, t, x)$S
      f11 <- margin_from_clayton_fit(fit1, 1L, t, x)$f
      S10 * f11
    },
    0,
    tau,
    rel.tol = rel.tol,
    subdivisions = subdivisions,
    stop.on.error = TRUE
  )$value

  S11_tau <- margin_from_clayton_fit(fit1, 1L, tau, x)$S
  S10_tau <- margin_from_clayton_fit(fit0, 1L, tau, x)$S

  Ws <- integrate(
    function(t) {
      m21 <- margin_from_clayton_fit(fit1, 2L, t, x)
      m20 <- margin_from_clayton_fit(fit0, 2L, t, x)
      c1 <- clayton_components(rep(S11_tau, length(t)), m21$S, fit1$theta)
      c0 <- clayton_components(rep(S10_tau, length(t)), m20$S, fit0$theta)
      c1$C * c0$Cv * m20$f
    },
    0,
    tau,
    rel.tol = rel.tol,
    subdivisions = subdivisions,
    stop.on.error = TRUE
  )$value

  Ls <- integrate(
    function(t) {
      m21 <- margin_from_clayton_fit(fit1, 2L, t, x)
      m20 <- margin_from_clayton_fit(fit0, 2L, t, x)
      c1 <- clayton_components(rep(S11_tau, length(t)), m21$S, fit1$theta)
      c0 <- clayton_components(rep(S10_tau, length(t)), m20$S, fit0$theta)
      c0$C * c1$Cv * m21$f
    },
    0,
    tau,
    rel.tol = rel.tol,
    subdivisions = subdivisions,
    stop.on.error = TRUE
  )$value

  W <- Wp + Ws
  L <- Lp + Ls

  if (!is.finite(W) || !is.finite(L) || W <= 0 || L <= 0) return(NA_real_)
  log(W / L)
}

suppressPackageStartupMessages(library(WR))



primary_change <- 2.0
secondary_change <- 0.50
eta_primary <- 1.00


class3_primary_rates <- function(x, a) {

  # Early period: time-varying linear effect plus persistent
  # treatment-specific nonlinear effect eta_primary * a * x^2.
  benefit_early <-
    0.30 * x +
    log(1.5) * a +
    0.50 * a * x +
    eta_primary * a * x^2

  # Late period: the same persistent quadratic treatment effect
  # remains present while the linear covariate effect changes.
  benefit_late <-
    -0.60 * x +
    log(1.2) * a -
    0.40 * a * x +
    eta_primary * a * x^2

  cbind(
    early = 0.20 * exp(-benefit_early),
    late = 0.20 * exp(-benefit_late)
  )
}


class3_secondary_rates <- function(x, a) {

  # Early-period covariate hazard coefficients:
  # control -0.10; treatment +0.15
  benefit_early <-
    0.10 * x +
    log(1.2) * a -
    0.25 * a * x

  # Late-period covariate hazard coefficients:
  # control -0.80; treatment -1.20
  benefit_late <-
    0.80 * x +
    log(1.4) * a +
    0.40 * a * x

  cbind(
    early = 2.00 * exp(-benefit_early),
    late = 2.00 * exp(-benefit_late)
  )
}


class3_margin <- function(t, x, a, outcome) {

  t <- as.numeric(t)

  if (outcome == 1L) {
    rr <- class3_primary_rates(x, a)
    cp <- primary_change
  } else if (outcome == 2L) {
    rr <- class3_secondary_rates(x, a)
    cp <- secondary_change
  } else {
    stop("outcome must be 1 or 2.")
  }

  if (length(x) != 1L || length(a) != 1L) {
    stop("class3_margin is intended for scalar x and a.")
  }

  r_early <- unname(rr[1L, "early"])
  r_late <- unname(rr[1L, "late"])

  cumulative_hazard <-
    r_early * pmin(t, cp) +
    r_late * pmax(t - cp, 0)

  hazard <- ifelse(
    t <= cp,
    r_early,
    r_late
  )

  survival <- exp(-cumulative_hazard)
  density <- hazard * survival

  list(
    hazard = hazard,
    cumulative_hazard = cumulative_hazard,
    survival = survival,
    density = density
  )
}


integrate_piecewise <- function(
    f,
    upper,
    cuts = numeric(0),
    rel.tol = 1e-10) {

  cuts <- sort(
    unique(
      cuts[
        is.finite(cuts) &
        cuts > 0 &
        cuts < upper
      ]
    )
  )

  endpoints <- c(0, cuts, upper)
  total <- 0

  for (j in seq_len(length(endpoints) - 1L)) {
    total <- total +
      integrate(
        f,
        lower = endpoints[j],
        upper = endpoints[j + 1L],
        rel.tol = rel.tol,
        subdivisions = 1000L
      )$value
  }

  total
}


class3_true_wr <- function(
    x,
    tau = 5,
    theta = 0.8571428571428571,
    rel.tol = 1e-10) {

  S1 <- function(t, a) {
    class3_margin(
      t = t,
      x = x,
      a = a,
      outcome = 1L
    )$survival
  }

  f1 <- function(t, a) {
    class3_margin(
      t = t,
      x = x,
      a = a,
      outcome = 1L
    )$density
  }

  S2 <- function(t, a) {
    class3_margin(
      t = t,
      x = x,
      a = a,
      outcome = 2L
    )$survival
  }

  f2 <- function(t, a) {
    class3_margin(
      t = t,
      x = x,
      a = a,
      outcome = 2L
    )$density
  }

  W_primary <- integrate_piecewise(
    function(t) {
      S1(t, 1L) * f1(t, 0L)
    },
    upper = tau,
    cuts = primary_change,
    rel.tol = rel.tol
  )

  L_primary <- integrate_piecewise(
    function(t) {
      S1(t, 0L) * f1(t, 1L)
    },
    upper = tau,
    cuts = primary_change,
    rel.tol = rel.tol
  )

  S11_tau <- S1(tau, 1L)
  S10_tau <- S1(tau, 0L)

  W_secondary <- integrate_piecewise(
    function(t) {
      S21 <- S2(t, 1L)
      S20 <- S2(t, 0L)

      c1 <- clayton_components(
        rep(S11_tau, length(t)),
        S21,
        theta
      )

      c0 <- clayton_components(
        rep(S10_tau, length(t)),
        S20,
        theta
      )

      c1$C * c0$Cv * f2(t, 0L)
    },
    upper = tau,
    cuts = secondary_change,
    rel.tol = rel.tol
  )

  L_secondary <- integrate_piecewise(
    function(t) {
      S21 <- S2(t, 1L)
      S20 <- S2(t, 0L)

      c1 <- clayton_components(
        rep(S11_tau, length(t)),
        S21,
        theta
      )

      c0 <- clayton_components(
        rep(S10_tau, length(t)),
        S20,
        theta
      )

      c0$C * c1$Cv * f2(t, 1L)
    },
    upper = tau,
    cuts = secondary_change,
    rel.tol = rel.tol
  )

  W <- W_primary + W_secondary
  L <- L_primary + L_secondary

  c(
    W_primary = W_primary,
    W_secondary = W_secondary,
    W = W,
    L_primary = L_primary,
    L_secondary = L_secondary,
    L = L,
    WR = W / L,
    logWR = log(W / L)
  )
}


rtruncnorm_centered <- function(n, a = -1, b = 1) {
  out <- numeric(n)
  k <- 0L

  while (k < n) {
    z <- rnorm(n - k)
    keep <- z >= a & z <= b
    m <- sum(keep)

    if (m > 0L) {
      out[(k + 1L):(k + m)] <- z[keep]
      k <- k + m
    }
  }

  out
}


inverse_piecewise_exponential <- function(
    U,
    rate_early,
    rate_late,
    change_point) {

  U <- pmin(
    pmax(
      as.numeric(U),
      .Machine$double.xmin
    ),
    1
  )

  rate_early <- as.numeric(rate_early)
  rate_late <- as.numeric(rate_late)

  H <- -log(U)
  H_change <- rate_early * change_point

  ifelse(
    H <= H_change,
    H / rate_early,
    change_point +
      (H - H_change) / rate_late
  )
}


simulate_class3_clayton <- function(
    n0 = 250,
    n1 = 250,
    theta = 0.8571428571428571,
    censor_lambda = 0.20,
    tmax = 8,
    seed = 20261006) {

  set.seed(seed)

  make_arm <- function(n, a) {
    x <- rtruncnorm_centered(n, -1, 1)

    cop <- copula::claytonCopula(param = theta, dim = 2)
    U <- copula::rCopula(n, cop)

    r1 <- class3_primary_rates(x, a)
    r2 <- class3_secondary_rates(x, a)

    T1 <- inverse_piecewise_exponential(
      U[, 1L],
      rate_early = r1[, "early"],
      rate_late = r1[, "late"],
      change_point = primary_change
    )

    T2 <- inverse_piecewise_exponential(
      U[, 2L],
      rate_early = r2[, "early"],
      rate_late = r2[, "late"],
      change_point = secondary_change
    )

    C_random <- rexp(n, rate = censor_lambda)
    C <- pmin(C_random, tmax)

    U1 <- pmin(T1, C)
    Delta1 <- as.integer(T1 <= C)

    U2 <- pmin(T2, T1, C)
    Delta2 <- as.integer(T2 <= pmin(T1, C))

    dat <- data.frame(
      x = x,
      A = a,
      U1 = U1,
      Delta1 = Delta1,
      U2 = U2,
      Delta2 = Delta2,
      T1_true = T1,
      T2_true = T2,
      C = C
    )

    attr(dat, "X_mat") <- matrix(
      x,
      ncol = 1L,
      dimnames = list(NULL, "x")
    )

    attr(dat, "tmax") <- tmax
    dat
  }

  list(
    dat0 = make_arm(n0, 0L),
    dat1 = make_arm(n1, 1L)
  )
}




kernel_weights_1d <- function(x, x0, h, eps = 1e-12) {
  u <- (x - x0) / h
  w <- dnorm(u)
  sw <- sum(w)
  if (!is.finite(sw) || sw < eps) return(rep(1 / length(x), length(x)))
  w / sw
}

rot_bandwidth_1d <- function(x) {
  sd(x) * length(x)^(-1 / 5)
}

weighted_km_surv <- function(time, status, weights, times, eps = 1e-8) {
  ord <- order(time)
  time <- time[ord]
  status <- status[ord]
  weights <- weights[ord]
  event_times <- sort(unique(time[status == 1]))
  if (length(event_times) == 0L) return(rep(1, length(times)))
  surv_vals <- numeric(length(event_times))
  S <- 1
  for (j in seq_along(event_times)) {
    tj <- event_times[j]
    risk <- sum(weights[time >= tj])
    event <- sum(weights[time == tj & status == 1])
    if (risk > eps) S <- S * (1 - event / risk)
    surv_vals[j] <- max(S, eps)
  }
  out <- vapply(times, function(tt) {
    idx <- which(event_times <= tt)
    if (length(idx) == 0L) 1 else surv_vals[max(idx)]
  }, numeric(1))
  pmax(out, eps)
}

Ghat_common <- function(t_vec, dat, w) {
  weighted_km_surv(
    time = dat$U1,
    status = 1 - dat$Delta1,
    weights = w,
    times = t_vec
  )
}

S1_local_raw <- function(t_vec, dat, w) {
  if (length(t_vec) == 0L) return(numeric(0))
  G <- Ghat_common(t_vec, dat, w)
  vapply(seq_along(t_vec), function(k) {
    sum(w * (dat$U1 > t_vec[k])) / G[k]
  }, numeric(1))
}

J_local_raw <- function(tau, u_vec, dat, w) {
  if (length(u_vec) == 0L) return(numeric(0))
  Gtau <- as.numeric(Ghat_common(tau, dat, w))
  vapply(u_vec, function(u) {
    sum(w * (dat$U1 > tau) * (dat$U2 > u)) / Gtau
  }, numeric(1))
}

wr_local_exact_1d <- function(x0, dat1, dat0, h, tau = 5) {
  w1 <- kernel_weights_1d(dat1$x, x0, h)
  w0 <- kernel_weights_1d(dat0$x, x0, h)
  t0 <- sort(unique(dat0$U1[dat0$Delta1 == 1L & dat0$U1 <= tau]))
  if (length(t0) > 0L) {
    G0 <- Ghat_common(t0, dat0, w0)
    dF0 <- vapply(t0, function(tt) sum(w0[dat0$Delta1 == 1L & dat0$U1 == tt]), numeric(1)) / G0
    Wp <- sum(S1_local_raw(t0, dat1, w1) * dF0)
  } else Wp <- 0
  t1 <- sort(unique(dat1$U1[dat1$Delta1 == 1L & dat1$U1 <= tau]))
  if (length(t1) > 0L) {
    G1 <- Ghat_common(t1, dat1, w1)
    dF1 <- vapply(t1, function(tt) sum(w1[dat1$Delta1 == 1L & dat1$U1 == tt]), numeric(1)) / G1
    Lp <- sum(S1_local_raw(t1, dat0, w0) * dF1)
  } else Lp <- 0
  G1tau <- as.numeric(Ghat_common(tau, dat1, w1))
  G0tau <- as.numeric(Ghat_common(tau, dat0, w0))
  u0 <- sort(unique(dat0$U2[dat0$U1 > tau & dat0$Delta2 == 1L & dat0$U2 <= tau]))
  if (length(u0) > 0L) {
    dJ0 <- vapply(u0, function(uu) sum(w0[dat0$U1 > tau & dat0$Delta2 == 1L & dat0$U2 == uu]), numeric(1)) / G0tau
    Ws <- sum(J_local_raw(tau, u0, dat1, w1) * dJ0)
  } else Ws <- 0
  u1 <- sort(unique(dat1$U2[dat1$U1 > tau & dat1$Delta2 == 1L & dat1$U2 <= tau]))
  if (length(u1) > 0L) {
    dJ1 <- vapply(u1, function(uu) sum(w1[dat1$U1 > tau & dat1$Delta2 == 1L & dat1$U2 == uu]), numeric(1)) / G1tau
    Ls <- sum(J_local_raw(tau, u1, dat0, w0) * dJ1)
  } else Ls <- 0
  W <- Wp + Ws
  L <- Lp + Ls
  if (!is.finite(W) || !is.finite(L) || W <= 0 || L <= 0) return(NA_real_)
  W / L
}

estimate_local_grid_1d <- function(dat1, dat0, x_grid, h, tau = 5) {
  vapply(x_grid, function(x0) wr_local_exact_1d(x0, dat1, dat0, h, tau), numeric(1))
}

select_bandwidth_local_1d <- function(dat1, dat0, x_grid, tau = 5, h_mult = c(0.7, 1, 1.3), c_ref = 1.25, tol = 1.10) {
  h0 <- rot_bandwidth_1d(c(dat1$x, dat0$x))
  h_vec <- h_mult * h0
  scores <- rep(Inf, length(h_vec))
  for (j in seq_along(h_vec)) {
    a <- estimate_local_grid_1d(dat1, dat0, x_grid, h_vec[j], tau)
    b <- estimate_local_grid_1d(dat1, dat0, x_grid, c_ref * h_vec[j], tau)
    good <- is.finite(a) & is.finite(b) & a > 0 & b > 0
    if (any(good)) scores[j] <- mean((log(a[good]) - log(b[good]))^2)
  }
  if (!any(is.finite(scores))) return(list(h_opt = h0, scores = scores, h_vec = h_vec))
  idx <- which(scores <= tol * min(scores, na.rm = TRUE))[1L]
  list(h_opt = h_vec[idx], scores = scores, h_vec = h_vec)
}

bootstrap_local_grid_1d <- function(dat1, dat0, x_grid, h, tau = 5, B = 200, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  wr_hat <- estimate_local_grid_1d(dat1, dat0, x_grid, h, tau)
  log_hat <- ifelse(is.finite(wr_hat) & wr_hat > 0, log(wr_hat), NA_real_)
  boot <- matrix(NA_real_, B, length(x_grid))
  for (b in seq_len(B)) {
    i1 <- sample.int(nrow(dat1), nrow(dat1), replace = TRUE)
    i0 <- sample.int(nrow(dat0), nrow(dat0), replace = TRUE)
    d1 <- dat1[i1, , drop = FALSE]
    d0 <- dat0[i0, , drop = FALSE]
    wrb <- try(estimate_local_grid_1d(d1, d0, x_grid, h, tau), silent = TRUE)
    if (!inherits(wrb, "try-error")) boot[b, ] <- ifelse(is.finite(wrb) & wrb > 0, log(wrb), NA_real_)
  }
  se <- apply(boot, 2, sd, na.rm = TRUE)
  list(logwr_hat = log_hat, se = se, boot_logwr = boot, n_success = colSums(is.finite(boot)))
}



make_long_for_pwreg <- function(T1, T2, C) {
  n <- length(T1)
  ID <- integer(0)
  time <- numeric(0)
  status <- integer(0)
  for (i in seq_len(n)) {
    first <- min(T1[i], T2[i], C[i])
    if (first == T2[i]) {
      ID <- c(ID, i)
      time <- c(time, T2[i])
      status <- c(status, 2L)
      second <- min(T1[i], C[i])
      ID <- c(ID, i)
      time <- c(time, second)
      status <- c(status, ifelse(second == T1[i], 1L, 0L))
    } else if (first == T1[i]) {
      ID <- c(ID, i)
      time <- c(time, T1[i])
      status <- c(status, 1L)
    } else {
      ID <- c(ID, i)
      time <- c(time, C[i])
      status <- c(status, 0L)
    }
  }
  data.frame(ID = ID, time = time, status = status)
}

make_maowang_data <- function(sim) {
  dat <- if (!is.null(sim$dat)) sim$dat else rbind(sim$dat0, sim$dat1)
  long <- make_long_for_pwreg(dat$T1_true, dat$T2_true, dat$C)
  Zsub <- cbind(x = dat$x, A = dat$A, Ax = dat$A * dat$x)
  list(ID = long$ID, time = long$time, status = long$status, Z = Zsub[long$ID, , drop = FALSE])
}

fit_maowang_grid_1d <- function(sim, x_grid) {
  d <- make_maowang_data(sim)
  fit <- try(WR::pwreg(ID = d$ID, time = d$time, status = d$status, Z = d$Z, rho = 0, eps = 1e-4, maxiter = 50), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  bhat <- fit$beta
  V <- fit$Var
  if (length(bhat) < 3L || any(dim(V) < c(3L, 3L))) return(NULL)
  theta <- c(bhat[2L], bhat[3L])
  Vtheta <- V[c(2L, 3L), c(2L, 3L), drop = FALSE]
  Xg <- cbind(1, x_grid)
  est <- as.vector(Xg %*% theta)
  se <- sqrt(pmax(diag(Xg %*% Vtheta %*% t(Xg)), 0))
  list(logwr_hat = est, se = se)
}

resample_arm_dat <- function(dat, idx) {
  out <- dat[idx, , drop = FALSE]
  X <- attr(dat, "X_mat")
  if (!is.null(X)) attr(out, "X_mat") <- X[idx, , drop = FALSE]
  attr(out, "tmax") <- attr(dat, "tmax")
  out
}

logwr_grid_from_clayton_fits <- function(fit1, fit0, x_grid, tau) {
  vapply(
    x_grid,
    function(xx) wr_from_clayton_fits(fit1, fit0, xx, tau = tau),
    numeric(1)
  )
}




############################################################
## CLASS-SPECIFIC DISPATCHERS
############################################################

simulate_dense_class <- function(class_id, seed, n0, n1) {
  if (class_id == 1L) {
    return(
      simulate_class1_gumbel(
        n0 = n0,
        n1 = n1,
        beta_x = 0.30,
        beta_A = log(1.5),
        gamma = 0.50,
        lambda1 = 0.20,
        lambda2 = 2.00,
        gumbel_theta = 2.0,
        censor_lambda = CENSOR_LAMBDA,
        tmax = TMAX,
        seed = seed
      )
    )
  }

  if (class_id == 2L) {
    return(
      simulate_class2_clayton(
        n0 = n0,
        n1 = n1,
        theta = THETA_CLAYTON,
        censor_lambda = CENSOR_LAMBDA,
        tmax = TMAX,
        seed = seed
      )
    )
  }

  if (class_id == 3L) {
    return(
      simulate_class3_clayton(
        n0 = n0,
        n1 = n1,
        theta = THETA_CLAYTON,
        censor_lambda = CENSOR_LAMBDA,
        tmax = TMAX,
        seed = seed
      )
    )
  }

  stop("Unknown class_id.")
}

truth_dense_class <- function(class_id, x_grid, tau = TAU) {
  if (class_id == 1L) {
    return(log(1.5) + 0.50 * x_grid)
  }

  if (class_id == 2L) {
    return(
      vapply(
        x_grid,
        function(xx) {
          unname(class2_true_wr(xx, tau = tau, theta = THETA_CLAYTON)["logWR"])
        },
        numeric(1)
      )
    )
  }

  if (class_id == 3L) {
    return(
      vapply(
        x_grid,
        function(xx) {
          unname(class3_true_wr(xx, tau = tau, theta = THETA_CLAYTON)["logWR"])
        },
        numeric(1)
      )
    )
  }

  stop("Unknown class_id.")
}

fit_ph_dense_class <- function(class_id, sim, x_grid, tau, m1, m2) {
  if (class_id == 1L) {
    f1 <- fit_gumbel_ph_sieve(sim$dat1, m1 = m1, m2 = m2)
    f0 <- fit_gumbel_ph_sieve(sim$dat0, m1 = m1, m2 = m2)

    if (f1$convergence != 0L || f0$convergence != 0L) {
      return(list(est = rep(NA_real_, length(x_grid)), conv1 = f1$convergence, conv0 = f0$convergence))
    }

    est <- logwr_grid_from_fits(f1, f0, x_grid, tau)

    return(list(est = est, conv1 = f1$convergence, conv0 = f0$convergence))
  }

  f1 <- fit_clayton_ph_sieve(sim$dat1, m1 = m1, m2 = m2)
  f0 <- fit_clayton_ph_sieve(sim$dat0, m1 = m1, m2 = m2)

  if (f1$convergence != 0L || f0$convergence != 0L) {
    return(list(est = rep(NA_real_, length(x_grid)), conv1 = f1$convergence, conv0 = f0$convergence))
  }

  est <- logwr_grid_from_clayton_fits(f1, f0, x_grid, tau)

  list(est = est, conv1 = f1$convergence, conv0 = f0$convergence)
}

############################################################
## ONE DENSE-GRID REPLICATION
############################################################

run_one_dense_rep <- function(r, class_id, settings) {
  G <- length(settings$x_dense)

  out <- list(
    r = r,
    mw = rep(NA_real_, G),
    local = rep(NA_real_, G),
    ph = rep(NA_real_, G),
    h_used = NA_real_,
    ph_conv1 = NA_integer_,
    ph_conv0 = NA_integer_,
    error = NA_character_
  )

  tryCatch({
    sim <- simulate_dense_class(
      class_id = class_id,
      seed = settings$seed + r,
      n0 = settings$n0,
      n1 = settings$n1
    )

    ## Mao--Wang: fit once, evaluate the fitted linear function densely.
    mw <- try(
      fit_maowang_grid_1d(sim, settings$x_dense),
      silent = TRUE
    )

    if (!inherits(mw, "try-error") && !is.null(mw)) {
      out$mw <- mw$logwr_hat
    }

    ## Local IPCW: select h on the original five-point grid, exactly as
    ## in the frozen simulation, then evaluate on the dense plotting grid.
    sel <- try(
      select_bandwidth_local_1d(
        sim$dat1,
        sim$dat0,
        settings$x_select,
        settings$tau,
        settings$h_mult,
        settings$c_ref,
        settings$tol
      ),
      silent = TRUE
    )

    if (!inherits(sel, "try-error") && is.finite(sel$h_opt)) {
      out$h_used <- sel$h_opt

      wr_loc <- try(
        estimate_local_grid_1d(
          sim$dat1,
          sim$dat0,
          settings$x_dense,
          sel$h_opt,
          settings$tau
        ),
        silent = TRUE
      )

      if (!inherits(wr_loc, "try-error")) {
        good <- is.finite(wr_loc) & wr_loc > 0
        out$local[good] <- log(wr_loc[good])
      }
    }

    ## PH--copula: fit once in each arm, then evaluate densely.
    ph <- try(
      fit_ph_dense_class(
        class_id = class_id,
        sim = sim,
        x_grid = settings$x_dense,
        tau = settings$tau,
        m1 = settings$m1,
        m2 = settings$m2
      ),
      silent = TRUE
    )

    if (!inherits(ph, "try-error") && !is.null(ph)) {
      out$ph <- ph$est
      out$ph_conv1 <- ph$conv1
      out$ph_conv0 <- ph$conv0
    }
  }, error = function(e) {
    out$error <<- conditionMessage(e)
  })

  out
}

############################################################
## SUMMARIZE ONE CLASS
############################################################

summarize_dense_class <- function(results, class_id, settings) {
  G <- length(settings$x_dense)

  get_mat <- function(name) {
    do.call(
      rbind,
      lapply(results, function(z) {
        if (is.null(z)) rep(NA_real_, G) else z[[name]]
      })
    )
  }

  mw <- get_mat("mw")
  local <- get_mat("local")
  ph <- get_mat("ph")

  truth <- truth_dense_class(class_id, settings$x_dense, settings$tau)

  mean_safe <- function(M) {
    vapply(seq_len(ncol(M)), function(j) mean(M[, j], na.rm = TRUE), numeric(1))
  }

  sd_safe <- function(M) {
    vapply(seq_len(ncol(M)), function(j) sd(M[, j], na.rm = TRUE), numeric(1))
  }

  n_safe <- function(M) {
    colSums(is.finite(M))
  }

  mw_mean <- mean_safe(mw)
  loc_mean <- mean_safe(local)
  ph_mean <- mean_safe(ph)

  mw_sd <- sd_safe(mw)
  loc_sd <- sd_safe(local)
  ph_sd <- sd_safe(ph)

  n_mw <- n_safe(mw)
  n_local <- n_safe(local)
  n_ph <- n_safe(ph)

  summary <- data.frame(
    class = class_id,
    x = settings$x_dense,
    truth = truth,
    mean_mw = mw_mean,
    mean_local = loc_mean,
    mean_ph = ph_mean,
    bias_mw = mw_mean - truth,
    bias_local = loc_mean - truth,
    bias_ph = ph_mean - truth,
    sd_mw = mw_sd,
    sd_local = loc_sd,
    sd_ph = ph_sd,
    mcse_mean_mw = mw_sd / sqrt(pmax(n_mw, 1)),
    mcse_mean_local = loc_sd / sqrt(pmax(n_local, 1)),
    mcse_mean_ph = ph_sd / sqrt(pmax(n_ph, 1)),
    n_mw = n_mw,
    n_local = n_local,
    n_ph = n_ph
  )

  h_used <- vapply(results, function(z) if (is.null(z)) NA_real_ else z$h_used, numeric(1))
  conv1 <- vapply(results, function(z) if (is.null(z)) NA_integer_ else z$ph_conv1, integer(1))
  conv0 <- vapply(results, function(z) if (is.null(z)) NA_integer_ else z$ph_conv0, integer(1))
  errors <- vapply(results, function(z) {
    if (is.null(z) || is.na(z$error)) "" else z$error
  }, character(1))

  diagnostics <- list(
    mean_local_bandwidth = mean(h_used, na.rm = TRUE),
    local_bandwidths = h_used,
    ph_convergence_arm1 = table(conv1, useNA = "ifany"),
    ph_convergence_arm0 = table(conv0, useNA = "ifany"),
    error_count = sum(nzchar(errors)),
    errors = errors[nzchar(errors)]
  )

  list(
    summary = summary,
    diagnostics = diagnostics,
    raw = list(mw = mw, local = local, ph = ph)
  )
}

############################################################
## PARALLEL CLASS RUNNER WITH CHECKPOINTING
############################################################

run_dense_class_parallel <- function(class_id, settings, n_cores, batch_size, resume = TRUE) {
  class_label <- paste0("Class", class_id)
  checkpoint_file <- paste0(OUTPUT_PREFIX, "_", class_label, "_checkpoint.rds")

  results <- vector("list", settings$R)

  if (resume && file.exists(checkpoint_file)) {
    old <- try(readRDS(checkpoint_file), silent = TRUE)

    if (!inherits(old, "try-error") &&
        is.list(old) &&
        identical(old$class_id, class_id) &&
        identical(old$settings, settings) &&
        length(old$results) == settings$R) {
      results <- old$results
      cat("Resuming", class_label, "from checkpoint.\n")
    }
  }

  pending <- which(vapply(results, is.null, logical(1)))

  if (length(pending) > 0L) {
    cat("\n============================================================\n")
    cat(class_label, "dense-grid point estimation\n")
    cat("R =", settings$R, " | grid points =", length(settings$x_dense), " | cores =", n_cores, "\n")
    cat("============================================================\n")

    cl <- parallel::makeCluster(n_cores, type = "PSOCK")
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages({
        library(copula)
        library(WR)
      })
      NULL
    })

    export_names <- c(
      "%||%",
      "rtruncnorm_centered", "as_X_matrix", "resample_arm_dat",
      "simulate_class1_gumbel",
      "bernstein_basis", "phi_from_eta", "baseline_eval",
      "gumbel_components", "theta_from_zeta", "zeta_from_theta",
      "unpack_gumbel_par", "make_gumbel_start", "loglik_gumbel_terminal",
      "fit_gumbel_ph_sieve", "margin_from_fit", "joint_surv_from_fit",
      "wr_from_gumbel_fits", "logwr_grid_from_fits",
      "class1_true_margin", "class1_true_wr_numeric",
      "clayton_components", "class2_rates", "class2_true_wr",
      "simulate_class2_clayton",
      "theta_from_zeta_clayton", "zeta_from_theta_clayton",
      "unpack_clayton_par", "make_clayton_start", "loglik_clayton_terminal",
      "fit_clayton_ph_sieve", "margin_from_clayton_fit", "wr_from_clayton_fits",
      "primary_change", "secondary_change", "eta_primary",
      "class3_primary_rates", "class3_secondary_rates", "class3_margin",
      "integrate_piecewise", "class3_true_wr", "inverse_piecewise_exponential",
      "simulate_class3_clayton",
      "kernel_weights_1d", "rot_bandwidth_1d", "weighted_km_surv",
      "Ghat_common", "S1_local_raw", "J_local_raw", "wr_local_exact_1d",
      "estimate_local_grid_1d", "select_bandwidth_local_1d",
      "make_long_for_pwreg", "make_maowang_data", "fit_maowang_grid_1d",
      "logwr_grid_from_clayton_fits",
      "simulate_dense_class", "fit_ph_dense_class", "run_one_dense_rep",
      "CENSOR_LAMBDA", "TMAX", "THETA_CLAYTON"
    )

    parallel::clusterExport(cl, export_names, envir = .GlobalEnv)

    batches <- split(pending, ceiling(seq_along(pending) / batch_size))

    for (idx in batches) {
      ans <- parallel::parLapply(
        cl,
        idx,
        run_one_dense_rep,
        class_id = class_id,
        settings = settings
      )

      for (k in seq_along(idx)) {
        results[[idx[k]]] <- ans[[k]]
      }

      saveRDS(
        list(class_id = class_id, settings = settings, results = results),
        checkpoint_file
      )

      completed <- sum(!vapply(results, is.null, logical(1)))
      cat("Completed", completed, "/", settings$R, "for", class_label, "\n")
      flush.console()
    }

    parallel::stopCluster(cl)
    on.exit(NULL, add = FALSE)
  }

  summarize_dense_class(results, class_id, settings)
}

############################################################
## RUN ALL THREE CLASSES
############################################################

time0 <- Sys.time()

cat("\n============================================================\n")
cat("DENSE-GRID THREE-CLASS POINTWISE SIMULATION\n")
cat("============================================================\n")
cat("Dense grid:", length(X_DENSE), "points from", min(X_DENSE), "to", max(X_DENSE), "\n")
cat("Bandwidth-selection grid:", paste(X_SELECT, collapse = ", "), "\n")
cat("No bootstrap.\n")
cat("R =", R_MC, " | n/arm =", N0, " | tau =", TAU, " | cores =", N_CORES, "\n")
cat("============================================================\n")

all_results <- vector("list", 3L)
names(all_results) <- c("Class1", "Class2", "Class3")

for (class_id in 1:3) {
  seed_name <- paste0("Class", class_id)

  settings <- list(
    R = R_MC,
    n0 = N0,
    n1 = N1,
    tau = TAU,
    m1 = M1,
    m2 = M2,
    seed = unname(CLASS_SEEDS[seed_name]),
    x_dense = X_DENSE,
    x_select = X_SELECT,
    h_mult = H_MULT,
    c_ref = C_REF,
    tol = BW_TOL
  )

  all_results[[class_id]] <- run_dense_class_parallel(
    class_id = class_id,
    settings = settings,
    n_cores = N_CORES,
    batch_size = BATCH_SIZE,
    resume = RESUME
  )

  saveRDS(all_results, paste0(OUTPUT_PREFIX, "_partial.rds"))
}

############################################################
## COMBINE SUMMARIES
############################################################

summary_all <- do.call(
  rbind,
  lapply(all_results, function(z) z$summary)
)
rownames(summary_all) <- NULL

class_names <- c(
  "PWF correct; PH--copula correct",
  "PWF fails; PH--copula correct",
  "PWF fails; PH margins misspecified"
)
summary_all$class_label <- class_names[summary_all$class]

## Keep the original five grid points for an internal check against
## the already-finalized five-point simulations.
idx5 <- vapply(
  X_SELECT,
  function(xx) which.min(abs(X_DENSE - xx)),
  integer(1)
)

original5 <- do.call(
  rbind,
  lapply(1:3, function(k) {
    z <- all_results[[k]]$summary[idx5, , drop = FALSE]
    z$class_label <- class_names[k]
    z
  })
)
rownames(original5) <- NULL

cat("\n\n============================================================\n")
cat("ORIGINAL FIVE-POINT GRID CHECK\n")
cat("These point estimates should reproduce the means from the frozen runs.\n")
cat("============================================================\n")
print(
  original5[, c("class", "x", "truth", "mean_mw", "mean_local", "mean_ph",
                 "n_mw", "n_local", "n_ph")],
  row.names = FALSE
)

cat("\n\n============================================================\n")
cat("DIAGNOSTICS\n")
cat("============================================================\n")
for (k in 1:3) {
  cat("\nClass", k, "\n")
  print(all_results[[k]]$diagnostics)
}

############################################################
## SAVE NUMERIC RESULTS
############################################################

final_object <- list(
  settings = list(
    R = R_MC,
    n0 = N0,
    n1 = N1,
    tau = TAU,
    tmax = TMAX,
    censor_lambda = CENSOR_LAMBDA,
    m1 = M1,
    m2 = M2,
    theta_clayton = THETA_CLAYTON,
    x_dense = X_DENSE,
    x_select = X_SELECT,
    h_mult = H_MULT,
    c_ref = C_REF,
    bw_tol = BW_TOL,
    class_seeds = CLASS_SEEDS,
    n_cores = N_CORES
  ),
  summary = summary_all,
  original5 = original5,
  classes = all_results
)

saveRDS(
  final_object,
  "pointwise_dense_grid_three_class_final.rds"
)

write.csv(
  summary_all,
  "pointwise_dense_grid_summary.csv",
  row.names = FALSE
)

write.csv(
  original5,
  "pointwise_dense_grid_original5.csv",
  row.names = FALSE
)

############################################################
## THREE-PANEL FIGURE
############################################################

plot_dense_figure <- function(filename, device = c("pdf", "png")) {
  device <- match.arg(device)

  if (device == "pdf") {
    pdf(filename, width = 10.5, height = 3.6, useDingbats = FALSE)
  } else {
    png(filename, width = 3150, height = 1080, res = 300)
  }

  oldpar <- par(no.readonly = TRUE)
  on.exit({
    par(oldpar)
    dev.off()
  }, add = TRUE)

  par(
    mfrow = c(1, 3),
    mar = c(4.2, 4.2, 2.6, 1.0),
    oma = c(0.5, 0.8, 0.2, 0.2),
    mgp = c(2.4, 0.8, 0),
    las = 1
  )

  panel_titles <- c("Class 1", "Class 2", "Class 3")

  for (k in 1:3) {
    d <- summary_all[summary_all$class == k, , drop = FALSE]

    ylim <- range(
      c(d$truth, d$mean_mw, d$mean_local, d$mean_ph),
      finite = TRUE
    )
    pad <- 0.08 * diff(ylim)
    if (!is.finite(pad) || pad <= 0) pad <- 0.05
    ylim <- ylim + c(-pad, pad)

    plot(
      d$x,
      d$truth,
      type = "l",
      lwd = 2.4,
      lty = 1,
      xlab = expression(x),
      ylab = if (k == 1L) expression(log~WR(x)) else "",
      ylim = ylim,
      main = panel_titles[k]
    )

    lines(d$x, d$mean_mw, lwd = 1.8, lty = 2)
    lines(d$x, d$mean_local, lwd = 1.8, lty = 3)
    lines(d$x, d$mean_ph, lwd = 1.8, lty = 4)

    if (k == 1L) {
      legend(
        "topleft",
        legend = c("Truth", "Mao--Wang", "Local IPCW", "PH--copula"),
        lty = c(1, 2, 3, 4),
        lwd = c(2.4, 1.8, 1.8, 1.8),
        bty = "n",
        cex = 0.78
      )
    }
  }
}

plot_dense_figure(
  "pointwise_dense_grid_figure.pdf",
  device = "pdf"
)

plot_dense_figure(
  "pointwise_dense_grid_figure.png",
  device = "png"
)

cat("\nSaved:\n")
cat("  pointwise_dense_grid_three_class_final.rds\n")
cat("  pointwise_dense_grid_summary.csv\n")
cat("  pointwise_dense_grid_original5.csv\n")
cat("  pointwise_dense_grid_figure.pdf\n")
cat("  pointwise_dense_grid_figure.png\n")

cat("\nTotal elapsed time:\n")
print(Sys.time() - time0)

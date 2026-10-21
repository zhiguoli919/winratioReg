############################################################
## PH--COPULA-ONLY 2D LOG-LINEARITY TEST SIMULATION
##
## PURPOSE
## -------
## Add PH--copula rejection rates to the existing 2D
## log-linearity test simulation WITHOUT rerunning Local IPCW.
##
## This script matches the finalized Local-IPCW 2D testing design:
##   - H0, H11 (mild), H12 (strong) linearity DGPs
##   - n0 = n1 = 300, 500, 800 in the final run
##   - 5 x 5 regular evaluation grid over [-1,1]^2
##   - Corr(X1,X2) = 0.30
##   - tau = 5, administrative censoring at tmax = 8
##   - censoring baseline rate = 0.05
##   - R = 500 Monte Carlo replications, B = 200 bootstraps
##   - same data-generating seed construction as the Local-IPCW code
##
## PH--copula model
## ----------------
## The event-time DGP uses a Clayton CDF copula followed by qexp(U, rate),
## so the corresponding survival copula is the 180-degree rotated Clayton.
## The likelihood below is written on the Clayton CDF scale.
##
## A COMMON 2D cubic PH basis is used under H0, H11, and H12:
##   z(x1,x2) = (x1, x2, x1^2, x2^2, x1*x2, x1^3, x2^3).
## This keeps the fitted PH working model fixed across hypotheses and contains
## the true marginal PH predictor in every 2D linearity-test DGP. Under H0,
## the unnecessary nonlinear coefficients are simply zero.
##
## Bootstrap test
## --------------
## Arm-specific Exp(1) weighted multiplier bootstrap. Each weighted PH--copula
## likelihood is refitted from the original MLE as a warm start, using the same
## Bernstein sieve degrees as the original fit. The test is based on the
## Q-weighted residual from the best linear projection of log WR on
## (1, x1, x2):
##
##   T_M,n = N {M lhat}' Q {M lhat},
##
## with centered multiplier statistic
##
##   T_M,n* = N {M(lhat* - lhat)}' Q {M(lhat* - lhat)}.
##
## Multiplying both observed and bootstrap statistics by the same N does not
## affect the bootstrap p-value or rejection decision.
##
## RUN MODES
## ---------
##   pilot : n=300 only; R=4;   B=10
##   check : n=300 only; R=32;  B=30
##   final : n=300,500,800; R=500; B=200
############################################################

time0 <- Sys.time()

suppressPackageStartupMessages({
  library(copula)
  library(MASS)
  library(parallel)
})

############################################################
## USER SETTINGS
############################################################

RUN_MODE <- "pilot"  # "pilot", "check", or "final"
if (!RUN_MODE %in% c("pilot", "check", "final")) {
  stop("RUN_MODE must be 'pilot', 'check', or 'final'.")
}

TAU <- 5
TMAX <- 8
ALPHA <- 0.05
CENSOR_RATE <- 0.05
RHO_X <- 0.30

GRID_2D_VALUES <- seq(-1, 1, length.out = 5)
GRID_2D <- as.matrix(
  expand.grid(
    x1 = GRID_2D_VALUES,
    x2 = GRID_2D_VALUES
  )
)

## Bernstein sieve degrees.
M1 <- 6L
M2 <- 6L

## Simpson grid used only to evaluate the smooth fitted PH--copula WR surface.
N_INT <- 201L
if (N_INT %% 2L != 1L) stop("N_INT must be odd for Simpson integration.")

## Higher-resolution grid for the one-time population-surface diagnostic.
N_INT_POP <- 2001L
if (N_INT_POP %% 2L != 1L) stop("N_INT_POP must be odd for Simpson integration.")

## Optimization settings. Original fits use multistart; multiplier fits use
## the original MLE as a warm start.
ORIG_MAXIT <- 1200L
BOOT_MAXIT <- 1200L
OPT_RELTOL <- 1e-7

N_AVAILABLE <- parallel::detectCores(logical = TRUE)
if (!is.finite(N_AVAILABLE) || N_AVAILABLE < 2L) N_AVAILABLE <- 2L

if (RUN_MODE == "pilot") {
  SAMPLE_SIZES <- 300L
  R_MC <- 4L
  B_BOOT <- 10L
  N_CORES <- max(1L, min(4L, N_AVAILABLE - 1L))
  BATCH_SIZE <- 1L
  PROGRESS_EVERY <- 1L
  RESUME <- FALSE
} else if (RUN_MODE == "check") {
  SAMPLE_SIZES <- 300L
  R_MC <- 32L
  B_BOOT <- 30L
  N_CORES <- max(1L, min(16L, N_AVAILABLE - 1L))
  BATCH_SIZE <- 1L
  PROGRESS_EVERY <- 4L
  RESUME <- FALSE
} else {
  SAMPLE_SIZES <- c(300L, 500L, 800L)
  R_MC <- 500L
  B_BOOT <- 200L
  N_CORES <- max(1L, min(16L, N_AVAILABLE - 1L))
  BATCH_SIZE <- 10L
  PROGRESS_EVERY <- 50L
  RESUME <- TRUE
}

## Matches the finalized 2D Local-IPCW simulation seed root.
BASE_SEED <- 20261101L
DESIGN_VERSION <- paste0(
  "phcopula_2d_loglinearity_claytonCDF_commonCubic_5x5_",
  RUN_MODE, "_v1"
)
OUTPUT_DIR <- paste0("phcopula_2d_loglinearity_", RUN_MODE, "_output")
CHECKPOINT_DIR <- file.path(OUTPUT_DIR, "checkpoints")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CHECKPOINT_DIR, showWarnings = FALSE, recursive = TRUE)

############################################################
## TERMINAL-EVENT OBSERVATION SCHEME
############################################################

make_observed_terminal_2d <- function(x1, x2, T1, T2, C, A) {
  U1 <- pmin(T1, C)
  Delta1 <- as.integer(T1 <= C)
  U2 <- pmin(T2, T1, C)
  Delta2 <- as.integer(T2 <= pmin(T1, C))

  data.frame(
    x1 = x1,
    x2 = x2,
    A = A,
    T1 = T1,
    T2 = T2,
    C = C,
    U1 = U1,
    Delta1 = Delta1,
    U2 = U2,
    Delta2 = Delta2
  )
}

############################################################
## EXACT 2D LINEARITY-TEST DGPs FROM THE LOCAL-IPCW PROGRAM
############################################################

make_linearity_dgp_2d <- function(
    scenario = c("H0", "mild", "strong")) {

  scenario <- match.arg(scenario)

  if (scenario == "H0") {
    arm1 <- list(
      eta_death_fun = function(x1, x2) {
        -0.30 + 0 * x1 + 0 * x2
      },
      eta_prog_fun = function(x1, x2) {
        -0.20 + 0.45 * x1 + 0.35 * x2
      },
      lambda_death = 0.06,
      lambda_prog = 0.22,
      theta = 1.2
    )

    arm0 <- list(
      eta_death_fun = function(x1, x2) {
        0 * x1 + 0 * x2
      },
      eta_prog_fun = function(x1, x2) {
        0.30 * x1 + 0.20 * x2
      },
      lambda_death = 0.06,
      lambda_prog = 0.22,
      theta = 1.2
    )
  }

  if (scenario == "mild") {
    arm1 <- list(
      eta_death_fun = function(x1, x2) {
        -0.30 +
          0.16 * x1^2 -
          0.12 * x2^2 +
          0.14 * x1 * x2 +
          0.08 * x1^3
      },
      eta_prog_fun = function(x1, x2) {
        -0.20 +
          0.45 * x1 +
          0.35 * x2 +
          0.28 * x1^2 -
          0.22 * x2^2 +
          0.25 * x1 * x2 +
          0.12 * x1^3 -
          0.10 * x2^3
      },
      lambda_death = 0.06,
      lambda_prog = 0.22,
      theta = 1.2
    )

    arm0 <- list(
      eta_death_fun = function(x1, x2) {
        0 * x1 + 0 * x2
      },
      eta_prog_fun = function(x1, x2) {
        0.30 * x1 + 0.20 * x2
      },
      lambda_death = 0.06,
      lambda_prog = 0.22,
      theta = 1.2
    )
  }

  if (scenario == "strong") {
    arm1 <- list(
      eta_death_fun = function(x1, x2) {
        -0.30 +
          0.32 * x1^2 -
          0.25 * x2^2 +
          0.30 * x1 * x2 +
          0.18 * x1^3 -
          0.12 * x2^3
      },
      eta_prog_fun = function(x1, x2) {
        -0.20 +
          0.45 * x1 +
          0.35 * x2 +
          0.55 * x1^2 -
          0.45 * x2^2 +
          0.50 * x1 * x2 +
          0.25 * x1^3 -
          0.22 * x2^3
      },
      lambda_death = 0.06,
      lambda_prog = 0.22,
      theta = 1.2
    )

    arm0 <- list(
      eta_death_fun = function(x1, x2) {
        0 * x1 + 0 * x2
      },
      eta_prog_fun = function(x1, x2) {
        0.30 * x1 + 0.20 * x2
      },
      lambda_death = 0.06,
      lambda_prog = 0.22,
      theta = 1.2
    )
  }

  list(arm1 = arm1, arm0 = arm0)
}

############################################################
## SIMULATE ONE ARM: EXACTLY MATCH THE 2D LOCAL-IPCW DGP
############################################################

simulate_arm_2d <- function(
    n,
    arm,
    A,
    seed,
    tmax = TMAX,
    censor_rate = CENSOR_RATE,
    rho_x = RHO_X) {

  set.seed(seed)

  Sigma <- matrix(
    c(1, rho_x, rho_x, 1),
    2,
    2
  )

  X <- MASS::mvrnorm(
    n,
    mu = c(0, 0),
    Sigma = Sigma
  )

  cop <- copula::claytonCopula(arm$theta, dim = 2)
  UV <- copula::rCopula(n, cop)

  eta1 <- arm$eta_death_fun(X[, 1], X[, 2])
  eta2 <- arm$eta_prog_fun(X[, 1], X[, 2])

  T1 <- qexp(
    UV[, 1],
    rate = arm$lambda_death * exp(eta1)
  )

  T2 <- qexp(
    UV[, 2],
    rate = arm$lambda_prog * exp(eta2)
  )

  ## Exact censoring model used for the 2D LINEARITY test in the
  ## finalized Local-IPCW simulation.
  base_censor_eta <- 0.15 * X[, 1] - 0.10 * X[, 2]
  censor_eta <- if (A == 1L) base_censor_eta else -base_censor_eta

  C <- pmin(
    rexp(n, rate = censor_rate * exp(censor_eta)),
    tmax
  )

  make_observed_terminal_2d(
    x1 = X[, 1],
    x2 = X[, 2],
    T1 = T1,
    T2 = T2,
    C = C,
    A = A
  )
}

############################################################
## 2D GRID WEIGHTS, PLANE PROJECTION, AND TEST STATISTIC
############################################################

grid_weights_2d <- function(x_grid) {
  rep(1 / nrow(x_grid), nrow(x_grid))
}

fit_weighted_plane <- function(y, x_grid, w = NULL) {
  if (is.null(w)) w <- grid_weights_2d(x_grid)
  if (length(y) != nrow(x_grid) || length(w) != nrow(x_grid)) return(NULL)
  if (any(!is.finite(y)) || any(!is.finite(w))) return(NULL)

  D <- cbind(1, x_grid[, 1], x_grid[, 2])
  DtQD <- crossprod(D, w * D)
  if (qr(DtQD)$rank < 3L) return(NULL)

  beta <- solve(DtQD, crossprod(D, w * y))
  fitted <- as.vector(D %*% beta)

  list(
    beta = as.numeric(beta),
    fitted = fitted,
    resid = as.numeric(y - fitted)
  )
}

ph_linearity_stat_2d <- function(logWR, x_grid, n_total) {
  if (any(!is.finite(logWR))) return(NULL)
  w <- grid_weights_2d(x_grid)
  fit <- fit_weighted_plane(logWR, x_grid, w)
  if (is.null(fit)) return(NULL)

  list(
    T = as.numeric(n_total * sum(w * fit$resid^2)),
    beta = fit$beta,
    fitted = fit$fitted,
    resid = fit$resid,
    weights = w
  )
}

############################################################
## PH COVARIATE BASIS
## Common across H0/H11/H12 to keep the fitted model fixed.
############################################################

ph_basis_2d_test <- function(x1, x2) {
  cbind(
    x1 = x1,
    x2 = x2,
    x1_sq = x1^2,
    x2_sq = x2^2,
    x1_x2 = x1 * x2,
    x1_cu = x1^3,
    x2_cu = x2^3
  )
}

prepare_ph_arm <- function(dat, tmax = TMAX, m1 = M1, m2 = M2) {
  out <- dat
  X <- ph_basis_2d_test(out$x1, out$x2)
  storage.mode(X) <- "double"

  attr(out, "X_mat") <- X
  attr(out, "tmax") <- tmax

  ## Precompute Bernstein basis matrices because event/censoring times do not
  ## change across multiplier-bootstrap refits.
  attr(out, "B1_m")  <- bernstein_basis(out$U1, m = m1,       lower = 0, upper = tmax)
  attr(out, "B1_m1") <- bernstein_basis(out$U1, m = m1 - 1L, lower = 0, upper = tmax)
  attr(out, "B2_m")  <- bernstein_basis(out$U2, m = m2,       lower = 0, upper = tmax)
  attr(out, "B2_m1") <- bernstein_basis(out$U2, m = m2 - 1L, lower = 0, upper = tmax)
  out
}

ph_grid_basis <- function(x_grid) {
  ph_basis_2d_test(x_grid[, 1], x_grid[, 2])
}

############################################################
## BERNSTEIN BASELINE HAZARD UTILITIES
############################################################

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
  Bm <- bernstein_basis(t, m = m, lower = lower, upper = upper)
  Bm1 <- bernstein_basis(t, m = m - 1L, lower = lower, upper = upper)
  dphi <- diff(phi)

  Lambda <- as.vector(Bm %*% phi)
  lambda <- as.vector((m / (upper - lower)) * (Bm1 %*% dphi))

  list(
    Lambda = pmax(Lambda, 0),
    lambda = pmax(lambda, .Machine$double.xmin)
  )
}

baseline_eval_precomputed <- function(dat, phi, margin = 1L, upper = TMAX) {
  m <- length(phi) - 1L
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

############################################################
## CLAYTON CDF COPULA COMPONENTS
############################################################

clayton_cdf_components <- function(u, v, theta) {
  if (!is.finite(theta) || theta <= 0) stop("Clayton theta must be positive.")

  eps <- 1e-12
  u <- pmin(pmax(u, eps), 1 - eps)
  v <- pmin(pmax(v, eps), 1 - eps)

  A <- pmax(u^(-theta) + v^(-theta) - 1, eps)
  C <- A^(-1 / theta)
  Cu <- u^(-theta - 1) * A^(-1 / theta - 1)
  Cv <- v^(-theta - 1) * A^(-1 / theta - 1)
  density <- (1 + theta) * u^(-theta - 1) * v^(-theta - 1) *
    A^(-1 / theta - 2)

  list(
    C = pmax(C, eps),
    Cu = pmin(pmax(Cu, eps), 1 - eps),
    Cv = pmin(pmax(Cv, eps), 1 - eps),
    density = pmax(density, eps)
  )
}

############################################################
## PH--COPULA PARAMETERIZATION
############################################################

theta_from_zeta_clayton <- function(zeta) exp(zeta)
zeta_from_theta_clayton <- function(theta) {
  if (theta <= 0) stop("Clayton theta must be positive.")
  log(theta)
}

unpack_ph_par <- function(par, p, m1, m2) {
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
    phi1 = phi_from_eta(eta1),
    phi2 = phi_from_eta(eta2)
  )
}

make_ph_start <- function(dat, m1, m2, theta0) {
  X <- attr(dat, "X_mat")
  p <- ncol(X)
  upper <- attr(dat, "tmax")

  lam1 <- max(sum(dat$Delta1) / max(sum(dat$U1), 1e-8), 1e-3)
  lam2 <- max(sum(dat$Delta2) / max(sum(dat$U2), 1e-8), 1e-3)
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

############################################################
## INDIVIDUAL CDF-CLAYTON LOG-LIKELIHOOD
############################################################

individual_loglik_clayton_cdf <- function(par, dat, m1 = M1, m2 = M2) {
  X <- attr(dat, "X_mat")
  p <- ncol(X)
  upper <- attr(dat, "tmax")
  pars <- unpack_ph_par(par, p, m1, m2)

  lp1 <- as.vector(X %*% pars$beta1)
  lp2 <- as.vector(X %*% pars$beta2)

  b1 <- baseline_eval_precomputed(dat, pars$phi1, margin = 1L, upper = upper)
  b2 <- baseline_eval_precomputed(dat, pars$phi2, margin = 2L, upper = upper)

  logS1 <- -exp(lp1) * b1$Lambda
  logS2 <- -exp(lp2) * b2$Lambda
  S1 <- exp(pmax(logS1, -745))
  S2 <- exp(pmax(logS2, -745))
  F1 <- 1 - S1
  F2 <- 1 - S2

  logf1 <- log(b1$lambda) + lp1 + logS1
  logf2 <- log(b2$lambda) + lp2 + logS2

  cc <- clayton_cdf_components(F1, F2, pars$theta)

  eps <- 1e-12
  d1 <- dat$Delta1
  d2 <- dat$Delta2

  both_events <- d1 == 1L & d2 == 1L
  only1 <- d1 == 1L & d2 == 0L
  only2 <- d1 == 0L & d2 == 1L
  both_cens <- d1 == 0L & d2 == 0L

  ## Joint survival under a CDF copula:
  ## J(t1,t2)=S1(t1)+S2(t2)-1+C(F1(t1),F2(t2)).
  J <- pmax(S1 + S2 - 1 + cc$C, eps)

  ll <- numeric(nrow(dat))
  ll[both_events] <-
    logf1[both_events] + logf2[both_events] + log(cc$density[both_events])
  ll[only1] <-
    logf1[only1] + log(pmax(1 - cc$Cu[only1], eps))
  ll[only2] <-
    logf2[only2] + log(pmax(1 - cc$Cv[only2], eps))
  ll[both_cens] <- log(J[both_cens])

  ll
}

############################################################
## PH--COPULA FITTING
############################################################

fit_phcopula_arm <- function(
    dat,
    weights = NULL,
    start = NULL,
    multistart = TRUE,
    maxit = ORIG_MAXIT,
    reltol = OPT_RELTOL) {

  n <- nrow(dat)
  if (is.null(weights)) weights <- rep(1, n)

  if (length(weights) != n || any(!is.finite(weights)) || any(weights <= 0)) {
    stop("Invalid likelihood weights.")
  }

  ## Common rescaling leaves the weighted-likelihood maximizer unchanged.
  weights <- weights / mean(weights)

  objective <- function(par) {
    ll <- individual_loglik_clayton_cdf(par, dat, m1 = M1, m2 = M2)
    if (any(!is.finite(ll))) return(1e100)
    val <- sum(weights * ll)
    if (!is.finite(val)) return(1e100)
    -val
  }

  if (!is.null(start)) {
    starts <- list(start)
  } else {
    starts <- lapply(c(0.5, 1.2, 2.0), function(th) {
      make_ph_start(dat, M1, M2, theta0 = th)
    })
  }

  if (!multistart) starts <- starts[1L]

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

  ## One Nelder-Mead fallback if the best BFGS fit did not converge.
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

  if (is.null(best)) return(NULL)

  p <- ncol(attr(dat, "X_mat"))
  pars <- unpack_ph_par(best$par, p, M1, M2)

  list(
    par_unconstr = best$par,
    convergence = best$convergence,
    negloglik = best$value,
    beta1 = pars$beta1,
    beta2 = pars$beta2,
    theta = pars$theta,
    phi1 = pars$phi1,
    phi2 = pars$phi2,
    p = p,
    m1 = M1,
    m2 = M2,
    lower = 0,
    upper = attr(dat, "tmax"),
    family = "clayton_cdf"
  )
}

############################################################
## FAST FITTED WIN-RATIO SURFACE ON THE 5 x 5 GRID
############################################################

simpson_rule <- function(tau, n = N_INT) {
  if (n %% 2L != 1L) stop("Simpson grid size must be odd.")

  t <- seq(0, tau, length.out = n)
  h <- tau / (n - 1L)
  w <- rep(2, n)
  w[seq(2L, n - 1L, by = 2L)] <- 4
  w[c(1L, n)] <- 1

  list(t = t, w = w * h / 3)
}

margin_grid <- function(fit, margin, Zg, t) {
  phi <- if (margin == 1L) fit$phi1 else fit$phi2
  beta <- if (margin == 1L) fit$beta1 else fit$beta2

  b <- baseline_eval(t, phi, lower = 0, upper = fit$upper)
  lp <- as.vector(Zg %*% beta)
  rr <- exp(lp)

  ## Rows = time; columns = covariate-grid points.
  cumhaz <- outer(b$Lambda, rr, "*")
  S <- exp(-cumhaz)
  f <- outer(b$lambda, rr, "*") * S

  list(S = S, f = f)
}

copula_matrix_clayton_cdf <- function(U, V, theta, component = c("C", "Cv")) {
  component <- match.arg(component)
  d <- dim(U)
  cc <- clayton_cdf_components(U, V, theta)
  matrix(cc[[component]], nrow = d[1], ncol = d[2])
}

logwr_surface_from_fits <- function(
    fit1,
    fit0,
    x_grid,
    tau = TAU,
    n_int = N_INT) {

  Zg <- ph_grid_basis(x_grid)
  rule <- simpson_rule(tau, n = n_int)
  tt <- rule$t
  ww <- rule$w

  m11 <- margin_grid(fit1, 1L, Zg, tt)
  m10 <- margin_grid(fit0, 1L, Zg, tt)
  m21 <- margin_grid(fit1, 2L, Zg, tt)
  m20 <- margin_grid(fit0, 2L, Zg, tt)

  ## Primary component.
  Wp <- colSums((m11$S * m10$f) * ww)
  Lp <- colSums((m10$S * m11$f) * ww)

  ## Secondary component conditional on both subjects surviving the primary
  ## event beyond tau.
  S11_tau <- m11$S[nrow(m11$S), ]
  S10_tau <- m10$S[nrow(m10$S), ]

  S11_mat <- matrix(
    S11_tau,
    nrow = length(tt),
    ncol = length(S11_tau),
    byrow = TRUE
  )
  S10_mat <- matrix(
    S10_tau,
    nrow = length(tt),
    ncol = length(S10_tau),
    byrow = TRUE
  )

  F11_tau <- 1 - S11_mat
  F10_tau <- 1 - S10_mat
  F21 <- 1 - m21$S
  F20 <- 1 - m20$S

  C1 <- copula_matrix_clayton_cdf(F11_tau, F21, fit1$theta, "C")
  C0 <- copula_matrix_clayton_cdf(F10_tau, F20, fit0$theta, "C")
  Cv1 <- copula_matrix_clayton_cdf(F11_tau, F21, fit1$theta, "Cv")
  Cv0 <- copula_matrix_clayton_cdf(F10_tau, F20, fit0$theta, "Cv")

  J1 <- pmax(S11_mat + m21$S - 1 + C1, 0)
  J0 <- pmax(S10_mat + m20$S - 1 + C0, 0)

  ## -dJ/dt = f2(t) {1 - C_v(F1(tau),F2(t))}.
  minus_dJ1 <- (1 - Cv1) * m21$f
  minus_dJ0 <- (1 - Cv0) * m20$f

  Ws <- colSums((J1 * minus_dJ0) * ww)
  Ls <- colSums((J0 * minus_dJ1) * ww)

  W <- Wp + Ws
  L <- Lp + Ls

  bad <- !is.finite(W) | !is.finite(L) | W <= 0 | L <= 0
  out <- log(W / L)
  out[bad] <- NA_real_
  out
}

############################################################
## MULTIPLIER-BOOTSTRAP 2D LOG-LINEARITY TEST
############################################################

phcopula_linearity_test_2d <- function(
    dat1,
    dat0,
    x_grid,
    tau,
    B,
    alpha,
    seed) {

  d1 <- prepare_ph_arm(dat1, tmax = TMAX, m1 = M1, m2 = M2)
  d0 <- prepare_ph_arm(dat0, tmax = TMAX, m1 = M1, m2 = M2)

  fit1 <- fit_phcopula_arm(
    d1,
    multistart = TRUE,
    maxit = ORIG_MAXIT,
    reltol = OPT_RELTOL
  )

  fit0 <- fit_phcopula_arm(
    d0,
    multistart = TRUE,
    maxit = ORIG_MAXIT,
    reltol = OPT_RELTOL
  )

  if (is.null(fit1) || is.null(fit0) ||
      fit1$convergence != 0L || fit0$convergence != 0L) {
    return(NULL)
  }

  logWR_hat <- logwr_surface_from_fits(
    fit1,
    fit0,
    x_grid = x_grid,
    tau = tau
  )

  if (any(!is.finite(logWR_hat))) return(NULL)

  n_total <- nrow(d1) + nrow(d0)
  obs <- ph_linearity_stat_2d(logWR_hat, x_grid, n_total = n_total)
  if (is.null(obs)) return(NULL)

  set.seed(seed)
  T_boot <- rep(NA_real_, B)

  for (b in seq_len(B)) {
    w1 <- rexp(nrow(d1), rate = 1)
    w0 <- rexp(nrow(d0), rate = 1)

    fb1 <- try(
      fit_phcopula_arm(
        d1,
        weights = w1,
        start = fit1$par_unconstr,
        multistart = FALSE,
        maxit = BOOT_MAXIT,
        reltol = OPT_RELTOL
      ),
      silent = TRUE
    )

    fb0 <- try(
      fit_phcopula_arm(
        d0,
        weights = w0,
        start = fit0$par_unconstr,
        multistart = FALSE,
        maxit = BOOT_MAXIT,
        reltol = OPT_RELTOL
      ),
      silent = TRUE
    )

    if (inherits(fb1, "try-error") || inherits(fb0, "try-error") ||
        is.null(fb1) || is.null(fb0) ||
        fb1$convergence != 0L || fb0$convergence != 0L) {
      next
    }

    logWR_b <- try(
      logwr_surface_from_fits(
        fb1,
        fb0,
        x_grid = x_grid,
        tau = tau
      ),
      silent = TRUE
    )

    if (inherits(logWR_b, "try-error") || any(!is.finite(logWR_b))) {
      next
    }

    ## Centered multiplier process required by the test theory.
    delta_logWR <- logWR_b - logWR_hat
    bst <- ph_linearity_stat_2d(
      delta_logWR,
      x_grid,
      n_total = n_total
    )

    if (!is.null(bst) && is.finite(bst$T)) {
      T_boot[b] <- bst$T
    }
  }

  T_boot_fin <- T_boot[is.finite(T_boot)]
  n_success <- length(T_boot_fin)

  ## Correct threshold for pilot/check/final runs.
  min_success <- if (B <= 20L) {
    max(3L, ceiling(0.5 * B))
  } else {
    max(20L, ceiling(0.5 * B))
  }

  if (n_success < min_success) {
    return(list(
      T_obs = obs$T,
      crit = NA_real_,
      p_value = NA_real_,
      reject = NA_real_,
      n_boot_success = n_success,
      theta_copula_arm1 = fit1$theta,
      theta_copula_arm0 = fit0$theta,
      plane_intercept = obs$beta[1],
      plane_x1 = obs$beta[2],
      plane_x2 = obs$beta[3]
    ))
  }

  crit <- as.numeric(
    quantile(T_boot_fin, probs = 1 - alpha, names = FALSE)
  )

  p_value <- mean(T_boot_fin >= obs$T)
  reject <- as.numeric(obs$T > crit)

  list(
    T_obs = obs$T,
    crit = crit,
    p_value = p_value,
    reject = reject,
    n_boot_success = n_success,
    theta_copula_arm1 = fit1$theta,
    theta_copula_arm0 = fit0$theta,
    plane_intercept = obs$beta[1],
    plane_x1 = obs$beta[2],
    plane_x2 = obs$beta[3]
  )
}

############################################################
## TRUE POPULATION WR SURFACE: ONE-TIME DIAGNOSTIC
############################################################

true_logwr_surface_2d <- function(
    scenario,
    x_grid = GRID_2D,
    tau = TAU,
    n_int = N_INT_POP) {

  dgp <- make_linearity_dgp_2d(scenario)
  rule <- simpson_rule(tau, n = n_int)
  tt <- rule$t
  ww <- rule$w
  K <- nrow(x_grid)

  x1 <- x_grid[, 1]
  x2 <- x_grid[, 2]

  rate11 <- dgp$arm1$lambda_death * exp(dgp$arm1$eta_death_fun(x1, x2))
  rate21 <- dgp$arm1$lambda_prog  * exp(dgp$arm1$eta_prog_fun(x1, x2))
  rate10 <- dgp$arm0$lambda_death * exp(dgp$arm0$eta_death_fun(x1, x2))
  rate20 <- dgp$arm0$lambda_prog  * exp(dgp$arm0$eta_prog_fun(x1, x2))

  S11 <- exp(-outer(tt, rate11, "*"))
  S21 <- exp(-outer(tt, rate21, "*"))
  S10 <- exp(-outer(tt, rate10, "*"))
  S20 <- exp(-outer(tt, rate20, "*"))

  f11 <- outer(rep(1, length(tt)), rate11, "*") * S11
  f21 <- outer(rep(1, length(tt)), rate21, "*") * S21
  f10 <- outer(rep(1, length(tt)), rate10, "*") * S10
  f20 <- outer(rep(1, length(tt)), rate20, "*") * S20

  Wp <- colSums((S11 * f10) * ww)
  Lp <- colSums((S10 * f11) * ww)

  S11_tau <- exp(-tau * rate11)
  S10_tau <- exp(-tau * rate10)
  S11_mat <- matrix(S11_tau, nrow = length(tt), ncol = K, byrow = TRUE)
  S10_mat <- matrix(S10_tau, nrow = length(tt), ncol = K, byrow = TRUE)

  F11_tau <- 1 - S11_mat
  F10_tau <- 1 - S10_mat
  F21 <- 1 - S21
  F20 <- 1 - S20

  C1 <- copula_matrix_clayton_cdf(F11_tau, F21, dgp$arm1$theta, "C")
  C0 <- copula_matrix_clayton_cdf(F10_tau, F20, dgp$arm0$theta, "C")
  Cv1 <- copula_matrix_clayton_cdf(F11_tau, F21, dgp$arm1$theta, "Cv")
  Cv0 <- copula_matrix_clayton_cdf(F10_tau, F20, dgp$arm0$theta, "Cv")

  J1 <- pmax(S11_mat + S21 - 1 + C1, 0)
  J0 <- pmax(S10_mat + S20 - 1 + C0, 0)

  minus_dJ1 <- (1 - Cv1) * f21
  minus_dJ0 <- (1 - Cv0) * f20

  Ws <- colSums((J1 * minus_dJ0) * ww)
  Ls <- colSums((J0 * minus_dJ1) * ww)

  W <- Wp + Ws
  L <- Lp + Ls
  out <- log(W / L)
  out[!is.finite(W) | !is.finite(L) | W <= 0 | L <= 0] <- NA_real_
  out
}

population_surface_check <- function() {
  scenarios <- c("H0", "mild", "strong")

  rows <- lapply(scenarios, function(sc) {
    ell <- true_logwr_surface_2d(sc)
    w <- grid_weights_2d(GRID_2D)
    fit <- fit_weighted_plane(ell, GRID_2D, w)

    data.frame(
      scenario = sc,
      plane_intercept = fit$beta[1],
      plane_x1 = fit$beta[2],
      plane_x2 = fit$beta[3],
      max_abs_nonlinear_resid = max(abs(fit$resid)),
      integrated_sq_nonlinearity = sum(w * fit$resid^2),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

############################################################
## SCENARIO DEFINITIONS
##
## local_scenario_row preserves rows 4/5/6 of the original 2D Local-IPCW
## scenario table so that the generated DATA use exactly the same seeds.
############################################################

scenario_table <- data.frame(
  scenario_id = c("L2_H0", "L2_MILD", "L2_STRONG"),
  hypothesis = c("H0", "H11", "H12"),
  scenario = c("H0", "mild", "strong"),
  local_scenario_row = c(4L, 5L, 6L),
  stringsAsFactors = FALSE
)

############################################################
## ONE MONTE CARLO REPLICATION
############################################################

run_one_rep_ph_test_2d <- function(r, scenario_def, settings) {
  time_start <- Sys.time()

  out <- list(
    r = r,
    reject = NA_real_,
    p_value = NA_real_,
    T_obs = NA_real_,
    crit = NA_real_,
    n_boot_success = NA_real_,
    copula_theta1 = NA_real_,
    copula_theta0 = NA_real_,
    plane_intercept = NA_real_,
    plane_x1 = NA_real_,
    plane_x2 = NA_real_,
    elapsed_min = NA_real_,
    error = NULL
  )

  dgp <- try(
    make_linearity_dgp_2d(scenario_def$scenario),
    silent = TRUE
  )

  if (inherits(dgp, "try-error")) {
    out$error <- "2D linearity DGP construction failed"
    return(out)
  }

  seed_base <- settings$seed + r

  dat1 <- try(
    simulate_arm_2d(
      n = settings$n1,
      arm = dgp$arm1,
      A = 1L,
      seed = seed_base + 1000L,
      tmax = settings$tmax,
      censor_rate = settings$censor_rate,
      rho_x = settings$rho_x
    ),
    silent = TRUE
  )

  dat0 <- try(
    simulate_arm_2d(
      n = settings$n0,
      arm = dgp$arm0,
      A = 0L,
      seed = seed_base + 2000L,
      tmax = settings$tmax,
      censor_rate = settings$censor_rate,
      rho_x = settings$rho_x
    ),
    silent = TRUE
  )

  if (inherits(dat1, "try-error") || inherits(dat0, "try-error")) {
    out$error <- "2D simulation failed"
    out$elapsed_min <- as.numeric(difftime(Sys.time(), time_start, units = "mins"))
    return(out)
  }

  tst <- try(
    phcopula_linearity_test_2d(
      dat1 = dat1,
      dat0 = dat0,
      x_grid = settings$grid_2d,
      tau = settings$tau,
      B = settings$B,
      alpha = settings$alpha,
      seed = seed_base + 500000L
    ),
    silent = TRUE
  )

  if (inherits(tst, "try-error") || is.null(tst)) {
    out$error <- "PH--copula fit/bootstrap test failed"
    out$elapsed_min <- as.numeric(difftime(Sys.time(), time_start, units = "mins"))
    return(out)
  }

  out$reject <- tst$reject
  out$p_value <- tst$p_value
  out$T_obs <- tst$T_obs
  out$crit <- tst$crit
  out$n_boot_success <- tst$n_boot_success
  out$copula_theta1 <- tst$theta_copula_arm1
  out$copula_theta0 <- tst$theta_copula_arm0
  out$plane_intercept <- tst$plane_intercept
  out$plane_x1 <- tst$plane_x1
  out$plane_x2 <- tst$plane_x2
  out$elapsed_min <- as.numeric(difftime(Sys.time(), time_start, units = "mins"))

  if (!is.finite(out$reject)) {
    out$error <- paste0(
      "Insufficient successful multiplier fits: ",
      out$n_boot_success, "/", settings$B
    )
  }

  out
}

############################################################
## SUMMARIZE ONE SCENARIO x SAMPLE-SIZE SETTING
############################################################

safe_mean <- function(x) {
  if (any(is.finite(x))) mean(x, na.rm = TRUE) else NA_real_
}

safe_median <- function(x) {
  if (any(is.finite(x))) median(x, na.rm = TRUE) else NA_real_
}

safe_min <- function(x) {
  if (any(is.finite(x))) min(x, na.rm = TRUE) else NA_real_
}

summarize_scenario <- function(results, scenario_def, settings) {
  get_num <- function(name) {
    vapply(
      results,
      function(z) {
        if (is.null(z) || is.null(z[[name]])) NA_real_ else as.numeric(z[[name]])
      },
      numeric(1)
    )
  }

  reject <- get_num("reject")
  p_value <- get_num("p_value")
  T_obs <- get_num("T_obs")
  crit <- get_num("crit")
  n_boot_success <- get_num("n_boot_success")
  ctheta1 <- get_num("copula_theta1")
  ctheta0 <- get_num("copula_theta0")
  intercept <- get_num("plane_intercept")
  slope1 <- get_num("plane_x1")
  slope2 <- get_num("plane_x2")
  elapsed_min <- get_num("elapsed_min")

  error_vec <- vapply(
    results,
    function(z) {
      if (is.null(z) || is.null(z$error)) "" else as.character(z$error)
    },
    character(1)
  )

  good <- is.finite(reject)
  n_valid <- sum(good)

  rejection_rate <- if (n_valid > 0L) mean(reject[good]) else NA_real_
  mcse_rejection <- if (n_valid > 0L && is.finite(rejection_rate)) {
    sqrt(rejection_rate * (1 - rejection_rate) / n_valid)
  } else {
    NA_real_
  }

  data.frame(
    scenario_id = scenario_def$scenario_id,
    hypothesis = scenario_def$hypothesis,
    scenario = scenario_def$scenario,
    n_per_arm = settings$n1,
    tau = settings$tau,
    R = settings$R,
    B = settings$B,
    rejection_rate = rejection_rate,
    mcse_rejection = mcse_rejection,
    mean_p_value = safe_mean(p_value),
    median_p_value = safe_median(p_value),
    mean_T_obs = safe_mean(T_obs),
    mean_crit = safe_mean(crit),
    mean_boot_success = safe_mean(n_boot_success),
    min_boot_success = safe_min(n_boot_success),
    mean_copula_theta_arm1 = safe_mean(ctheta1),
    mean_copula_theta_arm0 = safe_mean(ctheta0),
    mean_plane_intercept = safe_mean(intercept),
    mean_plane_x1 = safe_mean(slope1),
    mean_plane_x2 = safe_mean(slope2),
    n_valid = n_valid,
    error_count = sum(nzchar(error_vec)),
    mean_elapsed_min = safe_mean(elapsed_min),
    stringsAsFactors = FALSE
  )
}

############################################################
## PARALLEL RUNNER WITH CHECKPOINTING
############################################################

run_scenario_parallel <- function(
    scenario_def,
    settings,
    n_cores,
    batch_size,
    progress_every,
    resume = TRUE) {

  scenario_id <- scenario_def$scenario_id
  n_per_arm <- settings$n1

  checkpoint_file <- file.path(
    CHECKPOINT_DIR,
    paste0(scenario_id, "_n", n_per_arm, "_checkpoint.rds")
  )

  results <- vector("list", settings$R)

  if (resume && file.exists(checkpoint_file)) {
    old <- try(readRDS(checkpoint_file), silent = TRUE)

    if (!inherits(old, "try-error") &&
        is.list(old) &&
        identical(old$scenario_id, scenario_id) &&
        identical(old$design_version, DESIGN_VERSION) &&
        identical(old$settings, settings) &&
        length(old$results) == settings$R) {
      results <- old$results
      cat("Resuming", scenario_id, "n =", n_per_arm, "from checkpoint.\n")
    }
  }

  pending <- which(vapply(results, is.null, logical(1)))

  if (length(pending) == 0L) {
    cat(scenario_id, "n =", n_per_arm, "already complete.\n")
    return(list(
      summary = summarize_scenario(results, scenario_def, settings),
      results = results
    ))
  }

  cat("\n============================================================\n")
  cat(
    "Running", scenario_id,
    "|", scenario_def$hypothesis,
    "| n/arm =", n_per_arm, "\n"
  )
  cat(
    "PH--copula 2D log-linearity test | Clayton CDF | common cubic PH basis\n"
  )
  cat(
    "R =", settings$R,
    "| B =", settings$B,
    "| tau =", settings$tau,
    "| cores =", n_cores, "\n"
  )
  cat("============================================================\n")

  cl <- parallel::makeCluster(n_cores, type = "PSOCK")
  on.exit(parallel::stopCluster(cl), add = TRUE)

  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(copula)
      library(MASS)
    })
    NULL
  })

  export_names <- c(
    "TAU", "TMAX", "ALPHA", "CENSOR_RATE", "RHO_X",
    "GRID_2D_VALUES", "GRID_2D",
    "M1", "M2", "N_INT", "ORIG_MAXIT", "BOOT_MAXIT", "OPT_RELTOL",
    "make_observed_terminal_2d", "make_linearity_dgp_2d", "simulate_arm_2d",
    "grid_weights_2d", "fit_weighted_plane", "ph_linearity_stat_2d",
    "ph_basis_2d_test", "prepare_ph_arm", "ph_grid_basis",
    "bernstein_basis", "phi_from_eta", "baseline_eval", "baseline_eval_precomputed",
    "clayton_cdf_components", "theta_from_zeta_clayton", "zeta_from_theta_clayton",
    "unpack_ph_par", "make_ph_start", "individual_loglik_clayton_cdf",
    "fit_phcopula_arm", "simpson_rule", "margin_grid",
    "copula_matrix_clayton_cdf", "logwr_surface_from_fits",
    "phcopula_linearity_test_2d", "run_one_rep_ph_test_2d"
  )

  parallel::clusterExport(cl, export_names, envir = .GlobalEnv)

  batches <- split(pending, ceiling(seq_along(pending) / batch_size))
  last_reported <- sum(!vapply(results, is.null, logical(1)))

  for (idx in batches) {
    ans <- parallel::parLapplyLB(
      cl,
      idx,
      run_one_rep_ph_test_2d,
      scenario_def = scenario_def,
      settings = settings
    )

    for (k in seq_along(idx)) {
      results[[idx[k]]] <- ans[[k]]
    }

    saveRDS(
      list(
        scenario_id = scenario_id,
        design_version = DESIGN_VERSION,
        settings = settings,
        results = results
      ),
      checkpoint_file
    )

    completed <- sum(!vapply(results, is.null, logical(1)))
    report_now <-
      completed == settings$R ||
      floor(completed / progress_every) > floor(last_reported / progress_every)

    if (report_now) {
      cat(
        "Completed", completed, "/", settings$R,
        "for", scenario_id,
        "n =", n_per_arm, "\n"
      )
      flush.console()
      last_reported <- completed
    }
  }

  parallel::stopCluster(cl)
  on.exit(NULL, add = FALSE)

  list(
    summary = summarize_scenario(results, scenario_def, settings),
    results = results
  )
}

############################################################
## POPULATION-SURFACE CHECK
############################################################

pop_check <- population_surface_check()

cat("\n============================================================\n")
cat("POPULATION SURFACE CHECK FOR 2D LINEARITY DGPs\n")
cat("============================================================\n")
pop_print <- pop_check
num_pop <- vapply(pop_print, is.numeric, logical(1))
pop_print[num_pop] <- lapply(pop_print[num_pop], function(z) signif(z, 6))
print(pop_print, row.names = FALSE)

############################################################
## MASTER RUN TABLE
############################################################

run_table <- expand.grid(
  n_per_arm = SAMPLE_SIZES,
  scenario_row = seq_len(nrow(scenario_table)),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

run_table$scenario_id <- scenario_table$scenario_id[run_table$scenario_row]
run_table$hypothesis <- scenario_table$hypothesis[run_table$scenario_row]
run_table$scenario <- scenario_table$scenario[run_table$scenario_row]
run_table$local_scenario_row <- scenario_table$local_scenario_row[run_table$scenario_row]

hyp_order <- match(run_table$hypothesis, c("H0", "H11", "H12"))
run_table <- run_table[
  order(run_table$n_per_arm, hyp_order),
  ,
  drop = FALSE
]
rownames(run_table) <- NULL

############################################################
## RUN PH--COPULA 2D LOG-LINEARITY SIMULATION
############################################################

cat("\n============================================================\n")
cat("PH--COPULA 2D LOG-LINEARITY TEST SIMULATION\n")
cat("============================================================\n")
cat("Run mode:", RUN_MODE, "\n")
cat("Sample sizes per arm:", paste(SAMPLE_SIZES, collapse = ", "), "\n")
cat("R =", R_MC, "| B =", B_BOOT, "| tau =", TAU, "| alpha =", ALPHA, "\n")
cat("Grid: 5 x 5 =", nrow(GRID_2D), "points over [-1,1]^2\n")
cat("Corr(X1,X2) =", RHO_X, "\n")
cat("PH basis: x1 + x2 + x1^2 + x2^2 + x1*x2 + x1^3 + x2^3\n")
cat("Copula: Clayton CDF = rotated Clayton survival\n")
cat("Multiplier bootstrap: independent Exp(1) weights within each arm\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("============================================================\n\n")

all_results <- vector("list", nrow(run_table))
summary_list <- vector("list", nrow(run_table))
run_keys <- paste0(run_table$scenario_id, "_n", run_table$n_per_arm)
names(all_results) <- run_keys

for (ii in seq_len(nrow(run_table))) {
  rr <- run_table[ii, , drop = FALSE]
  def <- as.list(scenario_table[rr$scenario_row, , drop = FALSE])
  n_now <- as.integer(rr$n_per_arm)

  ## EXACT seed construction used by the finalized 2D Local-IPCW program.
  ## Its linearity scenarios occupied rows 4, 5, and 6.
  settings <- list(
    R = R_MC,
    B = B_BOOT,
    n0 = n_now,
    n1 = n_now,
    tau = TAU,
    tmax = TMAX,
    alpha = ALPHA,
    censor_rate = CENSOR_RATE,
    rho_x = RHO_X,
    grid_2d = GRID_2D,
    design_version = DESIGN_VERSION,
    seed = BASE_SEED + n_now * 1000000L +
      as.integer(rr$local_scenario_row) * 100000L
  )

  z <- run_scenario_parallel(
    scenario_def = def,
    settings = settings,
    n_cores = N_CORES,
    batch_size = BATCH_SIZE,
    progress_every = PROGRESS_EVERY,
    resume = RESUME
  )

  key <- paste0(rr$scenario_id, "_n", n_now)
  all_results[[key]] <- z$results
  summary_list[[ii]] <- z$summary

  saveRDS(
    list(
      design_version = DESIGN_VERSION,
      run_mode = RUN_MODE,
      run_table = run_table,
      population_check = pop_check,
      completed_results = all_results,
      completed_summaries = summary_list
    ),
    file.path(OUTPUT_DIR, "partial_state.rds")
  )
}

############################################################
## FINAL SUMMARY TABLES
############################################################

final_summary <- do.call(rbind, summary_list)
rownames(final_summary) <- NULL

hyp_ord <- match(final_summary$hypothesis, c("H0", "H11", "H12"))
final_summary <- final_summary[
  order(final_summary$n_per_arm, hyp_ord),
  ,
  drop = FALSE
]
rownames(final_summary) <- NULL

wide_rows <- lapply(SAMPLE_SIZES, function(n_now) {
  d <- final_summary[final_summary$n_per_arm == n_now, , drop = FALSE]

  get_rate <- function(h) {
    z <- d$rejection_rate[d$hypothesis == h]
    if (length(z) != 1L) NA_real_ else z
  }

  data.frame(
    dimension = "2D",
    method = "PH--copula",
    n_per_arm = n_now,
    H0 = get_rate("H0"),
    H11 = get_rate("H11"),
    H12 = get_rate("H12"),
    stringsAsFactors = FALSE
  )
})

table_wide <- do.call(rbind, wide_rows)
rownames(table_wide) <- NULL

cat("\n\n============================================================\n")
cat("PH--COPULA 2D LOG-LINEARITY TEST RESULTS\n")
cat("============================================================\n")
summary_print <- final_summary
num <- vapply(summary_print, is.numeric, logical(1))
summary_print[num] <- lapply(summary_print[num], function(z) round(z, 4))
print(summary_print, row.names = FALSE)

cat("\n\n============================================================\n")
cat("MANUSCRIPT-STYLE PH--COPULA 2D REJECTION-RATE TABLE\n")
cat("============================================================\n")
wide_print <- table_wide
num <- vapply(wide_print, is.numeric, logical(1))
wide_print[num] <- lapply(wide_print[num], function(z) round(z, 3))
print(wide_print, row.names = FALSE)

cat("\n\n============================================================\n")
cat("DIAGNOSTICS\n")
cat("============================================================\n")
for (ii in seq_len(nrow(final_summary))) {
  z <- final_summary[ii, , drop = FALSE]
  cat(
    "\n", z$scenario_id,
    " | n/arm = ", z$n_per_arm,
    " | ", z$hypothesis, "\n",
    sep = ""
  )
  cat("  n_valid               =", z$n_valid, "\n")
  cat("  error_count           =", z$error_count, "\n")
  cat("  mean_boot_success     =", round(z$mean_boot_success, 3), "\n")
  cat("  min_boot_success      =", z$min_boot_success, "\n")
  cat("  mean Clayton theta A1 =", round(z$mean_copula_theta_arm1, 4), "\n")
  cat("  mean Clayton theta A0 =", round(z$mean_copula_theta_arm0, 4), "\n")
}

############################################################
## SAVE
############################################################

summary_file <- file.path(
  OUTPUT_DIR,
  paste0("phcopula_2d_loglinearity_", RUN_MODE, "_summary.csv")
)
wide_file <- file.path(
  OUTPUT_DIR,
  paste0("phcopula_2d_loglinearity_", RUN_MODE, "_table_wide.csv")
)
pop_file <- file.path(
  OUTPUT_DIR,
  paste0("phcopula_2d_loglinearity_", RUN_MODE, "_population_check.csv")
)
rds_file <- file.path(
  OUTPUT_DIR,
  paste0("phcopula_2d_loglinearity_", RUN_MODE, "_results.rds")
)

write.csv(final_summary, summary_file, row.names = FALSE)
write.csv(table_wide, wide_file, row.names = FALSE)
write.csv(pop_check, pop_file, row.names = FALSE)

saveRDS(
  list(
    design_version = DESIGN_VERSION,
    run_mode = RUN_MODE,
    settings = list(
      R = R_MC,
      B = B_BOOT,
      sample_sizes = SAMPLE_SIZES,
      tau = TAU,
      tmax = TMAX,
      alpha = ALPHA,
      censor_rate = CENSOR_RATE,
      rho_x = RHO_X,
      grid_2d = GRID_2D,
      m1 = M1,
      m2 = M2,
      n_int = N_INT,
      ph_basis = "x1 + x2 + x1^2 + x2^2 + x1*x2 + x1^3 + x2^3",
      copula = "Clayton CDF / rotated Clayton survival",
      original_maxit = ORIG_MAXIT,
      bootstrap_maxit = BOOT_MAXIT,
      n_cores = N_CORES
    ),
    scenario_table = scenario_table,
    run_table = run_table,
    population_check = pop_check,
    summary = final_summary,
    table_wide = table_wide,
    results = all_results
  ),
  rds_file
)

cat("\nSaved files:\n")
cat("  ", normalizePath(summary_file, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("  ", normalizePath(wide_file, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("  ", normalizePath(pop_file, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("  ", normalizePath(rds_file, winslash = "/", mustWork = FALSE), "\n", sep = "")

cat("\nTotal elapsed time:\n")
print(Sys.time() - time0)

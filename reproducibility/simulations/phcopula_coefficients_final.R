############################################################
## PH--COPULA-ONLY COEFFICIENT SIMULATION FOR TABLE 2
##
## PURPOSE
## -------
## Add PH--copula columns to the existing Table 2 WITHOUT
## rerunning Mao--Wang or Local IPCW.
##
## The data-generating mechanisms, truth values, sample sizes,
## seeds, horizon, censoring, and 5 x 5 raw-covariate grid are
## exactly those in the finalized coefficient simulation.
##
## PH--copula working model used here
## ---------------------------------
## The marginal PH covariate basis is matched to the DGP while remaining
## correctly specified in every scenario:
##   S1/S3: z(x1,x2) = (x1, x2)
##   S2/S4: z(x1,x2) = (x1, x2, x1^2, x1*x2, x2^2)
## Thus S1/S3 avoid unnecessary quadratic nuisance terms, whereas S2/S4
## retain the nonlinear terms required by their generating hazards.
##
## Copula family is matched to the DGP:
##   S1/S3: Gumbel SURVIVAL copula.
##   S2/S4: the original DGP generated event times as qexp(U, rate)
##          with a Clayton CDF copula.  Therefore the corresponding
##          survival copula is the 180-degree rotated Clayton copula.
##          The code fits this model through the equivalent Clayton
##          CDF-copula likelihood and uses the implied joint survival
##          function when computing the win ratio.
##
## Inference
## ---------
## Arm-specific Exp(1) weighted multiplier bootstrap.
## In each bootstrap sample, the PH--copula model is refitted from
## the original MLE as a warm start.  The best linear projection of
## log WR on (1, x1, x2) is recomputed on the same 5 x 5 grid.
##
## Output columns for Table 2:
##   Bias, empirical SE, average bootstrap SE (SEE), and coverage (CP).
##
## IMPORTANT
## ---------
## This program does NOT run Local IPCW and does NOT run Mao--Wang.
##
## Final manuscript run uses RUN_MODE = "final": R=500, B=200,
## all four scenarios, with checkpoint/resume enabled.  Set RUN_MODE
## to "pilot" only if a smaller diagnostic rerun is desired.
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

RUN_MODE <- "final"  # "pilot" or "final"
if (!RUN_MODE %in% c("pilot", "final")) {
  stop("RUN_MODE must be 'pilot' or 'final'.")
}

N0 <- 250L
N1 <- 250L
TAU <- 8
TMAX <- 10
CENSOR_RATE <- 0.05
RHO_X <- 0.30
GRID_CONT <- seq(-0.75, 0.75, length.out = 5)
CI_ALPHA <- 0.05

## Bernstein sieve degrees, matching the PH--copula pointwise simulation.
M1 <- 6L
M2 <- 6L

## Simpson integration grid for computing the fitted win ratio surface.
## Must be odd.  201 is usually ample because all fitted functions are smooth.
N_INT <- 201L
if (N_INT %% 2L != 1L) stop("N_INT must be odd for Simpson integration.")

N_AVAILABLE <- parallel::detectCores(logical = TRUE)
if (!is.finite(N_AVAILABLE) || N_AVAILABLE < 2L) N_AVAILABLE <- 2L

if (RUN_MODE == "pilot") {
  ## Targeted diagnostic requested after the first PH--copula runs.
  ## Only S1/S3 are run here because the purpose is to check whether
  ## removing unnecessary quadratic PH terms improves finite-sample
  ## multiplier-bootstrap inference in the Gumbel scenarios.
  R_MC <- 50L
  B_BOOT <- 200L
  SCENARIOS_TO_RUN <- c("S1", "S3")
  N_CORES <- max(1L, min(16L, N_AVAILABLE - 1L))
  BATCH_SIZE <- 1L
  RESUME <- FALSE
} else {
  R_MC <- 500L
  B_BOOT <- 200L
  SCENARIOS_TO_RUN <- c("S1", "S2", "S3", "S4")
  N_CORES <- max(1L, min(16L, N_AVAILABLE - 1L))
  BATCH_SIZE <- 10L
  RESUME <- TRUE
}

## Original fits use multistart; bootstrap fits use a warm start.
ORIG_MAXIT <- 1200L
BOOT_MAXIT <- 1200L
OPT_RELTOL <- 1e-7

DESIGN_VERSION <- paste0("phcopula_coeff_scenario_basis_final_v3_", RUN_MODE)
OUTPUT_PREFIX <- paste0("coefficient_phcopula_table2_FINAL_revised_", RUN_MODE)

BASE_SEEDS <- c(
  S1 = 20261021L,
  S2 = 20261022L,
  S3 = 20261023L,
  S4 = 20261024L
)

`%||%` <- function(a, b) if (!is.null(a)) a else b

simulate_covariates_2d <- function(n, rho_x = 0.30) {
  Sigma <- matrix(c(1, rho_x, rho_x, 1), 2, 2)
  MASS::mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
}

############################################################
## SCENARIO CONFIGURATION
############################################################

scenario_config <- function(id) {
  id <- as.character(id)

  if (id == "S1") {
    return(list(
      id = "S1",
      label = "Constant WR; PWF holds",
      covariate_type = "two_continuous",
      x_grid = as.matrix(expand.grid(x1 = GRID_CONT, x2 = GRID_CONT)),
      theta_true = c(theta0 = log(1.5), theta1 = 0, theta2 = 0)
    ))
  }

  if (id == "S2") {
    return(list(
      id = "S2",
      label = "Approximately constant WR; PWF fails",
      covariate_type = "two_continuous",
      x_grid = as.matrix(expand.grid(x1 = GRID_CONT, x2 = GRID_CONT)),
      theta_true = c(
        theta0 = 0.4071280764,
        theta1 = -0.0392751643,
        theta2 = 0.0265045469
      )
    ))
  }

  if (id == "S3") {
    return(list(
      id = "S3",
      label = "Log-linear WR; PWF holds",
      covariate_type = "two_continuous",
      x_grid = as.matrix(expand.grid(x1 = GRID_CONT, x2 = GRID_CONT)),
      theta_true = c(theta0 = log(1.5), theta1 = 0.5, theta2 = 0.25)
    ))
  }

  if (id == "S4") {
    return(list(
      id = "S4",
      label = "Approximately log-linear WR; PWF fails (Candidate 6)",
      covariate_type = "two_continuous",
      x_grid = as.matrix(expand.grid(x1 = GRID_CONT, x2 = GRID_CONT)),
      theta_true = c(
        theta0 = 0.496292782,
        theta1 = 0.380788028,
        theta2 = 0.260591222
      )
    ))
  }

  stop("Unknown scenario id: ", id)
}

############################################################
## PWF-FAILURE DGPs
############################################################

make_constant_nonproportional_2d <- function() {
  arm1 <- list(
    eta_death_fun = function(x1, x2) {
      -0.44 +
        0.08 * x1 -
        0.06 * x2 +
        0.10 * x1^2 -
        0.08 * x1 * x2 -
        0.04 * x2^2
    },
    eta_prog_fun = function(x1, x2) {
      -0.17 -
        0.09 * x1 +
        0.07 * x2 -
        0.13 * x1^2 +
        0.10 * x1 * x2 +
        0.05 * x2^2
    },
    lambda_death = 0.20,
    lambda_prog = 0.22,
    theta = 1.2
  )

  arm0 <- list(
    eta_death_fun = function(x1, x2) {
      0.03 * x1 - 0.025 * x2
    },
    eta_prog_fun = function(x1, x2) {
      -0.03 * x1 + 0.025 * x2
    },
    lambda_death = 0.20,
    lambda_prog = 0.22,
    theta = 1.2
  )

  list(arm1 = arm1, arm0 = arm0)
}

make_candidate6 <- function() {
  arm1 <- list(
    eta_death_fun = function(x1, x2) {
      -0.55 - 0.45 * x1 - 0.30 * x2
    },
    eta_prog_fun = function(x1, x2) {
      -0.20 +
        0.05 * x1 +
        0.05 * x2 +
        0.45 * x1^2 -
        0.35 * x1 * x2
    },
    lambda_death = 0.20,
    lambda_prog = 0.22,
    theta = 1.2
  )

  arm0 <- list(
    eta_death_fun = function(x1, x2) {
      0.00
    },
    eta_prog_fun = function(x1, x2) {
      0.05 * x1 + 0.05 * x2
    },
    lambda_death = 0.20,
    lambda_prog = 0.22,
    theta = 1.2
  )

  list(arm1 = arm1, arm0 = arm0)
}

############################################################
## TERMINAL-EVENT OBSERVATION
############################################################

make_observed_terminal <- function(x1, x2, T1, T2, C, A) {
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
## SIMULATE PWF-HOLDS DATA
##
## Same Gumbel-Hougaard / Lehmann construction as the old
## coefficient simulation, now observed through the final
## terminal-event scheme.
############################################################

simulate_pwf_holds <- function(
    id, n0, n1, seed, tmax, censor_rate, rho_x = 0.30) {

  set.seed(seed)

  if (!id %in% c("S1", "S3")) {
    stop("PWF-holds simulator only supports S1/S3.")
  }

  n <- n0 + n1
  X <- simulate_covariates_2d(n, rho_x = rho_x)
  x1 <- X[, 1]
  x2 <- X[, 2]

  ## Fixed arm sizes with randomized treatment assignment.
  A <- sample(c(rep(0L, n0), rep(1L, n1)))

  beta1 <- 0.3
  beta2 <- -0.3
  betaA <- log(1.5)

  if (id == "S1") {
    gamma1 <- 0
    gamma2 <- 0
  } else {
    gamma1 <- 0.5
    gamma2 <- 0.25
  }

  eta <- beta1 * x1 + beta2 * x2 +
    betaA * A + gamma1 * A * x1 + gamma2 * A * x2

  cop <- copula::gumbelCopula(param = 2, dim = 2)
  UV <- copula::rCopula(n, cop)

  lambda1 <- 0.20
  lambda2 <- 2.00

  T1 <- -log(UV[, 1]) * exp(eta) / lambda1
  T2 <- -log(UV[, 2]) * exp(eta) / lambda2

  C <- pmin(rexp(n, rate = censor_rate), tmax)

  dat <- make_observed_terminal(x1, x2, T1, T2, C, A)

  list(
    dat = dat,
    dat0 = dat[dat$A == 0L, , drop = FALSE],
    dat1 = dat[dat$A == 1L, , drop = FALSE]
  )
}

############################################################
## SIMULATE PWF-FAILS DATA
############################################################

simulate_one_fail_arm <- function(n, arm, A, seed, tmax, censor_rate, rho_x = 0.3) {
  set.seed(seed)

  X <- simulate_covariates_2d(n, rho_x = rho_x)

  cop <- copula::claytonCopula(arm$theta, dim = 2)
  UV <- copula::rCopula(n, cop)

  eta1 <- arm$eta_death_fun(X[, 1], X[, 2])
  eta2 <- arm$eta_prog_fun(X[, 1], X[, 2])

  T1 <- qexp(UV[, 1], rate = arm$lambda_death * exp(eta1))
  T2 <- qexp(UV[, 2], rate = arm$lambda_prog * exp(eta2))

  C <- pmin(rexp(n, rate = censor_rate), tmax)

  make_observed_terminal(
    x1 = X[, 1],
    x2 = X[, 2],
    T1 = T1,
    T2 = T2,
    C = C,
    A = rep.int(A, n)
  )
}

simulate_pwf_fails <- function(
    id, n0, n1, seed, tmax, censor_rate, rho_x = 0.30) {
  if (id == "S2") {
    dgp <- make_constant_nonproportional_2d()
  } else if (id == "S4") {
    dgp <- make_candidate6()
  } else {
    stop("PWF-fails simulator only supports S2/S4.")
  }

  dat1 <- simulate_one_fail_arm(
    n = n1,
    arm = dgp$arm1,
    A = 1L,
    seed = seed + 10000L,
    tmax = tmax,
    censor_rate = censor_rate,
    rho_x = rho_x
  )

  dat0 <- simulate_one_fail_arm(
    n = n0,
    arm = dgp$arm0,
    A = 0L,
    seed = seed + 20000L,
    tmax = tmax,
    censor_rate = censor_rate,
    rho_x = rho_x
  )

  list(
    dat = rbind(dat0, dat1),
    dat0 = dat0,
    dat1 = dat1
  )
}

simulate_scenario <- function(
    id, n0, n1, seed, tmax, censor_rate, rho_x = 0.30) {

  if (id %in% c("S1", "S3")) {
    simulate_pwf_holds(
      id, n0, n1, seed, tmax, censor_rate, rho_x = rho_x
    )
  } else {
    simulate_pwf_fails(
      id, n0, n1, seed, tmax, censor_rate, rho_x = rho_x
    )
  }
}

############################################################
## BEST LINEAR PROJECTION ON THE RAW 5 x 5 GRID
############################################################

grid_weights_2d <- function(x_grid) {
  rep(1 / nrow(x_grid), nrow(x_grid))
}

fit_weighted_plane <- function(y, x_grid, w = NULL) {
  if (is.null(w)) w <- grid_weights_2d(x_grid)
  if (length(y) != nrow(x_grid)) return(NULL)
  if (any(!is.finite(y))) return(NULL)

  D <- cbind(1, x_grid[, 1], x_grid[, 2])
  DtQD <- crossprod(D, w * D)
  if (qr(DtQD)$rank < 3L) return(NULL)

  beta <- solve(DtQD, crossprod(D, w * y))
  names(beta) <- c("theta0", "theta1", "theta2")
  as.numeric(beta)
}

############################################################
## PH MODEL BASIS
##
## S1/S3 have truly linear marginal PH predictors, so use only
## (x1, x2).  S2/S4 genuinely require quadratic/interaction terms.
## This keeps the PH--copula working model correctly specified while
## avoiding unnecessary nuisance parameters in S1/S3.
############################################################

basis_for_scenario <- function(id) {
  if (id %in% c("S1", "S3")) "linear" else "quadratic"
}

ph_basis_2d <- function(x1, x2, basis_type = c("linear", "quadratic")) {
  basis_type <- match.arg(basis_type)

  if (basis_type == "linear") {
    return(cbind(
      x1 = x1,
      x2 = x2
    ))
  }

  cbind(
    x1 = x1,
    x2 = x2,
    x1_sq = x1^2,
    x1_x2 = x1 * x2,
    x2_sq = x2^2
  )
}

prepare_ph_arm <- function(
    dat,
    scenario_id,
    tmax = TMAX,
    m1 = M1,
    m2 = M2) {

  out <- dat
  basis_type <- basis_for_scenario(scenario_id)
  X <- ph_basis_2d(out$x1, out$x2, basis_type = basis_type)
  storage.mode(X) <- "double"

  attr(out, "X_mat") <- X
  attr(out, "basis_type") <- basis_type
  attr(out, "scenario_id") <- scenario_id
  attr(out, "tmax") <- tmax

  ## Precompute Bernstein bases at observed times.  These do not change
  ## across multiplier-bootstrap refits and substantially reduce runtime.
  attr(out, "B1_m") <- bernstein_basis(out$U1, m = m1, lower = 0, upper = tmax)
  attr(out, "B1_m1") <- bernstein_basis(out$U1, m = m1 - 1L, lower = 0, upper = tmax)
  attr(out, "B2_m") <- bernstein_basis(out$U2, m = m2, lower = 0, upper = tmax)
  attr(out, "B2_m1") <- bernstein_basis(out$U2, m = m2 - 1L, lower = 0, upper = tmax)
  out
}

ph_grid_basis <- function(x_grid, basis_type) {
  ph_basis_2d(
    x_grid[, 1],
    x_grid[, 2],
    basis_type = basis_type
  )
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
## COPULA COMPONENTS
############################################################

## Gumbel used on the SURVIVAL scale for S1/S3.
gumbel_components <- function(u, v, theta) {
  if (!is.finite(theta) || theta <= 1) stop("Gumbel theta must exceed 1.")
  eps <- 1e-12
  u <- pmin(pmax(u, eps), 1 - eps)
  v <- pmin(pmax(v, eps), 1 - eps)
  a <- -log(u)
  b <- -log(v)
  A <- pmax(a^theta + b^theta, eps)
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

## Clayton formulas used on the CDF scale for S2/S4.
## The implied survival copula is the 180-degree rotated Clayton copula.
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
## PARAMETERIZATION
############################################################

theta_from_zeta_gumbel <- function(zeta) 1 + exp(zeta)
zeta_from_theta_gumbel <- function(theta) {
  if (theta <= 1) stop("Gumbel theta must exceed 1.")
  log(theta - 1)
}

theta_from_zeta_clayton <- function(zeta) exp(zeta)
zeta_from_theta_clayton <- function(theta) {
  if (theta <= 0) stop("Clayton theta must be positive.")
  log(theta)
}

unpack_ph_par <- function(par, p, m1, m2, family) {
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

  theta <- if (family == "gumbel_survival") {
    theta_from_zeta_gumbel(zeta)
  } else {
    theta_from_zeta_clayton(zeta)
  }

  list(
    beta1 = beta1,
    beta2 = beta2,
    theta = theta,
    phi1 = phi_from_eta(eta1),
    phi2 = phi_from_eta(eta2)
  )
}

make_ph_start <- function(dat, m1, m2, family, theta0) {
  X <- attr(dat, "X_mat")
  p <- ncol(X)
  upper <- attr(dat, "tmax")
  lam1 <- max(sum(dat$Delta1) / max(sum(dat$U1), 1e-8), 1e-3)
  lam2 <- max(sum(dat$Delta2) / max(sum(dat$U2), 1e-8), 1e-3)
  inc1 <- rep(lam1 * upper / m1, m1)
  inc2 <- rep(lam2 * upper / m2, m2)

  zeta <- if (family == "gumbel_survival") {
    zeta_from_theta_gumbel(theta0)
  } else {
    zeta_from_theta_clayton(theta0)
  }

  c(rep(0, p), rep(0, p), zeta, log(inc1), log(inc2))
}

############################################################
## INDIVIDUAL LOG-LIKELIHOODS
############################################################

individual_loglik_gumbel_survival <- function(par, dat, m1 = M1, m2 = M2) {
  X <- attr(dat, "X_mat")
  p <- ncol(X)
  upper <- attr(dat, "tmax")
  pars <- unpack_ph_par(par, p, m1, m2, family = "gumbel_survival")

  lp1 <- as.vector(X %*% pars$beta1)
  lp2 <- as.vector(X %*% pars$beta2)
  b1 <- baseline_eval_precomputed(dat, pars$phi1, margin = 1L, upper = upper)
  b2 <- baseline_eval_precomputed(dat, pars$phi2, margin = 2L, upper = upper)

  logS1 <- -exp(lp1) * b1$Lambda
  logS2 <- -exp(lp2) * b2$Lambda
  S1 <- exp(pmax(logS1, -745))
  S2 <- exp(pmax(logS2, -745))
  logf1 <- log(b1$lambda) + lp1 + logS1
  logf2 <- log(b2$lambda) + lp2 + logS2
  cc <- gumbel_components(S1, S2, pars$theta)

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

individual_loglik_clayton_cdf <- function(par, dat, m1 = M1, m2 = M2) {
  X <- attr(dat, "X_mat")
  p <- ncol(X)
  upper <- attr(dat, "tmax")
  pars <- unpack_ph_par(par, p, m1, m2, family = "clayton_cdf")

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

  ## J(t1,t2) = P(T1>t1,T2>t2)
  J <- S1 + S2 - 1 + cc$C
  J <- pmax(J, eps)

  ll <- numeric(nrow(dat))
  ll[both_events] <- logf1[both_events] + logf2[both_events] + log(cc$density[both_events])
  ll[only1] <- logf1[only1] + log(pmax(1 - cc$Cu[only1], eps))
  ll[only2] <- logf2[only2] + log(pmax(1 - cc$Cv[only2], eps))
  ll[both_cens] <- log(J[both_cens])
  ll
}

############################################################
## PH--COPULA FITTING
############################################################

family_for_scenario <- function(id) {
  if (id %in% c("S1", "S3")) "gumbel_survival" else "clayton_cdf"
}

fit_phcopula_arm <- function(
    dat,
    family,
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
  ## Scaling all weights by a common constant does not change the maximizer.
  weights <- weights / mean(weights)

  llfun <- if (family == "gumbel_survival") {
    individual_loglik_gumbel_survival
  } else if (family == "clayton_cdf") {
    individual_loglik_clayton_cdf
  } else {
    stop("Unknown PH--copula family.")
  }

  objective <- function(par) {
    ll <- llfun(par, dat, m1 = M1, m2 = M2)
    if (any(!is.finite(ll))) return(1e100)
    val <- sum(weights * ll)
    if (!is.finite(val)) return(1e100)
    -val
  }

  if (!is.null(start)) {
    starts <- list(start)
  } else if (family == "gumbel_survival") {
    starts <- lapply(c(1.4, 2.0, 3.5), function(th) {
      make_ph_start(dat, M1, M2, family, theta0 = th)
    })
  } else {
    ## The PWF-failure DGP uses theta = 1.2; use several dispersed starts.
    starts <- lapply(c(0.5, 1.2, 2.0), function(th) {
      make_ph_start(dat, M1, M2, family, theta0 = th)
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

  ## One fallback if needed.
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
  pars <- unpack_ph_par(best$par, p, M1, M2, family)

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
    family = family,
    basis_type = attr(dat, "basis_type"),
    scenario_id = attr(dat, "scenario_id")
  )
}

############################################################
## FAST VECTORIZED WIN-RATIO EVALUATION ON THE 5 x 5 GRID
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
  ## Rows = time, columns = covariate grid points.
  cumhaz <- outer(b$Lambda, rr, "*")
  S <- exp(-cumhaz)
  f <- outer(b$lambda, rr, "*") * S
  list(S = S, f = f)
}

copula_matrix_gumbel <- function(U, V, theta, component = c("C", "Cv")) {
  component <- match.arg(component)
  d <- dim(U)
  cc <- gumbel_components(U, V, theta)
  matrix(cc[[component]], nrow = d[1], ncol = d[2])
}

copula_matrix_clayton_cdf <- function(U, V, theta, component = c("C", "Cv")) {
  component <- match.arg(component)
  d <- dim(U)
  cc <- clayton_cdf_components(U, V, theta)
  matrix(cc[[component]], nrow = d[1], ncol = d[2])
}

logwr_surface_from_fits <- function(fit1, fit0, raw_grid, tau = TAU, n_int = N_INT) {
  family <- fit1$family
  if (!identical(family, fit0$family)) {
    stop("Treatment arms use different copula families.")
  }
  if (!identical(fit1$basis_type, fit0$basis_type)) {
    stop("Treatment arms use different PH covariate bases.")
  }

  Zg <- ph_grid_basis(raw_grid, basis_type = fit1$basis_type)
  rule <- simpson_rule(tau, n = n_int)
  tt <- rule$t
  ww <- rule$w

  m11 <- margin_grid(fit1, 1L, Zg, tt)
  m10 <- margin_grid(fit0, 1L, Zg, tt)
  m21 <- margin_grid(fit1, 2L, Zg, tt)
  m20 <- margin_grid(fit0, 2L, Zg, tt)

  Wp <- colSums((m11$S * m10$f) * ww)
  Lp <- colSums((m10$S * m11$f) * ww)

  S11_tau <- m11$S[nrow(m11$S), ]
  S10_tau <- m10$S[nrow(m10$S), ]
  S11_mat <- matrix(S11_tau, nrow = length(tt), ncol = length(S11_tau), byrow = TRUE)
  S10_mat <- matrix(S10_tau, nrow = length(tt), ncol = length(S10_tau), byrow = TRUE)

  if (family == "gumbel_survival") {
    J1 <- copula_matrix_gumbel(S11_mat, m21$S, fit1$theta, "C")
    J0 <- copula_matrix_gumbel(S10_mat, m20$S, fit0$theta, "C")
    Cv1 <- copula_matrix_gumbel(S11_mat, m21$S, fit1$theta, "Cv")
    Cv0 <- copula_matrix_gumbel(S10_mat, m20$S, fit0$theta, "Cv")
    minus_dJ1 <- Cv1 * m21$f
    minus_dJ0 <- Cv0 * m20$f
  } else if (family == "clayton_cdf") {
    F11_tau <- 1 - S11_mat
    F10_tau <- 1 - S10_mat
    F21 <- 1 - m21$S
    F20 <- 1 - m20$S

    C1 <- copula_matrix_clayton_cdf(F11_tau, F21, fit1$theta, "C")
    C0 <- copula_matrix_clayton_cdf(F10_tau, F20, fit0$theta, "C")
    Cv1 <- copula_matrix_clayton_cdf(F11_tau, F21, fit1$theta, "Cv")
    Cv0 <- copula_matrix_clayton_cdf(F10_tau, F20, fit0$theta, "Cv")

    J1 <- S11_mat + m21$S - 1 + C1
    J0 <- S10_mat + m20$S - 1 + C0
    J1 <- pmax(J1, 0)
    J0 <- pmax(J0, 0)

    ## -dJ/dt = f2(t) * {1 - C_v(F1,F2)}.
    minus_dJ1 <- (1 - Cv1) * m21$f
    minus_dJ0 <- (1 - Cv0) * m20$f
  } else {
    stop("Unknown copula family in fitted object.")
  }

  Ws <- colSums((J1 * minus_dJ0) * ww)
  Ls <- colSums((J0 * minus_dJ1) * ww)

  W <- Wp + Ws
  L <- Lp + Ls
  bad <- !is.finite(W) | !is.finite(L) | W <= 0 | L <= 0
  out <- log(W / L)
  out[bad] <- NA_real_
  out
}

ph_projection_from_fits <- function(fit1, fit0, raw_grid, tau = TAU) {
  logwr <- logwr_surface_from_fits(fit1, fit0, raw_grid, tau = tau)
  if (any(!is.finite(logwr))) return(NULL)
  fit_weighted_plane(logwr, raw_grid)
}

############################################################
## ONE ORIGINAL PH--COPULA FIT + MULTIPLIER BOOTSTRAP
############################################################

fit_phcopula_coefficients <- function(
    dat1,
    dat0,
    raw_grid,
    scenario_id,
    tau,
    B,
    seed) {

  family <- family_for_scenario(scenario_id)
  d1 <- prepare_ph_arm(
    dat1, scenario_id = scenario_id,
    tmax = TMAX, m1 = M1, m2 = M2
  )
  d0 <- prepare_ph_arm(
    dat0, scenario_id = scenario_id,
    tmax = TMAX, m1 = M1, m2 = M2
  )

  fit1 <- fit_phcopula_arm(
    d1, family = family, multistart = TRUE,
    maxit = ORIG_MAXIT, reltol = OPT_RELTOL
  )
  fit0 <- fit_phcopula_arm(
    d0, family = family, multistart = TRUE,
    maxit = ORIG_MAXIT, reltol = OPT_RELTOL
  )

  if (is.null(fit1) || is.null(fit0) ||
      fit1$convergence != 0L || fit0$convergence != 0L) {
    return(NULL)
  }

  theta_hat <- ph_projection_from_fits(fit1, fit0, raw_grid, tau = tau)
  if (is.null(theta_hat) || any(!is.finite(theta_hat))) return(NULL)

  set.seed(seed)
  theta_boot <- matrix(NA_real_, nrow = B, ncol = 3L)
  colnames(theta_boot) <- c("theta0", "theta1", "theta2")

  for (b in seq_len(B)) {
    w1 <- rexp(nrow(d1), rate = 1)
    w0 <- rexp(nrow(d0), rate = 1)

    fb1 <- try(
      fit_phcopula_arm(
        d1,
        family = family,
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
        family = family,
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

    thb <- try(
      ph_projection_from_fits(fb1, fb0, raw_grid, tau = tau),
      silent = TRUE
    )
    if (!inherits(thb, "try-error") && !is.null(thb) && all(is.finite(thb))) {
      theta_boot[b, ] <- thb
    }
  }

  ok <- complete.cases(theta_boot)
  n_success <- sum(ok)
  min_success <- if (B <= 20L) {
    max(3L, ceiling(0.5 * B))
  } else {
    max(20L, ceiling(0.5 * B))
  }
  if (n_success < min_success) return(NULL)

  theta_se <- apply(theta_boot[ok, , drop = FALSE], 2, sd)

  list(
    theta_hat = theta_hat,
    theta_se = theta_se,
    theta_boot = theta_boot,
    n_success = n_success,
    family = family,
    theta_copula_arm1 = fit1$theta,
    theta_copula_arm0 = fit0$theta
  )
}

############################################################
## ONE MONTE CARLO REPLICATION
############################################################

run_one_rep_ph <- function(r, scenario_id, settings) {
  cfg <- scenario_config(scenario_id)
  out <- list(
    r = r,
    ph_theta = rep(NA_real_, 3L),
    ph_se = rep(NA_real_, 3L),
    ph_boot_success = NA_integer_,
    copula_theta1 = NA_real_,
    copula_theta0 = NA_real_,
    error = NULL
  )

  datsim <- try(
    simulate_scenario(
      id = scenario_id,
      n0 = settings$n0,
      n1 = settings$n1,
      seed = settings$seed + r,
      tmax = settings$tmax,
      censor_rate = settings$censor_rate,
      rho_x = settings$rho_x
    ),
    silent = TRUE
  )

  if (inherits(datsim, "try-error")) {
    out$error <- paste("simulation:", as.character(datsim))
    return(out)
  }

  z <- try(
    fit_phcopula_coefficients(
      dat1 = datsim$dat1,
      dat0 = datsim$dat0,
      raw_grid = cfg$x_grid,
      scenario_id = scenario_id,
      tau = settings$tau,
      B = settings$B,
      seed = settings$seed + 600000L + r
    ),
    silent = TRUE
  )

  if (inherits(z, "try-error") || is.null(z)) {
    out$error <- "PH--copula fit/bootstrap failed"
    return(out)
  }

  out$ph_theta <- as.numeric(z$theta_hat)
  out$ph_se <- as.numeric(z$theta_se)
  out$ph_boot_success <- z$n_success
  out$copula_theta1 <- z$theta_copula_arm1
  out$copula_theta0 <- z$theta_copula_arm0
  out
}

############################################################
## SUMMARIES
############################################################

summarize_ph_method <- function(theta_mat, se_mat, theta_true, alpha = CI_ALPHA) {
  zcrit <- qnorm(1 - alpha / 2)
  good <- apply(theta_mat, 1, function(z) all(is.finite(z))) &
    apply(se_mat, 1, function(z) all(is.finite(z)))

  th <- theta_mat[good, , drop = FALSE]
  ss <- se_mat[good, , drop = FALSE]

  if (nrow(th) == 0L) {
    return(data.frame(
      parameter = c("intercept", "x1_slope", "x2_slope"),
      truth = as.numeric(theta_true),
      EST_PH = NA_real_, BIAS_PH = NA_real_, SE_PH = NA_real_,
      SEE_PH = NA_real_, CP_PH = NA_real_, n_valid = 0L
    ))
  }

  est <- colMeans(th)
  bias <- est - theta_true
  emp_se <- apply(th, 2, sd)
  mean_se <- colMeans(ss)
  cp <- vapply(seq_len(3L), function(j) {
    mean(
      th[, j] - zcrit * ss[, j] <= theta_true[j] &
      th[, j] + zcrit * ss[, j] >= theta_true[j]
    )
  }, numeric(1))

  data.frame(
    parameter = c("intercept", "x1_slope", "x2_slope"),
    truth = as.numeric(theta_true),
    EST_PH = as.numeric(est),
    BIAS_PH = as.numeric(bias),
    SE_PH = as.numeric(emp_se),
    SEE_PH = as.numeric(mean_se),
    CP_PH = as.numeric(cp),
    n_valid = nrow(th)
  )
}

summarize_scenario_ph <- function(results, scenario_id) {
  cfg <- scenario_config(scenario_id)
  R <- length(results)
  theta <- matrix(NA_real_, R, 3L)
  se <- matrix(NA_real_, R, 3L)
  boot_success <- rep(NA_real_, R)
  ctheta1 <- rep(NA_real_, R)
  ctheta0 <- rep(NA_real_, R)
  errors <- character(0)

  for (r in seq_len(R)) {
    z <- results[[r]]
    if (is.null(z)) next
    theta[r, ] <- z$ph_theta
    se[r, ] <- z$ph_se
    boot_success[r] <- z$ph_boot_success
    ctheta1[r] <- z$copula_theta1
    ctheta0[r] <- z$copula_theta0
    if (!is.null(z$error)) errors <- c(errors, z$error)
  }

  colnames(theta) <- colnames(se) <- c("theta0", "theta1", "theta2")
  tab <- summarize_ph_method(theta, se, cfg$theta_true)
  tab$scenario <- cfg$label
  tab$family <- family_for_scenario(scenario_id)
  tab <- tab[, c(
    "scenario", "family", "parameter", "truth",
    "EST_PH", "BIAS_PH", "SE_PH", "SEE_PH", "CP_PH", "n_valid"
  )]

  list(
    scenario = cfg,
    table = tab,
    raw = list(
      ph_theta = theta,
      ph_se = se,
      ph_boot_success = boot_success,
      copula_theta1 = ctheta1,
      copula_theta0 = ctheta0,
      errors = errors
    ),
    diagnostics = list(
      n_valid_ph = sum(complete.cases(theta, se)),
      mean_boot_success = mean(boot_success, na.rm = TRUE),
      min_boot_success = if (any(is.finite(boot_success))) min(boot_success[is.finite(boot_success)]) else NA_real_,
      mean_copula_theta_arm1 = mean(ctheta1, na.rm = TRUE),
      mean_copula_theta_arm0 = mean(ctheta0, na.rm = TRUE),
      error_count = length(errors)
    )
  )
}

############################################################
## PARALLEL RUNNER WITH CHECKPOINTING
############################################################

run_scenario_parallel_ph <- function(
    scenario_id,
    settings,
    n_cores,
    batch_size,
    resume = TRUE) {

  cfg <- scenario_config(scenario_id)
  checkpoint_file <- paste0(OUTPUT_PREFIX, "_", scenario_id, "_checkpoint.rds")
  results <- vector("list", settings$R)

  if (resume && file.exists(checkpoint_file)) {
    old <- try(readRDS(checkpoint_file), silent = TRUE)
    if (!inherits(old, "try-error") && is.list(old) &&
        identical(old$scenario_id, scenario_id) &&
        identical(old$design_version, DESIGN_VERSION) &&
        identical(old$settings, settings) &&
        length(old$results) == settings$R) {
      results <- old$results
      cat("Resuming", scenario_id, "from checkpoint.\n")
    }
  }

  pending <- which(vapply(results, is.null, logical(1)))
  if (length(pending) == 0L) {
    cat(scenario_id, "already complete.\n")
    return(summarize_scenario_ph(results, scenario_id))
  }

  cat("\n============================================================\n")
  cat(cfg$label, "\n")
  cat("PH--copula only | family =", family_for_scenario(scenario_id),
      "| R =", settings$R, "| B =", settings$B,
      "| n/arm =", settings$n0, "| tau =", settings$tau,
      "| cores =", n_cores, "\n")
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
    "GRID_CONT", "TMAX", "TAU", "M1", "M2", "N_INT",
    "ORIG_MAXIT", "BOOT_MAXIT", "OPT_RELTOL", "CI_ALPHA",
    "simulate_covariates_2d", "scenario_config",
    "make_constant_nonproportional_2d", "make_candidate6",
    "make_observed_terminal", "simulate_pwf_holds",
    "simulate_one_fail_arm", "simulate_pwf_fails", "simulate_scenario",
    "grid_weights_2d", "fit_weighted_plane",
    "basis_for_scenario", "ph_basis_2d", "prepare_ph_arm", "ph_grid_basis",
    "bernstein_basis", "phi_from_eta", "baseline_eval", "baseline_eval_precomputed",
    "gumbel_components", "clayton_cdf_components",
    "theta_from_zeta_gumbel", "zeta_from_theta_gumbel",
    "theta_from_zeta_clayton", "zeta_from_theta_clayton",
    "unpack_ph_par", "make_ph_start",
    "individual_loglik_gumbel_survival", "individual_loglik_clayton_cdf",
    "family_for_scenario", "fit_phcopula_arm",
    "simpson_rule", "margin_grid", "copula_matrix_gumbel",
    "copula_matrix_clayton_cdf", "logwr_surface_from_fits",
    "ph_projection_from_fits", "fit_phcopula_coefficients",
    "run_one_rep_ph"
  )

  parallel::clusterExport(cl, export_names, envir = .GlobalEnv)

  batches <- split(pending, ceiling(seq_along(pending) / batch_size))
  for (idx in batches) {
    ans <- parallel::parLapply(
      cl,
      idx,
      run_one_rep_ph,
      scenario_id = scenario_id,
      settings = settings
    )

    for (k in seq_along(idx)) results[[idx[k]]] <- ans[[k]]

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
    cat("Completed", completed, "/", settings$R, "for", scenario_id, "\n")
    flush.console()
  }

  parallel::stopCluster(cl)
  on.exit(NULL, add = FALSE)
  summarize_scenario_ph(results, scenario_id)
}

############################################################
## RUN ALL FOUR SCENARIOS
############################################################

cat("\n============================================================\n")
cat("PH--COPULA-ONLY COEFFICIENT SIMULATION FOR TABLE 2\n")
cat("============================================================\n")
cat("Run mode:", RUN_MODE, "\n")
cat("Common raw grid:", paste(round(GRID_CONT, 3), collapse = ", "), "x itself\n")
cat("PH covariate basis: S1/S3 = x1, x2; S2/S4 = x1, x2, x1^2, x1*x2, x2^2\n")
cat("Scenarios in this run:", paste(SCENARIOS_TO_RUN, collapse = ", "), "\n")
cat("S1/S3: Gumbel survival copula\n")
cat("S2/S4: Clayton CDF copula = rotated-Clayton survival copula\n")
cat("R =", R_MC, "| B =", B_BOOT, "| n/arm =", N0,
    "| tau =", TAU, "| cores =", N_CORES, "\n")
cat("============================================================\n\n")

settings <- list(
  R = R_MC,
  B = B_BOOT,
  n0 = N0,
  n1 = N1,
  tau = TAU,
  tmax = TMAX,
  censor_rate = CENSOR_RATE,
  rho_x = RHO_X,
  m1 = M1,
  m2 = M2,
  n_int = N_INT,
  seed = NA_integer_
)

all_results <- vector("list", length(SCENARIOS_TO_RUN))
names(all_results) <- SCENARIOS_TO_RUN

for (sid in names(all_results)) {
  settings$seed <- unname(BASE_SEEDS[sid])
  all_results[[sid]] <- run_scenario_parallel_ph(
    scenario_id = sid,
    settings = settings,
    n_cores = N_CORES,
    batch_size = BATCH_SIZE,
    resume = RESUME
  )
  saveRDS(all_results, paste0(OUTPUT_PREFIX, "_partial.rds"))
}

############################################################
## FINAL TABLE AND DIAGNOSTICS
############################################################

final_ph_table <- do.call(rbind, lapply(all_results, function(z) z$table))
rownames(final_ph_table) <- NULL

cat("\n\n============================================================\n")
cat("PH--COPULA COLUMNS FOR TABLE 2\n")
cat("============================================================\n")
print_table <- final_ph_table
num <- vapply(print_table, is.numeric, logical(1))
print_table[num] <- lapply(print_table[num], function(z) round(z, 4))
print(print_table, row.names = FALSE)

cat("\n\n============================================================\n")
cat("DIAGNOSTICS\n")
cat("============================================================\n")
for (sid in names(all_results)) {
  cat("\n", sid, ": ", all_results[[sid]]$scenario$label, "\n", sep = "")
  print(all_results[[sid]]$diagnostics)
}

############################################################
## SAVE
############################################################

final_output <- list(
  settings = list(
    run_mode = RUN_MODE,
    R = R_MC,
    B = B_BOOT,
    n0 = N0,
    n1 = N1,
    tau = TAU,
    tmax = TMAX,
    censor_rate = CENSOR_RATE,
    rho_x = RHO_X,
    grid_cont = GRID_CONT,
    m1 = M1,
    m2 = M2,
    n_int = N_INT,
    ph_basis_S1_S3 = "x1 + x2",
    ph_basis_S2_S4 = "x1 + x2 + x1^2 + x1*x2 + x2^2",
    copula_S1_S3 = "Gumbel survival",
    copula_S2_S4 = "Clayton CDF / rotated Clayton survival",
    design_version = DESIGN_VERSION
  ),
  scenarios = lapply(SCENARIOS_TO_RUN, scenario_config),
  final_ph_table = final_ph_table,
  results = all_results
)

rds_file <- paste0(OUTPUT_PREFIX, ".rds")
csv_file <- paste0(OUTPUT_PREFIX, "_table.csv")

saveRDS(final_output, rds_file)
write.csv(final_ph_table, csv_file, row.names = FALSE)

cat("\nSaved:\n")
cat("  ", rds_file, "\n", sep = "")
cat("  ", csv_file, "\n", sep = "")
cat("\nTotal elapsed time:\n")
print(Sys.time() - time0)

############################################################
## FINAL LOCAL-IPCW COEFFICIENT SIMULATION WITH n^(-1/5)
##
## PURPOSE
## -------
## Rerun ONLY the Local IPCW portion of the four final
## coefficient-estimation scenarios using the inference-oriented
## two-dimensional bandwidth
##
##     h_j = sd(X_j) * n^(-1/5),  j = 1,2,
##
## where n is the pooled sample size used for smoothing.
##
## The bandwidth is computed ONCE from the original pooled sample
## within each Monte Carlo replication and is held fixed throughout
## all B bootstrap resamples.
##
## Four scenarios:
##   S1. Constant WR; PWF holds
##   S2. Approximately constant fixed-horizon WR; PWF fails
##   S3. Log-linear WR; PWF holds
##   S4. Approximately log-linear fixed-horizon WR; PWF fails
##
## Final design:
##   - all four scenarios use two continuous covariates
##   - Corr(X1, X2) = 0.30
##   - common 5 x 5 evaluation grid on [-0.75, 0.75]^2
##   - n0 = n1 = 250
##   - tau = 8
##   - tmax = 10
##   - R = 500
##   - B = 200
##   - exact event-time Local IPCW estimator
##   - normal bootstrap-SE confidence intervals
##   - PSOCK parallelization
##   - checkpointing every batch
##
## IMPORTANT
## ---------
## This program does NOT rerun Mao-Wang.
## Keep the Mao-Wang columns from the previous finalized
## all-continuous coefficient simulation unchanged, and replace
## only the Local IPCW columns by the output from this program.
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

R_MC <- 500L
B_BOOT <- 200L

N0 <- 250L
N1 <- 250L

TAU <- 8
TMAX <- 10
CENSOR_RATE <- 0.05
RHO_X <- 0.30

GRID_CONT <- seq(-0.75, 0.75, length.out = 5)
CI_ALPHA <- 0.05

N_AVAILABLE <- parallel::detectCores(logical = TRUE)
if (!is.finite(N_AVAILABLE) || N_AVAILABLE < 2L) N_AVAILABLE <- 2L

N_CORES <- max(1L, min(16L, N_AVAILABLE - 1L))
BATCH_SIZE <- 50L
RESUME <- TRUE

DESIGN_VERSION <- "all_continuous_local_n15_v3"
OUTPUT_PREFIX <- "coefficient_local_n15_all4_final"

BASE_SEEDS <- c(
  S1 = 20261021L,
  S2 = 20261022L,
  S3 = 20261023L,
  S4 = 20261024L
)

############################################################
## COVARIATE HELPER
############################################################

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
## LOCAL IPCW: KERNEL AND CENSORING SURVIVAL
############################################################

kernel_weights_2d <- function(X, x0, h, eps = 1e-12) {
  u1 <- (X[, 1] - x0[1]) / h[1]
  u2 <- (X[, 2] - x0[2]) / h[2]

  w <- dnorm(u1) * dnorm(u2)
  sw <- sum(w)

  if (!is.finite(sw) || sw < eps) return(rep(1 / nrow(X), nrow(X)))
  w / sw
}


weighted_km_surv <- function(time, status, weights, times, eps = 1e-6) {
  if (length(times) == 0L) return(numeric(0))

  ord <- order(time)
  time <- time[ord]
  status <- status[ord]
  weights <- weights[ord]

  event_times <- sort(unique(time[status == 1L]))

  if (length(event_times) == 0L) return(rep(1, length(times)))

  surv_vals <- numeric(length(event_times))
  S <- 1

  for (j in seq_along(event_times)) {
    tj <- event_times[j]
    risk <- sum(weights[time >= tj])
    event <- sum(weights[time == tj & status == 1L])

    if (is.finite(risk) && risk > eps) {
      S <- S * (1 - event / risk)
    }

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
    status = 1L - dat$Delta1,
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

############################################################
## LOCAL IPCW: EXACT EVENT-TIME STIELTJES ESTIMATOR
############################################################

wr_local_exact_2d <- function(x0, dat1, dat0, h, tau) {
  X1 <- as.matrix(dat1[, c("x1", "x2")])
  X0 <- as.matrix(dat0[, c("x1", "x2")])

  w1 <- kernel_weights_2d(X1, x0, h)
  w0 <- kernel_weights_2d(X0, x0, h)

  ##########################################################
  ## Primary component
  ##########################################################

  t0 <- sort(unique(dat0$U1[dat0$Delta1 == 1L & dat0$U1 <= tau]))

  if (length(t0) > 0L) {
    G0 <- Ghat_common(t0, dat0, w0)
    dF0 <- vapply(t0, function(tt) {
      sum(w0[dat0$Delta1 == 1L & dat0$U1 == tt])
    }, numeric(1)) / G0

    Wp <- sum(S1_local_raw(t0, dat1, w1) * dF0)
  } else {
    Wp <- 0
  }

  t1 <- sort(unique(dat1$U1[dat1$Delta1 == 1L & dat1$U1 <= tau]))

  if (length(t1) > 0L) {
    G1 <- Ghat_common(t1, dat1, w1)
    dF1 <- vapply(t1, function(tt) {
      sum(w1[dat1$Delta1 == 1L & dat1$U1 == tt])
    }, numeric(1)) / G1

    Lp <- sum(S1_local_raw(t1, dat0, w0) * dF1)
  } else {
    Lp <- 0
  }

  ##########################################################
  ## Secondary component
  ##########################################################

  G1tau <- as.numeric(Ghat_common(tau, dat1, w1))
  G0tau <- as.numeric(Ghat_common(tau, dat0, w0))

  u0 <- sort(unique(
    dat0$U2[
      dat0$U1 > tau &
        dat0$Delta2 == 1L &
        dat0$U2 <= tau
    ]
  ))

  if (length(u0) > 0L) {
    dJ0 <- vapply(u0, function(uu) {
      sum(w0[
        dat0$U1 > tau &
          dat0$Delta2 == 1L &
          dat0$U2 == uu
      ])
    }, numeric(1)) / G0tau

    Ws <- sum(J_local_raw(tau, u0, dat1, w1) * dJ0)
  } else {
    Ws <- 0
  }

  u1 <- sort(unique(
    dat1$U2[
      dat1$U1 > tau &
        dat1$Delta2 == 1L &
        dat1$U2 <= tau
    ]
  ))

  if (length(u1) > 0L) {
    dJ1 <- vapply(u1, function(uu) {
      sum(w1[
        dat1$U1 > tau &
          dat1$Delta2 == 1L &
          dat1$U2 == uu
      ])
    }, numeric(1)) / G1tau

    Ls <- sum(J_local_raw(tau, u1, dat0, w0) * dJ1)
  } else {
    Ls <- 0
  }

  W <- Wp + Ws
  L <- Lp + Ls

  if (!is.finite(W) || !is.finite(L) || W <= 0 || L <= 0) return(NA_real_)
  W / L
}

estimate_local_surface_2d <- function(dat1, dat0, x_grid, h, tau) {
  apply(x_grid, 1, function(x0) {
    wr_local_exact_2d(
      x0 = as.numeric(x0),
      dat1 = dat1,
      dat0 = dat0,
      h = h,
      tau = tau
    )
  })
}

############################################################
## PROJECTION
############################################################

grid_weights_2d <- function(x_grid) {
  rep(1 / nrow(x_grid), nrow(x_grid))
}

fit_weighted_plane <- function(y, x_grid, w = NULL) {
  if (is.null(w)) w <- grid_weights_2d(x_grid)
  if (length(y) != nrow(x_grid)) return(NULL)
  if (any(!is.finite(y))) return(NULL)

  X <- cbind(1, x_grid[, 1], x_grid[, 2])
  XtWX <- crossprod(X, w * X)

  if (qr(XtWX)$rank < 3L) return(NULL)

  beta <- solve(XtWX, crossprod(X, w * y))
  fitted <- as.vector(X %*% beta)

  names(beta) <- c("theta0", "theta1", "theta2")

  list(
    beta = as.numeric(beta),
    fitted = fitted,
    resid = y - fitted
  )
}

fit_loglinear_from_surface <- function(WR_hat, x_grid) {
  if (any(!is.finite(WR_hat)) || any(WR_hat <= 0)) return(NULL)
  fit <- fit_weighted_plane(log(WR_hat), x_grid)
  if (is.null(fit)) return(NULL)
  theta <- fit$beta
  names(theta) <- c("theta0", "theta1", "theta2")
  theta
}


############################################################
## INFERENCE-ORIENTED TWO-DIMENSIONAL BANDWIDTH
##
## h_j = sd(X_j) * n^(-1/5)
############################################################

bandwidth_n15_2d <- function(X) {
  n <- nrow(X)

  h <- apply(X, 2, sd) * n^(-1 / 5)

  h[!is.finite(h) | h <= 0] <- n^(-1 / 5)

  h
}


############################################################
## LOCAL BOOTSTRAP FOR COEFFICIENT INFERENCE
############################################################

bootstrap_local_coeff_2d <- function(
    dat1,
    dat0,
    x_grid,
    h,
    tau,
    B,
    seed) {

  set.seed(seed)

  wr_hat <- estimate_local_surface_2d(dat1, dat0, x_grid, h, tau)
  theta_hat <- fit_loglinear_from_surface(wr_hat, x_grid)

  if (is.null(theta_hat)) return(NULL)

  theta_boot <- matrix(NA_real_, nrow = B, ncol = 3L)
  colnames(theta_boot) <- c("theta0", "theta1", "theta2")

  n1 <- nrow(dat1)
  n0 <- nrow(dat0)

  for (b in seq_len(B)) {
    i1 <- sample.int(n1, n1, replace = TRUE)
    i0 <- sample.int(n0, n0, replace = TRUE)

    d1 <- dat1[i1, , drop = FALSE]
    d0 <- dat0[i0, , drop = FALSE]

    wrb <- try(
      estimate_local_surface_2d(d1, d0, x_grid, h, tau),
      silent = TRUE
    )

    if (inherits(wrb, "try-error")) next

    thb <- fit_loglinear_from_surface(wrb, x_grid)
    if (!is.null(thb) && all(is.finite(thb))) {
      theta_boot[b, ] <- thb
    }
  }

  ok <- complete.cases(theta_boot)
  n_success <- sum(ok)

  min_success <- if (B <= 20L) {
    max(3L, ceiling(0.5 * B))
  } else {
    max(20L, floor(0.5 * B))
  }

  if (n_success < min_success) return(NULL)

  theta_se <- apply(theta_boot[ok, , drop = FALSE], 2, sd)

  list(
    theta_hat = theta_hat,
    theta_se = theta_se,
    n_success = n_success,
    boot_theta = theta_boot
  )
}


############################################################
## ONE MONTE CARLO REPLICATION
############################################################

run_one_rep_local_n15 <- function(r, scenario_id, settings) {
  cfg <- scenario_config(scenario_id)

  out <- list(
    r = r,
    loc_theta = rep(NA_real_, 3L),
    loc_se = rep(NA_real_, 3L),
    h = rep(NA_real_, 2L),
    local_boot_success = NA_integer_,
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

  Xall <- as.matrix(rbind(
    datsim$dat1[, c("x1", "x2")],
    datsim$dat0[, c("x1", "x2")]
  ))

  h <- bandwidth_n15_2d(Xall)

  if (any(!is.finite(h)) || any(h <= 0)) {
    out$error <- "invalid n^(-1/5) bandwidth"
    return(out)
  }

  out$h <- h

  loc <- try(
    bootstrap_local_coeff_2d(
      dat1 = datsim$dat1,
      dat0 = datsim$dat0,
      x_grid = cfg$x_grid,
      h = h,
      tau = settings$tau,
      B = settings$B,
      seed = settings$seed + 500000L + r
    ),
    silent = TRUE
  )

  if (inherits(loc, "try-error") || is.null(loc)) {
    out$error <- "local bootstrap failed"
    return(out)
  }

  out$loc_theta <- as.numeric(loc$theta_hat)
  out$loc_se <- as.numeric(loc$theta_se)
  out$local_boot_success <- loc$n_success

  out
}

############################################################
## SUMMARIES
############################################################

summarize_local_method <- function(theta_mat, se_mat, theta_true, alpha = 0.05) {
  zcrit <- qnorm(1 - alpha / 2)

  good <- apply(theta_mat, 1, function(z) all(is.finite(z))) &
    apply(se_mat, 1, function(z) all(is.finite(z)))

  th <- theta_mat[good, , drop = FALSE]
  ss <- se_mat[good, , drop = FALSE]

  if (nrow(th) == 0L) {
    return(data.frame(
      parameter = c("intercept", "x1_slope", "x2_slope"),
      truth = as.numeric(theta_true),
      EST_Local = NA_real_,
      BIAS_Local = NA_real_,
      SE_Local = NA_real_,
      SEE_Local = NA_real_,
      CP_Local = NA_real_,
      n_valid = 0L
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
    EST_Local = as.numeric(est),
    BIAS_Local = as.numeric(bias),
    SE_Local = as.numeric(emp_se),
    SEE_Local = as.numeric(mean_se),
    CP_Local = as.numeric(cp),
    n_valid = nrow(th)
  )
}

summarize_scenario_local <- function(results, scenario_id, alpha = 0.05) {
  cfg <- scenario_config(scenario_id)
  R <- length(results)

  loc_theta <- matrix(NA_real_, R, 3L)
  loc_se <- matrix(NA_real_, R, 3L)
  h <- matrix(NA_real_, R, 2L)
  boot_success <- rep(NA_real_, R)
  errors <- character(0)

  for (r in seq_len(R)) {
    z <- results[[r]]
    if (is.null(z)) next

    loc_theta[r, ] <- z$loc_theta
    loc_se[r, ] <- z$loc_se
    h[r, ] <- z$h
    boot_success[r] <- z$local_boot_success

    if (!is.null(z$error)) errors <- c(errors, z$error)
  }

  colnames(loc_theta) <- colnames(loc_se) <- c("theta0", "theta1", "theta2")
  colnames(h) <- c("h1", "h2")

  tab <- summarize_local_method(
    theta_mat = loc_theta,
    se_mat = loc_se,
    theta_true = cfg$theta_true,
    alpha = alpha
  )

  tab$scenario <- cfg$label
  tab <- tab[, c(
    "scenario", "parameter", "truth",
    "EST_Local", "BIAS_Local", "SE_Local", "SEE_Local",
    "CP_Local", "n_valid"
  )]

  list(
    scenario = cfg,
    table = tab,
    raw = list(
      loc_theta = loc_theta,
      loc_se = loc_se,
      h = h,
      local_boot_success = boot_success,
      errors = errors
    ),
    diagnostics = list(
      n_valid_local = sum(complete.cases(loc_theta, loc_se)),
      mean_h1 = mean(h[, 1], na.rm = TRUE),
      mean_h2 = mean(h[, 2], na.rm = TRUE),
      mean_local_boot_success = mean(boot_success, na.rm = TRUE),
      min_local_boot_success = if (any(is.finite(boot_success))) {
        min(boot_success[is.finite(boot_success)])
      } else {
        NA_real_
      },
      error_count = length(errors)
    )
  )
}

############################################################
## PARALLEL RUNNER WITH CHECKPOINTING
############################################################

run_scenario_parallel_local <- function(
    scenario_id,
    settings,
    n_cores,
    batch_size,
    resume = TRUE) {

  cfg <- scenario_config(scenario_id)

  checkpoint_file <- paste0(
    OUTPUT_PREFIX, "_", scenario_id, "_checkpoint.rds"
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
      cat("Resuming", scenario_id, "from checkpoint.\n")
    }
  }

  pending <- which(vapply(results, is.null, logical(1)))

  if (length(pending) == 0L) {
    cat(scenario_id, "already complete.\n")
    return(
      summarize_scenario_local(
        results = results,
        scenario_id = scenario_id,
        alpha = CI_ALPHA
      )
    )
  }

  cat("\n============================================================\n")
  cat(cfg$label, "\n")
  cat(
    "Local IPCW only | n^(-1/5) bandwidth |",
    "R =", settings$R,
    "B =", settings$B,
    "n/arm =", settings$n0,
    "tau =", settings$tau,
    "cores =", n_cores,
    "\n"
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
    "GRID_CONT",
    "simulate_covariates_2d",
    "scenario_config",
    "make_constant_nonproportional_2d",
    "make_candidate6",
    "make_observed_terminal",
    "simulate_pwf_holds",
    "simulate_one_fail_arm",
    "simulate_pwf_fails",
    "simulate_scenario",
    "kernel_weights_2d",
    "bandwidth_n15_2d",
    "weighted_km_surv",
    "Ghat_common",
    "S1_local_raw",
    "J_local_raw",
    "wr_local_exact_2d",
    "estimate_local_surface_2d",
    "grid_weights_2d",
    "fit_weighted_plane",
    "fit_loglinear_from_surface",
    "bootstrap_local_coeff_2d",
    "run_one_rep_local_n15"
  )

  parallel::clusterExport(
    cl,
    export_names,
    envir = .GlobalEnv
  )

  batches <- split(
    pending,
    ceiling(seq_along(pending) / batch_size)
  )

  for (idx in batches) {
    ans <- parallel::parLapply(
      cl,
      idx,
      run_one_rep_local_n15,
      scenario_id = scenario_id,
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

    cat(
      "Completed",
      completed,
      "/",
      settings$R,
      "for",
      scenario_id,
      "\n"
    )

    flush.console()
  }

  summarize_scenario_local(
    results = results,
    scenario_id = scenario_id,
    alpha = CI_ALPHA
  )
}

############################################################
## RUN ALL FOUR SCENARIOS
############################################################

cat("\n============================================================\n")
cat("FINAL LOCAL-IPCW COEFFICIENT RERUN\n")
cat("============================================================\n")
cat("Design version:", DESIGN_VERSION, "\n")
cat("Bandwidth: sd(X_j) * n^(-1/5)\n")
cat("Bandwidth computed once in original sample and fixed in bootstrap.\n")
cat("All scenarios use two continuous covariates.\n")
cat("Covariate correlation rho =", RHO_X, "\n")
cat("Common 5 x 5 grid:", paste(round(GRID_CONT, 3), collapse = ", "), "\n")
cat(
  "R =", R_MC,
  "B =", B_BOOT,
  "n/arm =", N0,
  "tau =", TAU,
  "cores =", N_CORES,
  "\n"
)
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
  design_version = DESIGN_VERSION,
  bandwidth_rule = "sd(X_j) * n^(-1/5)",
  seed = NA_integer_
)

all_results <- vector("list", 4L)
names(all_results) <- c("S1", "S2", "S3", "S4")

for (sid in names(all_results)) {
  settings$seed <- unname(BASE_SEEDS[sid])

  all_results[[sid]] <- run_scenario_parallel_local(
    scenario_id = sid,
    settings = settings,
    n_cores = N_CORES,
    batch_size = BATCH_SIZE,
    resume = RESUME
  )

  saveRDS(
    all_results,
    paste0(OUTPUT_PREFIX, "_partial.rds")
  )
}

############################################################
## COMBINE AND PRINT FINAL LOCAL TABLE
############################################################

final_local_table <- do.call(
  rbind,
  lapply(all_results, function(z) z$table)
)

rownames(final_local_table) <- NULL

cat("\n\n============================================================\n")
cat("FINAL LOCAL-IPCW COEFFICIENT TABLE: n^(-1/5)\n")
cat("============================================================\n")

print_table <- final_local_table
num <- vapply(print_table, is.numeric, logical(1))
print_table[num] <- lapply(
  print_table[num],
  function(z) round(z, 4)
)

print(
  print_table,
  row.names = FALSE
)

############################################################
## DIAGNOSTICS
############################################################

cat("\n\n============================================================\n")
cat("DIAGNOSTICS\n")
cat("============================================================\n")

for (sid in names(all_results)) {
  cat(
    "\n",
    sid,
    ": ",
    all_results[[sid]]$scenario$label,
    "\n",
    sep = ""
  )

  print(
    all_results[[sid]]$diagnostics
  )
}

############################################################
## INTERNAL REPRODUCTION CHECK FOR SCENARIO 3
############################################################

s3 <- final_local_table[
  final_local_table$scenario == "Log-linear WR; PWF holds",
  ,
  drop = FALSE
]

cat("\n\n============================================================\n")
cat("SCENARIO 3 REPRODUCTION CHECK\n")
cat("============================================================\n")
cat(
  "Expected from the separate n^(-1/5) bootstrap diagnostic:\n",
  "  intercept: EST about 0.41259, bias about 0.00712, CP about 0.952\n",
  "  x1 slope : EST about 0.47426, bias about -0.02574, CP about 0.950\n",
  "  x2 slope : EST about 0.22845, bias about -0.02155, CP about 0.960\n",
  sep = ""
)
cat("\nCurrent combined-run Scenario 3 rows:\n")

s3_print <- s3
num_s3 <- vapply(s3_print, is.numeric, logical(1))
s3_print[num_s3] <- lapply(
  s3_print[num_s3],
  function(z) round(z, 5)
)

print(
  s3_print,
  row.names = FALSE
)

############################################################
## SAVE FINAL RESULTS
############################################################

final_output <- list(
  settings = list(
    R = R_MC,
    B = B_BOOT,
    n0 = N0,
    n1 = N1,
    tau = TAU,
    tmax = TMAX,
    censor_rate = CENSOR_RATE,
    rho_x = RHO_X,
    n_cores = N_CORES,
    grid_cont = GRID_CONT,
    bandwidth_rule = "sd(X_j) * n^(-1/5)",
    design_version = DESIGN_VERSION
  ),
  scenarios = lapply(
    c("S1", "S2", "S3", "S4"),
    scenario_config
  ),
  final_local_table = final_local_table,
  results = all_results
)

rds_file <- paste0(
  OUTPUT_PREFIX,
  ".rds"
)

csv_file <- paste0(
  OUTPUT_PREFIX,
  "_table.csv"
)

saveRDS(
  final_output,
  rds_file
)

write.csv(
  final_local_table,
  csv_file,
  row.names = FALSE
)

cat("\nSaved:\n")
cat("  ", rds_file, "\n", sep = "")
cat("  ", csv_file, "\n", sep = "")

cat("\nTotal elapsed time:\n")
print(
  Sys.time() - time0
)


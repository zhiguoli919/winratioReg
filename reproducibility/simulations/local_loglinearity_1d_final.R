############################################################
## REVISED 1D TEST SIMULATIONS
##
## Constancy and log-linearity tests for the covariate-specific
## log win ratio with ONE continuous covariate.
##
## CURRENT exact-event-time Local IPCW estimator:
##   - latent T2 retained
##   - U1 = min(T1, C)
##   - Delta1 = I(T1 <= C)
##   - U2 = min(T2, T1, C)
##   - Delta2 = I(T2 <= min(T1, C))
##   - one common censoring survival G_a within each arm
##   - exact event-time Stieltjes estimator
##   - no artificial time-integration grid
##
## REVISED 1D DESIGN TO EXAMINE POWER
## ------------------------------------------------------------
## Sample sizes:
##   n0 = n1 = 300, 500, 800
##
## Monte Carlo / bootstrap:
##   R = 500
##   B = 200
##   alpha = 0.05
##
## Tests:
##   Constancy: H0, H11 (mild), H12 (strong)
##   Linearity: H0, H11 (mild), H12 (strong)
##
## Test grid:
##   30 equally spaced points over [-1, 1]
##
## Event/censoring design:
##   tau = 5
##   tmax = 8
##   censoring baseline rate = 0.05
##
## Bandwidth:
##   h = 1.06 * sd(X) * N^(-1/5),
##   where N is the pooled sample size across the two arms.
##   The bandwidth is computed once in the original sample and
##   held fixed throughout the bootstrap.
##
## Parallelization:
##   - PSOCK cluster, suitable for Windows
##   - Monte Carlo replications parallelized
##   - checkpoints every 25 replications
##   - simple progress messages every 50 replications
##
## OUTPUT DIRECTORY:
##   test_simulation_1d_revised_output/
##
## To resume after interruption:
##   run the same script again with RESUME <- TRUE.
############################################################

time0 <- Sys.time()

suppressPackageStartupMessages({
  library(copula)
  library(parallel)
})

############################################################
## CENTRAL SETTINGS
############################################################

R_FINAL <- 500L
B_BOOT <- 200L

SAMPLE_SIZES <- c(
  300L,
  500L,
  800L
)

TAU <- 5
TMAX <- 8
ALPHA <- 0.05
CENSOR_RATE <- 0.05

GRID_1D <- seq(
  -1,
  1,
  length.out = 30
)

available_cores <- parallel::detectCores(logical = TRUE)
if (!is.finite(available_cores) || available_cores < 2L) {
  available_cores <- 2L
}
N_CORES <- max(1L, min(16L, available_cores - 1L))

BATCH_SIZE <- 25L
PROGRESS_EVERY <- 50L
RESUME <- TRUE

OUTPUT_DIR <- "test_simulation_1d_30grid_n300_500_800_output"
CHECKPOINT_DIR <- file.path(OUTPUT_DIR, "checkpoints")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CHECKPOINT_DIR, showWarnings = FALSE, recursive = TRUE)

DESIGN_VERSION <-
  "revised_1d_tests_exact_event_commonG_30grid_n300_500_800_v1"
BASE_SEED <- 20261105L

############################################################
## 1D KERNEL AND BANDWIDTH
############################################################

kernel_weights_1d <- function(x, x0, h, eps = 1e-12) {
  u <- (x - x0) / h
  w <- dnorm(u)
  sw <- sum(w)

  if (!is.finite(sw) || sw < eps) {
    return(rep(1 / length(x), length(x)))
  }

  w / sw
}
bandwidth_1d <- function(x) {
  n <- length(x)
  h <- 1.06 * sd(x) * n^(-1 / 5)

  if (!is.finite(h) || h <= 0) {
    h <- n^(-1 / 5)
  }

  h
}

############################################################
## COMMON CENSORING SURVIVAL
############################################################

weighted_km_surv <- function(
    time,
    status,
    weights,
    times,
    eps = 1e-8) {

  ord <- order(time)

  time <- time[ord]
  status <- status[ord]
  weights <- weights[ord]

  event_times <- sort(unique(time[status == 1L]))

  if (length(event_times) == 0L) {
    return(rep(1, length(times)))
  }

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

  out <- vapply(
    times,
    function(tt) {
      idx <- which(event_times <= tt)
      if (length(idx) == 0L) {
        1
      } else {
        surv_vals[max(idx)]
      }
    },
    numeric(1)
  )

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
  if (length(t_vec) == 0L) {
    return(numeric(0))
  }

  G <- Ghat_common(t_vec, dat, w)

  vapply(
    seq_along(t_vec),
    function(k) {
      sum(w * (dat$U1 > t_vec[k])) / G[k]
    },
    numeric(1)
  )
}

J_local_raw <- function(tau, u_vec, dat, w) {
  if (length(u_vec) == 0L) {
    return(numeric(0))
  }

  Gtau <- as.numeric(
    Ghat_common(tau, dat, w)
  )

  vapply(
    u_vec,
    function(u) {
      sum(
        w *
          (dat$U1 > tau) *
          (dat$U2 > u)
      ) / Gtau
    },
    numeric(1)
  )
}

############################################################
## CURRENT TERMINAL-EVENT OBSERVATION SCHEME: 1D
############################################################

make_observed_terminal_1d <- function(
    x,
    T1,
    T2,
    C,
    A) {

  U1 <- pmin(T1, C)
  Delta1 <- as.integer(T1 <= C)

  U2 <- pmin(T2, T1, C)
  Delta2 <- as.integer(
    T2 <= pmin(T1, C)
  )

  data.frame(
    x = x,
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
## DGPs: 1D CONSTANCY TEST
############################################################

make_constancy_dgp_1d <- function(
    scenario = c("H0", "mild", "strong")) {

  scenario <- match.arg(scenario)

  if (scenario == "H0") {
    arm1 <- list(
      eta_death_fun = function(x) {
        -0.35 + 0.25 * x
      },
      eta_prog_fun = function(x) {
        -0.25 + 0.15 * x
      },
      lambda_death = 0.10,
      lambda_prog = 0.18,
      theta = 1.2
    )

    arm0 <- list(
      eta_death_fun = function(x) {
        0.25 * x
      },
      eta_prog_fun = function(x) {
        0.15 * x
      },
      lambda_death = 0.10,
      lambda_prog = 0.18,
      theta = 1.2
    )
  }

  if (scenario == "mild") {
    arm1 <- list(
      eta_death_fun = function(x) {
        -0.35 + 0.05 * x + 0.04 * x^2
      },
      eta_prog_fun = function(x) {
        -0.25 + 0.04 * x^2
      },
      lambda_death = 0.10,
      lambda_prog = 0.18,
      theta = 1.2
    )

    arm0 <- list(
      eta_death_fun = function(x) {
        0.25 * x
      },
      eta_prog_fun = function(x) {
        0.15 * x
      },
      lambda_death = 0.10,
      lambda_prog = 0.18,
      theta = 1.2
    )
  }

  if (scenario == "strong") {
    arm1 <- list(
      eta_death_fun = function(x) {
        -0.35 - 0.15 * x + 0.08 * x^2
      },
      eta_prog_fun = function(x) {
        -0.25 - 0.20 * x + 0.08 * x^2
      },
      lambda_death = 0.10,
      lambda_prog = 0.18,
      theta = 1.2
    )

    arm0 <- list(
      eta_death_fun = function(x) {
        0.25 * x
      },
      eta_prog_fun = function(x) {
        0.15 * x
      },
      lambda_death = 0.10,
      lambda_prog = 0.18,
      theta = 1.2
    )
  }

  list(
    arm1 = arm1,
    arm0 = arm0
  )
}

############################################################
## DGPs: 1D LINEARITY TEST
############################################################

make_linearity_dgp_1d <- function(
    scenario = c("H0", "mild", "strong")) {

  scenario <- match.arg(scenario)

  if (scenario == "H0") {
    arm1 <- list(
      eta_death_fun = function(x) {
        -0.25 - 0.20 * x
      },
      eta_prog_fun = function(x) {
        -0.18 - 0.15 * x
      },
      lambda_death = 0.10,
      lambda_prog = 0.18,
      theta = 1.2
    )

    arm0 <- list(
      eta_death_fun = function(x) {
        0
      },
      eta_prog_fun = function(x) {
        0
      },
      lambda_death = 0.10,
      lambda_prog = 0.18,
      theta = 1.2
    )
  }

  if (scenario == "mild") {
    arm1 <- list(
      eta_death_fun = function(x) {
        -0.25 -
          0.20 * x +
          0.25 * x^2 -
          0.12 * x^3
      },
      eta_prog_fun = function(x) {
        -0.18 -
          0.15 * x +
          0.38 * x^2 -
          0.22 * x^3
      },
      lambda_death = 0.10,
      lambda_prog = 0.18,
      theta = 1.2
    )

    arm0 <- list(
      eta_death_fun = function(x) {
        0
      },
      eta_prog_fun = function(x) {
        0
      },
      lambda_death = 0.10,
      lambda_prog = 0.18,
      theta = 1.2
    )
  }

  if (scenario == "strong") {
    arm1 <- list(
      eta_death_fun = function(x) {
        -0.25 -
          0.20 * x +
          0.45 * x^2 -
          0.25 * x^3
      },
      eta_prog_fun = function(x) {
        -0.18 -
          0.15 * x +
          0.70 * x^2 -
          0.45 * x^3
      },
      lambda_death = 0.10,
      lambda_prog = 0.18,
      theta = 1.2
    )

    arm0 <- list(
      eta_death_fun = function(x) {
        0
      },
      eta_prog_fun = function(x) {
        0
      },
      lambda_death = 0.10,
      lambda_prog = 0.18,
      theta = 1.2
    )
  }

  list(
    arm1 = arm1,
    arm0 = arm0
  )
}
\n############################################################
## DGP DISPATCH: 1D ONLY
############################################################

get_dgp_1d <- function(test_type, scenario) {
  if (test_type == "constancy") {
    return(make_constancy_dgp_1d(scenario))
  }
  if (test_type == "linearity") {
    return(make_linearity_dgp_1d(scenario))
  }
  stop("Unknown 1D test_type.")
}

############################################################
## SIMULATE ONE ARM: 1D
############################################################

simulate_arm_1d <- function(
    n,
    arm,
    A,
    seed,
    tmax = TMAX,
    censor_rate = CENSOR_RATE) {

  set.seed(seed)

  x <- rnorm(n)

  cop <- claytonCopula(
    arm$theta,
    dim = 2
  )

  UV <- rCopula(
    n,
    cop
  )

  T1 <- qexp(
    UV[, 1],
    rate = arm$lambda_death *
      exp(arm$eta_death_fun(x))
  )

  T2 <- qexp(
    UV[, 2],
    rate = arm$lambda_prog *
      exp(arm$eta_prog_fun(x))
  )

  censor_eta <- if (A == 1L) {
    0.20 * x
  } else {
    -0.20 * x
  }

  C <- pmin(
    rexp(
      n,
      rate = censor_rate *
        exp(censor_eta)
    ),
    tmax
  )

  make_observed_terminal_1d(
    x = x,
    T1 = T1,
    T2 = T2,
    C = C,
    A = A
  )
}

############################################################
## EXACT EVENT-TIME LOCAL IPCW: 1D
############################################################

wr_local_exact_1d <- function(
    x0,
    dat1,
    dat0,
    h,
    tau) {

  w1 <- kernel_weights_1d(
    dat1$x,
    x0,
    h
  )

  w0 <- kernel_weights_1d(
    dat0$x,
    x0,
    h
  )

  ##########################################################
  ## Primary component
  ##########################################################

  t0 <- sort(
    unique(
      dat0$U1[
        dat0$Delta1 == 1L &
          dat0$U1 <= tau
      ]
    )
  )

  if (length(t0) > 0L) {
    G0 <- Ghat_common(
      t0,
      dat0,
      w0
    )

    dF0 <- vapply(
      t0,
      function(tt) {
        sum(
          w0[
            dat0$Delta1 == 1L &
              dat0$U1 == tt
          ]
        )
      },
      numeric(1)
    ) / G0

    Wp <- sum(
      S1_local_raw(
        t0,
        dat1,
        w1
      ) * dF0
    )
  } else {
    Wp <- 0
  }

  t1 <- sort(
    unique(
      dat1$U1[
        dat1$Delta1 == 1L &
          dat1$U1 <= tau
      ]
    )
  )

  if (length(t1) > 0L) {
    G1 <- Ghat_common(
      t1,
      dat1,
      w1
    )

    dF1 <- vapply(
      t1,
      function(tt) {
        sum(
          w1[
            dat1$Delta1 == 1L &
              dat1$U1 == tt
          ]
        )
      },
      numeric(1)
    ) / G1

    Lp <- sum(
      S1_local_raw(
        t1,
        dat0,
        w0
      ) * dF1
    )
  } else {
    Lp <- 0
  }

  ##########################################################
  ## Secondary component
  ##########################################################

  G1tau <- as.numeric(
    Ghat_common(
      tau,
      dat1,
      w1
    )
  )

  G0tau <- as.numeric(
    Ghat_common(
      tau,
      dat0,
      w0
    )
  )

  u0 <- sort(
    unique(
      dat0$U2[
        dat0$U1 > tau &
          dat0$Delta2 == 1L &
          dat0$U2 <= tau
      ]
    )
  )

  if (length(u0) > 0L) {
    dJ0 <- vapply(
      u0,
      function(uu) {
        sum(
          w0[
            dat0$U1 > tau &
              dat0$Delta2 == 1L &
              dat0$U2 == uu
          ]
        )
      },
      numeric(1)
    ) / G0tau

    Ws <- sum(
      J_local_raw(
        tau,
        u0,
        dat1,
        w1
      ) * dJ0
    )
  } else {
    Ws <- 0
  }

  u1 <- sort(
    unique(
      dat1$U2[
        dat1$U1 > tau &
          dat1$Delta2 == 1L &
          dat1$U2 <= tau
      ]
    )
  )

  if (length(u1) > 0L) {
    dJ1 <- vapply(
      u1,
      function(uu) {
        sum(
          w1[
            dat1$U1 > tau &
              dat1$Delta2 == 1L &
              dat1$U2 == uu
          ]
        )
      },
      numeric(1)
    ) / G1tau

    Ls <- sum(
      J_local_raw(
        tau,
        u1,
        dat0,
        w0
      ) * dJ1
    )
  } else {
    Ls <- 0
  }

  W <- Wp + Ws
  L <- Lp + Ls

  if (
    !is.finite(W) ||
      !is.finite(L) ||
      W <= 0 ||
      L <= 0
  ) {
    return(NA_real_)
  }

  W / L
}

estimate_wr_1d <- function(
    dat1,
    dat0,
    x_grid,
    h,
    tau) {

  vapply(
    x_grid,
    function(x0) {
      wr_local_exact_1d(
        x0 = x0,
        dat1 = dat1,
        dat0 = dat0,
        h = h,
        tau = tau
      )
    },
    numeric(1)
  )
}

############################################################
## 1D GRID WEIGHTS AND LINE PROJECTION
############################################################

trapz_weights <- function(x) {
  J <- length(x)

  if (J < 2L) {
    return(1)
  }

  w <- numeric(J)

  w[1] <- (x[2] - x[1]) / 2
  w[J] <- (x[J] - x[J - 1]) / 2

  if (J > 2L) {
    for (j in 2:(J - 1L)) {
      w[j] <-
        (x[j + 1L] - x[j - 1L]) / 2
    }
  }

  w / sum(w)
}
fit_weighted_line <- function(
    y,
    x,
    w) {

  X <- cbind(1, x)

  XtWX <- crossprod(
    X,
    w * X
  )

  if (qr(XtWX)$rank < 2L) {
    return(NULL)
  }

  beta <- solve(
    XtWX,
    crossprod(
      X,
      w * y
    )
  )

  fitted <- as.vector(
    X %*% beta
  )

  list(
    beta = beta,
    fitted = fitted
  )
}

############################################################
## 1D LACK-OF-FIT STATISTIC
############################################################

test_stat_1d <- function(
    logWR,
    x_grid,
    test_type) {

  if (any(!is.finite(logWR))) {
    return(NULL)
  }

  w <- trapz_weights(
    x_grid
  )

  if (test_type == "constancy") {
    fitted <- rep(
      sum(w * logWR),
      length(logWR)
    )
  } else {
    fit <- fit_weighted_line(
      logWR,
      x_grid,
      w
    )

    if (is.null(fit)) {
      return(NULL)
    }

    fitted <- fit$fitted
  }

  resid <- logWR - fitted

  list(
    D = sum(w * resid^2),
    fitted = fitted,
    resid = resid
  )
}

############################################################
## CENTERED BOOTSTRAP TEST: 1D
############################################################

bootstrap_test_1d <- function(
    dat1,
    dat0,
    test_type,
    x_grid,
    h,
    tau,
    B_boot,
    alpha,
    seed) {

  set.seed(seed)

  WR_hat <- estimate_wr_1d(
    dat1 = dat1,
    dat0 = dat0,
    x_grid = x_grid,
    h = h,
    tau = tau
  )

  if (
    any(!is.finite(WR_hat)) ||
      any(WR_hat <= 0)
  ) {
    return(NULL)
  }

  logWR_hat <- log(
    WR_hat
  )

  obs <- test_stat_1d(
    logWR_hat,
    x_grid,
    test_type
  )

  if (is.null(obs)) {
    return(NULL)
  }

  D_boot <- rep(
    NA_real_,
    B_boot
  )

  n1 <- nrow(dat1)
  n0 <- nrow(dat0)

  for (b in seq_len(B_boot)) {
    id1 <- sample.int(
      n1,
      n1,
      replace = TRUE
    )

    id0 <- sample.int(
      n0,
      n0,
      replace = TRUE
    )

    dat1_b <- dat1[
      id1,
      ,
      drop = FALSE
    ]

    dat0_b <- dat0[
      id0,
      ,
      drop = FALSE
    ]

    WR_b <- try(
      estimate_wr_1d(
        dat1 = dat1_b,
        dat0 = dat0_b,
        x_grid = x_grid,
        h = h,
        tau = tau
      ),
      silent = TRUE
    )

    if (
      inherits(WR_b, "try-error") ||
        any(!is.finite(WR_b)) ||
        any(WR_b <= 0)
    ) {
      next
    }

    diff_log <-
      log(WR_b) -
      logWR_hat

    boot_stat <- test_stat_1d(
      diff_log,
      x_grid,
      test_type
    )

    if (!is.null(boot_stat)) {
      D_boot[b] <- boot_stat$D
    }
  }

  D_boot_fin <- D_boot[
    is.finite(D_boot)
  ]

  min_success <- max(
    20L,
    ceiling(0.80 * B_boot)
  )

  if (length(D_boot_fin) < min_success) {
    return(
      list(
        D_obs = obs$D,
        crit = NA_real_,
        p_value = NA_real_,
        reject = NA_real_,
        n_boot_success = length(D_boot_fin)
      )
    )
  }

  crit <- as.numeric(
    quantile(
      D_boot_fin,
      probs = 1 - alpha,
      names = FALSE
    )
  )

  p_value <- mean(
    D_boot_fin >= obs$D
  )

  reject <- as.numeric(
    obs$D > crit
  )

  list(
    D_obs = obs$D,
    crit = crit,
    p_value = p_value,
    reject = reject,
    n_boot_success = length(D_boot_fin)
  )
}
\n############################################################
## 1D SCENARIO DEFINITIONS
############################################################

scenario_table <- data.frame(
  scenario_id = c(
    "C1_H0", "C1_MILD", "C1_STRONG",
    "L1_H0", "L1_MILD", "L1_STRONG"
  ),
  dimension = rep("1D", 6),
  test_type = c(rep("constancy", 3), rep("linearity", 3)),
  scenario = rep(c("H0", "mild", "strong"), 2),
  stringsAsFactors = FALSE
)

############################################################
## ONE MONTE CARLO REPLICATION: 1D ONLY
############################################################

run_one_rep <- function(r, scenario_def, settings) {
  time_start <- Sys.time()
  out <- list(
    r = r,
    reject = NA_real_,
    p_value = NA_real_,
    D_obs = NA_real_,
    crit = NA_real_,
    h1 = NA_real_,
    h2 = NA_real_,
    n_boot_success = NA_real_,
    elapsed_min = NA_real_,
    error = NULL
  )

  test_type <- scenario_def$test_type
  scenario <- scenario_def$scenario

  dgp <- try(get_dgp_1d(test_type, scenario), silent = TRUE)
  if (inherits(dgp, "try-error")) {
    out$error <- "1D DGP construction failed"
    return(out)
  }

  seed_base <- settings$seed + r

  dat1 <- try(
    simulate_arm_1d(
      n = settings$n1,
      arm = dgp$arm1,
      A = 1L,
      seed = seed_base + 1000L,
      tmax = settings$tmax,
      censor_rate = settings$censor_rate
    ),
    silent = TRUE
  )

  dat0 <- try(
    simulate_arm_1d(
      n = settings$n0,
      arm = dgp$arm0,
      A = 0L,
      seed = seed_base + 2000L,
      tmax = settings$tmax,
      censor_rate = settings$censor_rate
    ),
    silent = TRUE
  )

  if (inherits(dat1, "try-error") || inherits(dat0, "try-error")) {
    out$error <- "1D simulation failed"
    return(out)
  }

  h <- bandwidth_1d(c(dat1$x, dat0$x))
  out$h1 <- h

  tst <- try(
    bootstrap_test_1d(
      dat1 = dat1,
      dat0 = dat0,
      test_type = test_type,
      x_grid = settings$grid_1d,
      h = h,
      tau = settings$tau,
      B_boot = settings$B,
      alpha = settings$alpha,
      seed = seed_base + 500000L
    ),
    silent = TRUE
  )

  if (inherits(tst, "try-error") || is.null(tst)) {
    out$error <- "bootstrap test failed"
    out$elapsed_min <- as.numeric(difftime(Sys.time(), time_start, units = "mins"))
    return(out)
  }

  out$reject <- tst$reject
  out$p_value <- tst$p_value
  out$D_obs <- tst$D_obs
  out$crit <- tst$crit
  out$n_boot_success <- tst$n_boot_success
  out$elapsed_min <- as.numeric(difftime(Sys.time(), time_start, units = "mins"))
  out
}

############################################################
## SUMMARIZE ONE SCENARIO
############################################################

summarize_scenario <- function(
    results,
    scenario_def,
    settings) {

  get_num <- function(name) {
    vapply(
      results,
      function(z) {
        if (is.null(z)) {
          NA_real_
        } else {
          as.numeric(z[[name]])
        }
      },
      numeric(1)
    )
  }

  reject <- get_num("reject")
  p_value <- get_num("p_value")
  D_obs <- get_num("D_obs")
  crit <- get_num("crit")
  h1 <- get_num("h1")
  h2 <- get_num("h2")
  n_boot_success <- get_num("n_boot_success")
  elapsed_min <- get_num("elapsed_min")

  error_vec <- vapply(
    results,
    function(z) {
      if (
        is.null(z) ||
          is.null(z$error)
      ) {
        ""
      } else {
        z$error
      }
    },
    character(1)
  )

  good <- is.finite(
    reject
  )

  n_valid <- sum(good)

  rejection_rate <- if (n_valid > 0L) {
    mean(
      reject[good]
    )
  } else {
    NA_real_
  }

  mcse_rejection <- if (
    n_valid > 0L &&
      is.finite(rejection_rate)
  ) {
    sqrt(
      rejection_rate *
        (1 - rejection_rate) /
        n_valid
    )
  } else {
    NA_real_
  }

  data.frame(
    scenario_id = scenario_def$scenario_id,
    dimension = scenario_def$dimension,
    test = scenario_def$test_type,
    scenario = scenario_def$scenario,
    n_per_arm = settings$n1,
    tau = settings$tau,
    R = settings$R,
    B = settings$B,
    rejection_rate = rejection_rate,
    mcse_rejection = mcse_rejection,
    mean_p_value = mean(
      p_value,
      na.rm = TRUE
    ),
    median_p_value = median(
      p_value,
      na.rm = TRUE
    ),
    mean_D_obs = mean(
      D_obs,
      na.rm = TRUE
    ),
    mean_crit = mean(
      crit,
      na.rm = TRUE
    ),
    mean_h1 = mean(
      h1,
      na.rm = TRUE
    ),
    mean_h2 = if (
      any(is.finite(h2))
    ) {
      mean(
        h2,
        na.rm = TRUE
      )
    } else {
      NA_real_
    },
    mean_boot_success = mean(
      n_boot_success,
      na.rm = TRUE
    ),
    min_boot_success = if (
      any(is.finite(n_boot_success))
    ) {
      min(
        n_boot_success,
        na.rm = TRUE
      )
    } else {
      NA_real_
    },
    n_valid = n_valid,
    error_count = sum(
      nzchar(error_vec)
    ),
    mean_elapsed_min = mean(
      elapsed_min,
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
}
\n############################################################
## PARALLEL SCENARIO RUNNER WITH CHECKPOINTING
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
    if (
      !inherits(old, "try-error") &&
      is.list(old) &&
      identical(old$scenario_id, scenario_id) &&
      identical(old$design_version, DESIGN_VERSION) &&
      identical(old$settings, settings) &&
      length(old$results) == settings$R
    ) {
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
    "Running", scenario_id, "|", scenario_def$test_type, "|",
    scenario_def$scenario, "| n/arm =", n_per_arm, "\n"
  )
  cat(
    "R =", settings$R, "| B =", settings$B,
    "| tau =", settings$tau, "| cores =", n_cores, "\n"
  )
  cat("============================================================\n")

  cl <- parallel::makeCluster(n_cores, type = "PSOCK")
  on.exit(parallel::stopCluster(cl), add = TRUE)

  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages(library(copula))
    NULL
  })

  export_names <- c(
    "kernel_weights_1d",
    "bandwidth_1d",
    "weighted_km_surv",
    "Ghat_common",
    "S1_local_raw",
    "J_local_raw",
    "make_observed_terminal_1d",
    "make_constancy_dgp_1d",
    "make_linearity_dgp_1d",
    "get_dgp_1d",
    "simulate_arm_1d",
    "wr_local_exact_1d",
    "estimate_wr_1d",
    "trapz_weights",
    "fit_weighted_line",
    "test_stat_1d",
    "bootstrap_test_1d",
    "run_one_rep",
    "TAU",
    "TMAX",
    "CENSOR_RATE"
  )

  parallel::clusterExport(cl, export_names, envir = .GlobalEnv)

  batches <- split(pending, ceiling(seq_along(pending) / batch_size))
  last_reported <- sum(!vapply(results, is.null, logical(1)))

  for (idx in batches) {
    ans <- parallel::parLapplyLB(
      cl,
      idx,
      run_one_rep,
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
        "for", scenario_id, "n =", n_per_arm, "\n"
      )
      flush.console()
      last_reported <- completed
    }
  }

  list(
    summary = summarize_scenario(results, scenario_def, settings),
    results = results
  )
}

############################################################
## RUN REVISED 1D SIMULATION
############################################################

cat("\n============================================================\n")
cat("REVISED 1D TEST SIMULATIONS: EXACT-EVENT LOCAL IPCW\n")
cat("============================================================\n")
cat("Sample sizes per arm:", paste(SAMPLE_SIZES, collapse = ", "), "\n")
cat("R =", R_FINAL, "| B =", B_BOOT, "| tau =", TAU, "| alpha =", ALPHA, "\n")
cat("1D test grid:", length(GRID_1D), "equally spaced points on [-1,1]\n")
cat("1D bandwidth: 1.06 * sd(X) * pooled_N^(-1/5)\n")
cat("Bandwidth held fixed throughout each bootstrap.\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("============================================================\n\n")

run_table <- expand.grid(
  n_per_arm = SAMPLE_SIZES,
  scenario_row = seq_len(nrow(scenario_table)),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
run_table$test_type <- scenario_table$test_type[run_table$scenario_row]
run_table$scenario <- scenario_table$scenario[run_table$scenario_row]
run_table$scenario_id <- scenario_table$scenario_id[run_table$scenario_row]

test_order <- match(run_table$test_type, c("constancy", "linearity"))
scenario_order <- match(run_table$scenario, c("H0", "mild", "strong"))
run_table <- run_table[
  order(test_order, scenario_order, run_table$n_per_arm),
  ,
  drop = FALSE
]
rownames(run_table) <- NULL

all_results <- vector("list", nrow(run_table))
summary_list <- vector("list", nrow(run_table))
run_keys <- paste0(run_table$scenario_id, "_n", run_table$n_per_arm)
names(all_results) <- run_keys

for (ii in seq_len(nrow(run_table))) {
  rr <- run_table[ii, ]
  def <- as.list(scenario_table[rr$scenario_row, , drop = FALSE])
  n_now <- as.integer(rr$n_per_arm)

  settings <- list(
    R = R_FINAL,
    B = B_BOOT,
    n0 = n_now,
    n1 = n_now,
    tau = TAU,
    tmax = TMAX,
    alpha = ALPHA,
    censor_rate = CENSOR_RATE,
    grid_1d = GRID_1D,
    design_version = DESIGN_VERSION,
    seed = BASE_SEED + n_now * 1000000L + rr$scenario_row * 100000L
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
      run_table = run_table,
      completed_results = all_results,
      completed_summaries = summary_list
    ),
    file.path(OUTPUT_DIR, "partial_state.rds")
  )
}

final_summary <- do.call(rbind, summary_list)
rownames(final_summary) <- NULL
final_summary$hypothesis <- ifelse(
  final_summary$scenario == "H0",
  "H0",
  ifelse(final_summary$scenario == "mild", "H11", "H12")
)

test_ord <- match(final_summary$test, c("constancy", "linearity"))
hyp_ord <- match(final_summary$hypothesis, c("H0", "H11", "H12"))
final_summary <- final_summary[
  order(final_summary$n_per_arm, test_ord, hyp_ord),
  ,
  drop = FALSE
]
rownames(final_summary) <- NULL

compact <- final_summary[
  ,
  c(
    "dimension", "n_per_arm", "test", "hypothesis",
    "rejection_rate", "mcse_rejection", "mean_p_value",
    "mean_D_obs", "mean_crit", "mean_h1", "mean_h2",
    "mean_boot_success", "min_boot_success", "n_valid", "error_count"
  )
]

wide_rows <- list()
for (j in seq_along(SAMPLE_SIZES)) {
  n_now <- SAMPLE_SIZES[j]
  d <- final_summary[final_summary$n_per_arm == n_now, , drop = FALSE]

  get_rate <- function(test_now, hyp_now) {
    z <- d$rejection_rate[d$test == test_now & d$hypothesis == hyp_now]
    if (length(z) != 1L) return(NA_real_)
    z
  }

  wide_rows[[j]] <- data.frame(
    dimension = "1D",
    n_per_arm = n_now,
    const_H0 = get_rate("constancy", "H0"),
    const_H11 = get_rate("constancy", "H11"),
    const_H12 = get_rate("constancy", "H12"),
    linear_H0 = get_rate("linearity", "H0"),
    linear_H11 = get_rate("linearity", "H11"),
    linear_H12 = get_rate("linearity", "H12"),
    stringsAsFactors = FALSE
  )
}

table_wide <- do.call(rbind, wide_rows)
rownames(table_wide) <- NULL

cat("\n\n============================================================\n")
cat("REVISED 1D TEST-SIMULATION RESULTS\n")
cat("============================================================\n")
compact_print <- compact
numeric_cols <- vapply(compact_print, is.numeric, logical(1))
compact_print[numeric_cols] <- lapply(
  compact_print[numeric_cols],
  function(z) round(z, 4)
)
print(compact_print, row.names = FALSE)

cat("\n\n============================================================\n")
cat("MANUSCRIPT-STYLE 1D REJECTION-RATE TABLE\n")
cat("============================================================\n")
table_wide_print <- table_wide
numeric_cols <- vapply(table_wide_print, is.numeric, logical(1))
table_wide_print[numeric_cols] <- lapply(
  table_wide_print[numeric_cols],
  function(z) round(z, 3)
)
print(table_wide_print, row.names = FALSE)

summary_file <- file.path(OUTPUT_DIR, "revised_1d_test_simulation_summary.csv")
compact_file <- file.path(OUTPUT_DIR, "revised_1d_test_simulation_compact.csv")
wide_file <- file.path(OUTPUT_DIR, "revised_1d_test_simulation_table_wide.csv")
rds_file <- file.path(OUTPUT_DIR, "revised_1d_test_simulation_results.rds")

write.csv(final_summary, summary_file, row.names = FALSE)
write.csv(compact, compact_file, row.names = FALSE)
write.csv(table_wide, wide_file, row.names = FALSE)

saveRDS(
  list(
    design_version = DESIGN_VERSION,
    settings = list(
      R = R_FINAL,
      B = B_BOOT,
      sample_sizes = SAMPLE_SIZES,
      tau = TAU,
      tmax = TMAX,
      alpha = ALPHA,
      censor_rate = CENSOR_RATE,
      grid_1d = GRID_1D,
      n_cores = N_CORES
    ),
    run_table = run_table,
    scenario_table = scenario_table,
    summary = final_summary,
    compact = compact,
    table_wide = table_wide,
    results = all_results
  ),
  rds_file
)

cat("\nSaved files:\n")
cat("  ", normalizePath(summary_file, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("  ", normalizePath(compact_file, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("  ", normalizePath(wide_file, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("  ", normalizePath(rds_file, winslash = "/", mustWork = FALSE), "\n", sep = "")

cat("\nTotal elapsed time:\n")
print(Sys.time() - time0)

############################################################
## winratioReg: Stage 1 Local IPCW numerical validation
##
## Purpose
## -------
## Compare the package implementation of Local IPCW with direct
## reference implementations copied from the finalized manuscript
## simulation algorithms.
##
## The script checks BOTH:
##   1. one continuous covariate (final 1D log-linearity design), and
##   2. two continuous covariates (final coefficient-estimation design).
##
## It compares:
##   - rule-of-thumb bandwidths,
##   - pointwise win-ratio estimates,
##   - log win-ratio estimates,
##   - best linear projection coefficients,
##   - unscaled lack-of-fit criterion D,
##   - manuscript-scaled statistic N * prod(h) * D.
##
## Run from the root of the winratioReg RStudio project:
##   source("inst/validation/local_ipcw_validation.R")
############################################################

options(stringsAsFactors = FALSE)

cat("\n============================================================\n")
cat("winratioReg LOCAL IPCW NUMERICAL VALIDATION\n")
cat("============================================================\n\n")

## ---------------------------------------------------------
## 0. PRE-FLIGHT CHECKS
## ---------------------------------------------------------

if (!file.exists("DESCRIPTION")) {
  stop(
    "Run this script from the root of the winratioReg project (the folder containing DESCRIPTION).",
    call. = FALSE
  )
}

needed <- c("devtools", "copula", "MASS")
missing_pkgs <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0L) {
  stop(
    "Please install the following package(s) first: ",
    paste(missing_pkgs, collapse = ", "),
    call. = FALSE
  )
}

devtools::load_all(".", quiet = TRUE)

## Tight enough to detect a coding discrepancy, while allowing harmless
## floating-point differences between algebraically identical calculations.
TOL <- 1e-10

## ---------------------------------------------------------
## 1. REFERENCE FUNCTIONS FROM FINAL LOCAL-IPCW SCRIPTS
## ---------------------------------------------------------

ref_kernel_weights_1d <- function(x, x0, h, eps = 1e-12) {
  u <- (x - x0) / h
  w <- stats::dnorm(u)
  sw <- sum(w)
  if (!is.finite(sw) || sw < eps) return(rep(1 / length(x), length(x)))
  w / sw
}

ref_kernel_weights_2d <- function(X, x0, h, eps = 1e-12) {
  u1 <- (X[, 1] - x0[1]) / h[1]
  u2 <- (X[, 2] - x0[2]) / h[2]
  w <- stats::dnorm(u1) * stats::dnorm(u2)
  sw <- sum(w)
  if (!is.finite(sw) || sw < eps) return(rep(1 / nrow(X), nrow(X)))
  w / sw
}

ref_weighted_km_surv <- function(time, status, weights, times, eps = 1e-8) {
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
    if (is.finite(risk) && risk > eps) S <- S * (1 - event / risk)
    surv_vals[j] <- max(S, eps)
  }

  out <- vapply(times, function(tt) {
    idx <- which(event_times <= tt)
    if (length(idx) == 0L) 1 else surv_vals[max(idx)]
  }, numeric(1))

  pmax(out, eps)
}

ref_Ghat_common <- function(t_vec, dat, w) {
  ref_weighted_km_surv(
    time = dat$U1,
    status = 1L - dat$Delta1,
    weights = w,
    times = t_vec
  )
}

ref_S1_local_raw <- function(t_vec, dat, w) {
  if (length(t_vec) == 0L) return(numeric(0))
  G <- ref_Ghat_common(t_vec, dat, w)
  vapply(seq_along(t_vec), function(k) {
    sum(w * (dat$U1 > t_vec[k])) / G[k]
  }, numeric(1))
}

ref_J_local_raw <- function(tau, u_vec, dat, w) {
  if (length(u_vec) == 0L) return(numeric(0))
  Gtau <- as.numeric(ref_Ghat_common(tau, dat, w))
  vapply(u_vec, function(u) {
    sum(w * (dat$U1 > tau) * (dat$U2 > u)) / Gtau
  }, numeric(1))
}

ref_wr_from_weights <- function(dat1, dat0, w1, w0, tau) {
  ## Primary component
  t0 <- sort(unique(dat0$U1[dat0$Delta1 == 1L & dat0$U1 <= tau]))
  if (length(t0) > 0L) {
    G0 <- ref_Ghat_common(t0, dat0, w0)
    dF0 <- vapply(t0, function(tt) {
      sum(w0[dat0$Delta1 == 1L & dat0$U1 == tt])
    }, numeric(1)) / G0
    Wp <- sum(ref_S1_local_raw(t0, dat1, w1) * dF0)
  } else {
    Wp <- 0
  }

  t1 <- sort(unique(dat1$U1[dat1$Delta1 == 1L & dat1$U1 <= tau]))
  if (length(t1) > 0L) {
    G1 <- ref_Ghat_common(t1, dat1, w1)
    dF1 <- vapply(t1, function(tt) {
      sum(w1[dat1$Delta1 == 1L & dat1$U1 == tt])
    }, numeric(1)) / G1
    Lp <- sum(ref_S1_local_raw(t1, dat0, w0) * dF1)
  } else {
    Lp <- 0
  }

  ## Secondary component
  G1tau <- as.numeric(ref_Ghat_common(tau, dat1, w1))
  G0tau <- as.numeric(ref_Ghat_common(tau, dat0, w0))

  u0 <- sort(unique(dat0$U2[
    dat0$U1 > tau & dat0$Delta2 == 1L & dat0$U2 <= tau
  ]))
  if (length(u0) > 0L) {
    dJ0 <- vapply(u0, function(uu) {
      sum(w0[
        dat0$U1 > tau & dat0$Delta2 == 1L & dat0$U2 == uu
      ])
    }, numeric(1)) / G0tau
    Ws <- sum(ref_J_local_raw(tau, u0, dat1, w1) * dJ0)
  } else {
    Ws <- 0
  }

  u1 <- sort(unique(dat1$U2[
    dat1$U1 > tau & dat1$Delta2 == 1L & dat1$U2 <= tau
  ]))
  if (length(u1) > 0L) {
    dJ1 <- vapply(u1, function(uu) {
      sum(w1[
        dat1$U1 > tau & dat1$Delta2 == 1L & dat1$U2 == uu
      ])
    }, numeric(1)) / G1tau
    Ls <- sum(ref_J_local_raw(tau, u1, dat0, w0) * dJ1)
  } else {
    Ls <- 0
  }

  W <- Wp + Ws
  L <- Lp + Ls
  if (!is.finite(W) || !is.finite(L) || W <= 0 || L <= 0) return(NA_real_)
  W / L
}

ref_wr_1d <- function(x0, dat1, dat0, h, tau) {
  w1 <- ref_kernel_weights_1d(dat1$x, x0, h)
  w0 <- ref_kernel_weights_1d(dat0$x, x0, h)
  ref_wr_from_weights(dat1, dat0, w1, w0, tau)
}

ref_estimate_1d <- function(dat1, dat0, grid, h, tau) {
  vapply(grid, ref_wr_1d, numeric(1), dat1 = dat1, dat0 = dat0, h = h, tau = tau)
}

ref_wr_2d <- function(x0, dat1, dat0, h, tau) {
  X1 <- as.matrix(dat1[, c("x1", "x2")])
  X0 <- as.matrix(dat0[, c("x1", "x2")])
  w1 <- ref_kernel_weights_2d(X1, x0, h)
  w0 <- ref_kernel_weights_2d(X0, x0, h)
  ref_wr_from_weights(dat1, dat0, w1, w0, tau)
}

ref_estimate_2d <- function(dat1, dat0, grid, h, tau) {
  apply(grid, 1L, function(x0) {
    ref_wr_2d(as.numeric(x0), dat1, dat0, h, tau)
  })
}

ref_trapz_weights <- function(x) {
  J <- length(x)
  if (J < 2L) return(1)
  w <- numeric(J)
  w[1L] <- (x[2L] - x[1L]) / 2
  w[J] <- (x[J] - x[J - 1L]) / 2
  if (J > 2L) {
    for (j in 2L:(J - 1L)) w[j] <- (x[j + 1L] - x[j - 1L]) / 2
  }
  w / sum(w)
}

ref_project_1d <- function(logWR, x) {
  w <- ref_trapz_weights(x)
  X <- cbind(1, x)
  beta <- solve(crossprod(X, w * X), crossprod(X, w * logWR))
  fitted <- as.vector(X %*% beta)
  list(beta = as.numeric(beta), D = sum(w * (logWR - fitted)^2), weights = w)
}

ref_project_2d <- function(logWR, grid) {
  w <- rep(1 / nrow(grid), nrow(grid))
  X <- cbind(1, grid[, 1], grid[, 2])
  beta <- solve(crossprod(X, w * X), crossprod(X, w * logWR))
  fitted <- as.vector(X %*% beta)
  list(beta = as.numeric(beta), D = sum(w * (logWR - fitted)^2), weights = w)
}

## ---------------------------------------------------------
## 2. FINAL 1D LINEARITY-DESIGN DGP: FIXED DATASET
## ---------------------------------------------------------

make_observed_1d <- function(x, T1, T2, C, A) {
  data.frame(
    x = x,
    A = A,
    U1 = pmin(T1, C),
    Delta1 = as.integer(T1 <= C),
    U2 = pmin(T2, T1, C),
    Delta2 = as.integer(T2 <= pmin(T1, C))
  )
}

simulate_1d_arm <- function(n, A, seed, scenario = "mild",
                            tmax = 8, censor_rate = 0.05) {
  if (!scenario %in% c("H0", "mild", "strong")) stop("Unknown 1D scenario.")

  if (scenario == "H0") {
    eta1 <- function(x) -0.25 - 0.20 * x
    eta2 <- function(x) -0.18 - 0.15 * x
  } else if (scenario == "mild") {
    eta1 <- function(x) -0.25 - 0.20 * x + 0.25 * x^2 - 0.12 * x^3
    eta2 <- function(x) -0.18 - 0.15 * x + 0.38 * x^2 - 0.22 * x^3
  } else {
    eta1 <- function(x) -0.25 - 0.20 * x + 0.45 * x^2 - 0.25 * x^3
    eta2 <- function(x) -0.18 - 0.15 * x + 0.70 * x^2 - 0.45 * x^3
  }

  ## Final design has no covariate effect in the control-arm event hazards.
  if (A == 0L) {
    eta1 <- function(x) 0 * x
    eta2 <- function(x) 0 * x
  }

  set.seed(seed)
  x <- stats::rnorm(n)
  cop <- copula::claytonCopula(1.2, dim = 2)
  UV <- copula::rCopula(n, cop)
  T1 <- stats::qexp(UV[, 1], rate = 0.10 * exp(eta1(x)))
  T2 <- stats::qexp(UV[, 2], rate = 0.18 * exp(eta2(x)))
  censor_eta <- if (A == 1L) 0.20 * x else -0.20 * x
  C <- pmin(stats::rexp(n, rate = censor_rate * exp(censor_eta)), tmax)
  make_observed_1d(x, T1, T2, C, A)
}

cat("1D validation dataset ...\n")
N1D <- 300L
TAU1D <- 5
GRID1D <- seq(-1, 1, length.out = 30)
D1_1D <- simulate_1d_arm(N1D, 1L, seed = 202611051L, scenario = "mild")
D0_1D <- simulate_1d_arm(N1D, 0L, seed = 202611052L, scenario = "mild")
DAT1D <- rbind(D0_1D, D1_1D)

H1D_REF <- 1.06 * stats::sd(DAT1D$x) * nrow(DAT1D)^(-1 / 5)
WR1D_REF <- ref_estimate_1d(D1_1D, D0_1D, GRID1D, H1D_REF, TAU1D)
P1D_REF <- ref_project_1d(log(WR1D_REF), GRID1D)
T1D_REF <- nrow(DAT1D) * H1D_REF * P1D_REF$D

FIT1D <- winratio_local(
  data = DAT1D,
  arm = "A",
  time1 = "U1",
  status1 = "Delta1",
  time2 = "U2",
  status2 = "Delta2",
  covariates = "x",
  tau = TAU1D,
  treatment = 1,
  control = 0,
  grid = GRID1D,
  bandwidth = "rot",
  B = 0
)

BLP1D_PKG <- blp(FIT1D, type = "linear")$estimate
D1D_PKG <- sum(FIT1D$projection_weights * FIT1D$projection$linear$resid^2)
T1D_PKG <- nrow(DAT1D) * prod(FIT1D$bandwidth) * D1D_PKG

## ---------------------------------------------------------
## 3. FINAL 2D COEFFICIENT-DESIGN DGP (S3): FIXED DATASET
## ---------------------------------------------------------

simulate_2d_s3 <- function(n0 = 250L, n1 = 250L, seed = 20261023L,
                           tmax = 10, censor_rate = 0.05, rho = 0.30) {
  set.seed(seed)
  n <- n0 + n1
  Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
  X <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
  x1 <- X[, 1]
  x2 <- X[, 2]
  A <- sample(c(rep(0L, n0), rep(1L, n1)))

  eta <- 0.3 * x1 - 0.3 * x2 + log(1.5) * A +
    0.5 * A * x1 + 0.25 * A * x2

  cop <- copula::gumbelCopula(param = 2, dim = 2)
  UV <- copula::rCopula(n, cop)

  T1 <- -log(UV[, 1]) * exp(eta) / 0.20
  T2 <- -log(UV[, 2]) * exp(eta) / 2.00
  C <- pmin(stats::rexp(n, rate = censor_rate), tmax)

  dat <- data.frame(
    x1 = x1,
    x2 = x2,
    A = A,
    U1 = pmin(T1, C),
    Delta1 = as.integer(T1 <= C),
    U2 = pmin(T2, T1, C),
    Delta2 = as.integer(T2 <= pmin(T1, C))
  )

  list(
    dat = dat,
    dat1 = dat[dat$A == 1L, , drop = FALSE],
    dat0 = dat[dat$A == 0L, , drop = FALSE]
  )
}

cat("2D validation dataset ...\n")
TAU2D <- 8
GRID2D_VALUES <- seq(-0.75, 0.75, length.out = 5)
GRID2D <- as.matrix(expand.grid(x1 = GRID2D_VALUES, x2 = GRID2D_VALUES))
SIM2D <- simulate_2d_s3()
X2D <- as.matrix(SIM2D$dat[, c("x1", "x2")])
H2D_REF <- apply(X2D, 2, stats::sd) * nrow(X2D)^(-1 / 5)
WR2D_REF <- ref_estimate_2d(SIM2D$dat1, SIM2D$dat0, GRID2D, H2D_REF, TAU2D)
P2D_REF <- ref_project_2d(log(WR2D_REF), GRID2D)
T2D_REF <- nrow(SIM2D$dat) * prod(H2D_REF) * P2D_REF$D

FIT2D <- winratio_local(
  data = SIM2D$dat,
  arm = "A",
  time1 = "U1",
  status1 = "Delta1",
  time2 = "U2",
  status2 = "Delta2",
  covariates = c("x1", "x2"),
  tau = TAU2D,
  treatment = 1,
  control = 0,
  grid = GRID2D,
  bandwidth = "rot",
  B = 0
)

BLP2D_PKG <- blp(FIT2D, type = "linear")$estimate
D2D_PKG <- sum(FIT2D$projection_weights * FIT2D$projection$linear$resid^2)
T2D_PKG <- nrow(SIM2D$dat) * prod(FIT2D$bandwidth) * D2D_PKG

## ---------------------------------------------------------
## 4. COMPARISON TABLE
## ---------------------------------------------------------

max_abs <- function(x, y) max(abs(as.numeric(x) - as.numeric(y)))

results <- data.frame(
  check = c(
    "1D bandwidth",
    "1D WR grid",
    "1D log-WR grid",
    "1D BLP coefficients",
    "1D lack-of-fit D",
    "1D scaled statistic N*h*D",
    "2D bandwidth",
    "2D WR grid",
    "2D log-WR grid",
    "2D BLP coefficients",
    "2D lack-of-fit D",
    "2D scaled statistic N*h1*h2*D"
  ),
  max_abs_difference = c(
    max_abs(H1D_REF, FIT1D$bandwidth),
    max_abs(WR1D_REF, FIT1D$estimates$WR),
    max_abs(log(WR1D_REF), FIT1D$estimates$logWR),
    max_abs(P1D_REF$beta, BLP1D_PKG),
    abs(P1D_REF$D - D1D_PKG),
    abs(T1D_REF - T1D_PKG),
    max_abs(H2D_REF, FIT2D$bandwidth),
    max_abs(WR2D_REF, FIT2D$estimates$WR),
    max_abs(log(WR2D_REF), FIT2D$estimates$logWR),
    max_abs(P2D_REF$beta, BLP2D_PKG),
    abs(P2D_REF$D - D2D_PKG),
    abs(T2D_REF - T2D_PKG)
  ),
  tolerance = TOL,
  stringsAsFactors = FALSE
)

results$pass <- results$max_abs_difference <= results$tolerance

cat("\n============================================================\n")
cat("VALIDATION RESULTS\n")
cat("============================================================\n")
print(results, row.names = FALSE, digits = 12)

cat("\nReference/package BLP coefficients:\n")
cat("  1D reference:", paste(signif(P1D_REF$beta, 12), collapse = ", "), "\n")
cat("  1D package:  ", paste(signif(BLP1D_PKG, 12), collapse = ", "), "\n")
cat("  2D reference:", paste(signif(P2D_REF$beta, 12), collapse = ", "), "\n")
cat("  2D package:  ", paste(signif(BLP2D_PKG, 12), collapse = ", "), "\n")

out_dir <- file.path("inst", "validation")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_file <- file.path(out_dir, "local_ipcw_validation_results.csv")
utils::write.csv(results, out_file, row.names = FALSE)
cat("\nSaved comparison table to:\n  ", normalizePath(out_file, winslash = "/", mustWork = FALSE), "\n", sep = "")

if (all(results$pass)) {
  cat("\n============================================================\n")
  cat("PASS: Local IPCW package reproduces the finalized reference algorithms\n")
  cat("      within tolerance ", format(TOL, scientific = TRUE), ".\n", sep = "")
  cat("============================================================\n\n")
} else {
  cat("\n============================================================\n")
  cat("FAIL: At least one numerical comparison exceeds tolerance.\n")
  cat("      Please send the printed table to ChatGPT before changing code.\n")
  cat("============================================================\n\n")
  stop("Local IPCW numerical validation failed.", call. = FALSE)
}


############################################################
## C9633 UPDATED LOCAL IPCW ANALYSIS: WR(tumor size), PS = 0
##
## Primary event: death
## Secondary event: recurrence/progression before death
## Horizon: tau = 5 years
##
## Key updates relative to the old code:
##   1. ONE common arm-specific local censoring survival G_a
##      estimated from the primary-event follow-up time.
##   2. Exact event-time evaluation of the win/loss integrals.
##   3. No artificial numerical time grid for the integrals.
##   4. 21 covariate grid points from the empirical 5th to 95th percentiles.
##   5. Stability-based bandwidth selector.
##   6. Within-arm nonparametric bootstrap with centered/basic
##      pointwise confidence intervals on the log-WR scale.
##
## This script runs the PS = 0 subgroup analysis only.
############################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
})

############################################################
## 1. SETTINGS
############################################################

file_path <- Sys.getenv("C9633_DATA", unset = "C9633.csv.xlsx")

tau <- 5

ps_value <- 0

B_boot <- 500

seed_boot <- 20261025

grid_length <- 21
grid_lower_q <- 0.05
grid_upper_q <- 0.95

h_mult <- c(0.6, 0.8, 1.0, 1.2, 1.5)
c_ref <- 1.25
tol_bw <- 1.10

alpha <- 0.05

############################################################
## 2. READ AND CONSTRUCT C9633 ENDPOINTS
############################################################

dat_raw <- read_excel(file_path)

dat <- dat_raw %>%
  rename_with(tolower) %>%
  mutate(
    treat = as.numeric(treat),
    age = as.numeric(age),
    tsize = as.numeric(tsize),
    ps = as.numeric(ps),
    survtime = as.numeric(survtime),
    survcens = as.numeric(survcens),
    failtime = as.numeric(failtime),
    failcens = as.numeric(failcens)
  )

required_vars <- c(
  "treat", "tsize", "ps",
  "survtime", "survcens",
  "failtime", "failcens"
)

missing_vars <- setdiff(required_vars, names(dat))

if (length(missing_vars) > 0L) {
  stop(
    "Missing variables: ",
    paste(missing_vars, collapse = ", ")
  )
}

dat <- dat %>%
  mutate(
    death_time = survtime,
    death_event = as.integer(survcens == 1),

    prog_event = as.integer(failcens == 1),

    prog_time = case_when(
      failcens == 1 & !is.na(failtime) ~ failtime,
      failcens %in% c(0, 2) ~ survtime,
      TRUE ~ NA_real_
    )
  ) %>%
  mutate(
    prog_time = ifelse(
      !is.na(prog_time) & !is.na(death_time),
      pmin(prog_time, death_time),
      prog_time
    )
  )

analysis_dat <- dat %>%
  filter(ps == ps_value) %>%
  transmute(
    A = treat,
    x = as.numeric(tsize),
    Y1 = death_time,
    D1 = death_event,
    Y2 = prog_time,
    D2 = prog_event
  ) %>%
  filter(
    A %in% c(0, 1),
    is.finite(x),
    is.finite(Y1),
    is.finite(Y2),
    !is.na(D1),
    !is.na(D2)
  )

if (nrow(analysis_dat) == 0L) {
  stop("No complete observations are available for PS = 0.")
}

dat0 <- analysis_dat %>%
  filter(A == 0) %>%
  select(x, Y1, D1, Y2, D2)

dat1 <- analysis_dat %>%
  filter(A == 1) %>%
  select(x, Y1, D1, Y2, D2)

if (nrow(dat0) == 0L || nrow(dat1) == 0L) {
  stop("Both treatment arms must contain observations.")
}

cat("\n============================================================\n")
cat("UPDATED LOCAL IPCW ANALYSIS: WR(tumor size), PS = 0\n")
cat("============================================================\n")
cat("Total N:", nrow(analysis_dat), "\n")
cat("Control N:", nrow(dat0), "\n")
cat("Treatment N:", nrow(dat1), "\n")
cat("Performance status: PS = 0\n")
cat("Tumor size range:", range(analysis_dat$x), "\n")
cat("tau =", tau, "years\n")
cat("============================================================\n\n")

############################################################
## 3. KERNEL WEIGHTS
############################################################

K_gaussian <- function(u) {
  dnorm(u)
}

kernel_weights <- function(x, x0, h, eps = 1e-14) {

  w <- K_gaussian(
    (x - x0) / h
  )

  sw <- sum(w)

  if (!is.finite(sw) || sw <= eps) {
    return(
      rep(
        1 / length(x),
        length(x)
      )
    )
  }

  w / sw
}

############################################################
## 4. ONE COMMON LOCAL CENSORING SURVIVAL G_a
############################################################

weighted_G_common <- function(
    times,
    Y1,
    D1,
    w,
    eps = 1e-8) {

  times <- as.numeric(times)

  if (length(times) == 0L) {
    return(numeric(0))
  }

  censor <- as.integer(D1 == 0)

  censor_times <- sort(
    unique(
      Y1[censor == 1]
    )
  )

  if (length(censor_times) == 0L) {
    return(
      rep(
        1,
        length(times)
      )
    )
  }

  Gvals <- numeric(
    length(censor_times)
  )

  G <- 1

  for (j in seq_along(censor_times)) {

    tj <- censor_times[j]

    risk <- sum(
      w[Y1 >= tj]
    )

    dcen <- sum(
      w[
        Y1 == tj &
        censor == 1
      ]
    )

    if (is.finite(risk) && risk > eps) {
      frac <- dcen / risk
      frac <- min(
        max(frac, 0),
        1
      )
      G <- G * (1 - frac)
    }

    Gvals[j] <- max(
      G,
      eps
    )
  }

  out <- vapply(
    times,
    function(tt) {

      idx <- which(
        censor_times <= tt
      )

      if (length(idx) == 0L) {
        1
      } else {
        Gvals[max(idx)]
      }
    },
    numeric(1)
  )

  pmax(
    out,
    eps
  )
}

############################################################
## 5. LOCAL IPCW SURVIVAL COMPONENTS
############################################################

S1_local <- function(
    times,
    dat,
    w) {

  if (length(times) == 0L) {
    return(numeric(0))
  }

  G <- weighted_G_common(
    times = times,
    Y1 = dat$Y1,
    D1 = dat$D1,
    w = w
  )

  out <- vapply(
    times,
    function(tt) {
      sum(
        w *
        (dat$Y1 > tt)
      )
    },
    numeric(1)
  ) / G

  pmin(
    pmax(out, 0),
    1
  )
}


dF1_local <- function(
    event_times,
    dat,
    w) {

  if (length(event_times) == 0L) {
    return(numeric(0))
  }

  G <- weighted_G_common(
    times = event_times,
    Y1 = dat$Y1,
    D1 = dat$D1,
    w = w
  )

  mass <- vapply(
    event_times,
    function(tt) {
      sum(
        w *
        (dat$D1 == 1) *
        (dat$Y1 == tt)
      )
    },
    numeric(1)
  )

  mass / G
}


J_local <- function(
    tau,
    u,
    dat,
    w) {

  if (length(u) == 0L) {
    return(numeric(0))
  }

  Gtau <- weighted_G_common(
    times = tau,
    Y1 = dat$Y1,
    D1 = dat$D1,
    w = w
  )[1]

  out <- vapply(
    u,
    function(uu) {
      sum(
        w *
        (dat$Y1 > tau) *
        (dat$Y2 > uu)
      )
    },
    numeric(1)
  ) / Gtau

  pmin(
    pmax(out, 0),
    1
  )
}


minus_dJ_local <- function(
    tau,
    event_times,
    dat,
    w) {

  if (length(event_times) == 0L) {
    return(numeric(0))
  }

  Gtau <- weighted_G_common(
    times = tau,
    Y1 = dat$Y1,
    D1 = dat$D1,
    w = w
  )[1]

  mass <- vapply(
    event_times,
    function(uu) {
      sum(
        w *
        (dat$Y1 > tau) *
        (dat$D2 == 1) *
        (dat$Y2 == uu)
      )
    },
    numeric(1)
  )

  mass / Gtau
}

############################################################
## 6. EXACT EVENT-TIME LOCAL WIN-RATIO ESTIMATOR
############################################################

wr_local_exact <- function(
    x0,
    dat1,
    dat0,
    h,
    tau = 5) {

  w1 <- kernel_weights(
    dat1$x,
    x0,
    h
  )

  w0 <- kernel_weights(
    dat0$x,
    x0,
    h
  )

  e10 <- sort(
    unique(
      dat0$Y1[
        dat0$D1 == 1 &
        dat0$Y1 <= tau
      ]
    )
  )

  e11 <- sort(
    unique(
      dat1$Y1[
        dat1$D1 == 1 &
        dat1$Y1 <= tau
      ]
    )
  )

  if (length(e10) > 0L) {

    W_primary <- sum(
      S1_local(
        e10,
        dat1,
        w1
      ) *
      dF1_local(
        e10,
        dat0,
        w0
      )
    )

  } else {

    W_primary <- 0
  }

  if (length(e11) > 0L) {

    L_primary <- sum(
      S1_local(
        e11,
        dat0,
        w0
      ) *
      dF1_local(
        e11,
        dat1,
        w1
      )
    )

  } else {

    L_primary <- 0
  }

  e20 <- sort(
    unique(
      dat0$Y2[
        dat0$D2 == 1 &
        dat0$Y2 <= tau &
        dat0$Y1 > tau
      ]
    )
  )

  e21 <- sort(
    unique(
      dat1$Y2[
        dat1$D2 == 1 &
        dat1$Y2 <= tau &
        dat1$Y1 > tau
      ]
    )
  )

  if (length(e20) > 0L) {

    W_secondary <- sum(
      J_local(
        tau = tau,
        u = e20,
        dat = dat1,
        w = w1
      ) *
      minus_dJ_local(
        tau = tau,
        event_times = e20,
        dat = dat0,
        w = w0
      )
    )

  } else {

    W_secondary <- 0
  }

  if (length(e21) > 0L) {

    L_secondary <- sum(
      J_local(
        tau = tau,
        u = e21,
        dat = dat0,
        w = w0
      ) *
      minus_dJ_local(
        tau = tau,
        event_times = e21,
        dat = dat1,
        w = w1
      )
    )

  } else {

    L_secondary <- 0
  }

  W <- W_primary + W_secondary
  L <- L_primary + L_secondary

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


estimate_wr_curve_exact <- function(
    dat1,
    dat0,
    x_grid,
    h,
    tau = 5) {

  vapply(
    x_grid,
    function(x0) {
      wr_local_exact(
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
## 7. STABILITY-BASED BANDWIDTH SELECTION
############################################################

select_bandwidth_stability <- function(
    dat1,
    dat0,
    x_all,
    x_grid,
    tau,
    h_mult,
    c_ref = 1.25,
    tol = 1.10,
    verbose = TRUE) {

  n <- length(x_all)

  h0 <- sd(
    x_all,
    na.rm = TRUE
  ) *
    n^(-1 / 5)

  h_vec <- h0 *
    h_mult

  scores <- rep(
    Inf,
    length(h_vec)
  )

  for (k in seq_along(h_vec)) {

    h <- h_vec[k]

    wr_h <- estimate_wr_curve_exact(
      dat1 = dat1,
      dat0 = dat0,
      x_grid = x_grid,
      h = h,
      tau = tau
    )

    wr_ref <- estimate_wr_curve_exact(
      dat1 = dat1,
      dat0 = dat0,
      x_grid = x_grid,
      h = c_ref * h,
      tau = tau
    )

    if (
      all(is.finite(wr_h)) &&
      all(is.finite(wr_ref)) &&
      all(wr_h > 0) &&
      all(wr_ref > 0)
    ) {

      scores[k] <- mean(
        (
          log(wr_h) -
          log(wr_ref)
        )^2
      )
    }

    if (verbose) {
      cat(
        "Bandwidth candidate",
        k,
        "of",
        length(h_vec),
        ": h =",
        round(h, 5),
        "score =",
        signif(scores[k], 6),
        "\n"
      )
    }
  }

  if (!any(is.finite(scores))) {
    stop(
      "All bandwidth candidates failed."
    )
  }

  qmin <- min(
    scores,
    na.rm = TRUE
  )

  eligible <- which(
    is.finite(scores) &
    scores <= tol * qmin
  )

  if (length(eligible) == 0L) {
    stop(
      "No bandwidth satisfied the stability criterion."
    )
  }

  idx <- eligible[1]

  list(
    h0 = h0,
    h_vec = h_vec,
    scores = scores,
    selected_index = idx,
    h = h_vec[idx]
  )
}

############################################################
## 8. GRID AND BANDWIDTH
############################################################

x_all <- analysis_dat$x

x_grid <- seq(
  as.numeric(
    quantile(
      x_all,
      grid_lower_q,
      na.rm = TRUE
    )
  ),
  as.numeric(
    quantile(
      x_all,
      grid_upper_q,
      na.rm = TRUE
    )
  ),
  length.out = grid_length
)

cat(
  "\nGrid range:",
  range(x_grid),
  "\n\n"
)

bw <- select_bandwidth_stability(
  dat1 = dat1,
  dat0 = dat0,
  x_all = x_all,
  x_grid = x_grid,
  tau = tau,
  h_mult = h_mult,
  c_ref = c_ref,
  tol = tol_bw,
  verbose = TRUE
)

h <- bw$h

cat(
  "\nSelected bandwidth h =",
  h,
  "\n\n"
)

############################################################
## 9. ORIGINAL-SAMPLE CURVE
############################################################

WR_hat <- estimate_wr_curve_exact(
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
  stop(
    "Original local IPCW curve contains invalid values."
  )
}

logWR_hat <- log(
  WR_hat
)

############################################################
## 10. FIXED-GRID PROJECTION FUNCTIONS
############################################################

trapz_weights <- function(x) {

  J <- length(x)

  if (J < 2L) {
    stop(
      "At least two grid points are required."
    )
  }

  w <- numeric(J)

  w[1] <-
    (x[2] - x[1]) / 2

  w[J] <-
    (x[J] - x[J - 1]) / 2

  if (J > 2L) {
    for (j in 2:(J - 1L)) {
      w[j] <-
        (x[j + 1L] - x[j - 1L]) / 2
    }
  }

  w / sum(w)
}


project_constant <- function(
    logcurve,
    qweight) {

  alpha_hat <- sum(
    qweight *
    logcurve
  )

  fitted <- rep(
    alpha_hat,
    length(logcurve)
  )

  list(
    alpha = alpha_hat,
    fitted = fitted,
    resid = logcurve - fitted
  )
}


project_linear <- function(
    logcurve,
    x_grid,
    qweight) {

  X <- cbind(
    1,
    x_grid
  )

  WX <- qweight *
    X

  beta <- solve(
    t(X) %*% WX,
    t(X) %*% (
      qweight *
      logcurve
    )
  )

  fitted <- as.vector(
    X %*% beta
  )

  list(
    beta = as.vector(beta),
    fitted = fitted,
    resid = logcurve - fitted
  )
}

qweight <- trapz_weights(
  x_grid
)

const_obs <- project_constant(
  logWR_hat,
  qweight
)

linear_obs <- project_linear(
  logWR_hat,
  x_grid,
  qweight
)

############################################################
## 11. ONE BOOTSTRAP FOR CURVE, CIs, AND TESTS
############################################################

set.seed(
  seed_boot
)

boot_logWR <- matrix(
  NA_real_,
  nrow = B_boot,
  ncol = length(x_grid)
)

success <- logical(
  B_boot
)

for (b in seq_len(B_boot)) {

  id1 <- sample(
    seq_len(nrow(dat1)),
    size = nrow(dat1),
    replace = TRUE
  )

  id0 <- sample(
    seq_len(nrow(dat0)),
    size = nrow(dat0),
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

  wr_b <- try(
    estimate_wr_curve_exact(
      dat1 = dat1_b,
      dat0 = dat0_b,
      x_grid = x_grid,
      h = h,
      tau = tau
    ),
    silent = TRUE
  )

  if (
    !inherits(
      wr_b,
      "try-error"
    ) &&
    all(is.finite(wr_b)) &&
    all(wr_b > 0)
  ) {

    boot_logWR[
      b,
    ] <- log(
      wr_b
    )

    success[b] <- TRUE
  }

  if (
    b %% 50L == 0L ||
    b == B_boot
  ) {

    cat(
      "Completed bootstrap",
      b,
      "of",
      B_boot,
      "- successful:",
      sum(success),
      "\n"
    )

    flush.console()
  }
}

cat(
  "\nSuccessful bootstrap curves:",
  sum(success),
  "of",
  B_boot,
  "\n\n"
)

if (sum(success) < 20L) {
  stop(
    "Too few successful bootstrap curves."
  )
}

boot_diff <- sweep(
  boot_logWR,
  2,
  logWR_hat,
  FUN = "-"
)

q025 <- apply(
  boot_diff,
  2,
  quantile,
  probs = alpha / 2,
  na.rm = TRUE,
  names = FALSE
)

q975 <- apply(
  boot_diff,
  2,
  quantile,
  probs = 1 - alpha / 2,
  na.rm = TRUE,
  names = FALSE
)

basic_log_lower <-
  logWR_hat -
  q975

basic_log_upper <-
  logWR_hat -
  q025

SE_logWR <- apply(
  boot_logWR,
  2,
  sd,
  na.rm = TRUE
)

zcrit <- qnorm(
  1 - alpha / 2
)

normal_log_lower <-
  logWR_hat -
  zcrit *
  SE_logWR

normal_log_upper <-
  logWR_hat +
  zcrit *
  SE_logWR

############################################################
## 12. CONSTANCY AND LINEARITY TESTS
############################################################

n_total <- nrow(
  analysis_dat
)

scale_test <-
  n_total *
  h

D_const_obs <- sum(
  qweight *
  const_obs$resid^2
)

D_linear_obs <- sum(
  qweight *
  linear_obs$resid^2
)

T_const_obs <-
  scale_test *
  D_const_obs

T_linear_obs <-
  scale_test *
  D_linear_obs

T_const_boot <- rep(
  NA_real_,
  B_boot
)

T_linear_boot <- rep(
  NA_real_,
  B_boot
)

alpha_boot <- rep(
  NA_real_,
  B_boot
)

beta_boot <- matrix(
  NA_real_,
  nrow = B_boot,
  ncol = 2
)

for (b in which(success)) {

  log_b <- boot_logWR[
    b,
  ]

  diff_b <-
    log_b -
    logWR_hat

  const_diff <- project_constant(
    diff_b,
    qweight
  )

  linear_diff <- project_linear(
    diff_b,
    x_grid,
    qweight
  )

  T_const_boot[b] <-
    scale_test *
    sum(
      qweight *
      const_diff$resid^2
    )

  T_linear_boot[b] <-
    scale_test *
    sum(
      qweight *
      linear_diff$resid^2
    )

  alpha_boot[b] <-
    project_constant(
      log_b,
      qweight
    )$alpha

  beta_boot[
    b,
  ] <- project_linear(
    log_b,
    x_grid,
    qweight
  )$beta
}

p_constancy <- mean(
  T_const_boot[
    is.finite(T_const_boot)
  ] >=
    T_const_obs
)

p_linearity <- mean(
  T_linear_boot[
    is.finite(T_linear_boot)
  ] >=
    T_linear_obs
)

############################################################
## 13. CENTERED BOOTSTRAP CI FOR CONSTANT PROJECTION
############################################################

alpha_hat <- const_obs$alpha

alpha_diff <-
  alpha_boot -
  alpha_hat

alpha_diff <- alpha_diff[
  is.finite(alpha_diff)
]

alpha_q <- quantile(
  alpha_diff,
  probs = c(
    alpha / 2,
    1 - alpha / 2
  ),
  na.rm = TRUE,
  names = FALSE
)

constant_logWR_CI <- c(
  alpha_hat -
    alpha_q[2],
  alpha_hat -
    alpha_q[1]
)

constant_WR_hat <- exp(
  alpha_hat
)

constant_WR_CI <- exp(
  constant_logWR_CI
)

############################################################
## 14. RESULTS
############################################################

results <- data.frame(
  Tumor_size = x_grid,
  logWR = logWR_hat,
  SE_logWR = SE_logWR,
  WR = WR_hat,
  Normal_lower = exp(normal_log_lower),
  Normal_upper = exp(normal_log_upper),
  Basic_lower = exp(basic_log_lower),
  Basic_upper = exp(basic_log_upper)
)

summary_results <- data.frame(
  Analysis = "WR(tumor size), PS = 0",
  N = n_total,
  N_control = nrow(dat0),
  N_treatment = nrow(dat1),
  tau = tau,
  h0 = bw$h0,
  h = h,
  Constant_logWR = alpha_hat,
  Constant_WR = constant_WR_hat,
  Constant_WR_lower = constant_WR_CI[1],
  Constant_WR_upper = constant_WR_CI[2],
  Linear_intercept = linear_obs$beta[1],
  Linear_slope = linear_obs$beta[2],
  P_constancy = p_constancy,
  P_linearity = p_linearity,
  Bootstrap_success = sum(success)
)

cat(
  "\n============================================================\n"
)

cat(
  "FINAL UPDATED LOCAL IPCW RESULTS: WR(tumor size), PS = 0\n"
)

cat(
  "============================================================\n\n"
)

print(
  results,
  row.names = FALSE
)

cat(
  "\nSummary:\n\n"
)

print(
  summary_results,
  row.names = FALSE
)

############################################################
## 15. PLOT
############################################################

ylim_plot <- range(
  c(
    results$Basic_lower,
    results$Basic_upper,
    1
  ),
  finite = TRUE
)

plot(
  results$Tumor_size,
  results$WR,
  type = "l",
  lwd = 2,
  ylim = ylim_plot,
  xlab = "Tumor size",
  ylab = "Win ratio",
  main = "Updated local IPCW: WR(tumor size), PS = 0"
)

lines(
  results$Tumor_size,
  results$Basic_lower,
  lty = 2
)

lines(
  results$Tumor_size,
  results$Basic_upper,
  lty = 2
)

abline(
  h = 1,
  lty = 3
)

rug(
  analysis_dat$x,
  ticksize = 0.02
)

legend(
  "topleft",
  legend = c(
    "Local IPCW",
    "Pointwise 95% CI",
    "WR = 1"
  ),
  lty = c(
    1,
    2,
    3
  ),
  lwd = c(
    2,
    1,
    1
  ),
  bty = "n"
)

############################################################
## 16. SAVE RESULTS
############################################################

out_dir <- dirname(
  file_path
)

prefix <- "C9633_local_IPCW_tsize_PS0_updated"

write.csv(
  results,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_curve.csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  summary_results,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_summary.csv"
    )
  ),
  row.names = FALSE
)

saveRDS(
  list(
    analysis = "WR(tumor size), PS = 0",
    tau = tau,
    x_grid = x_grid,
    h = h,
    bandwidth = bw,
    dat0 = dat0,
    dat1 = dat1,
    results = results,
    summary = summary_results,
    boot_logWR = boot_logWR,
    success = success,
    T_const_boot = T_const_boot,
    T_linear_boot = T_linear_boot
  ),
  file.path(
    out_dir,
    paste0(
      prefix,
      "_full_results.rds"
    )
  )
)

res_tsize_ps0 <- list(
  dat0 = dat0,
  dat1 = dat1,
  x_grid = x_grid,
  h = h,
  wr_result = results,
  test_summary = summary_results,
  bandwidth = bw,
  bootstrap_logWR = boot_logWR
)

cat(
  "\nAnalysis complete.\n"
)

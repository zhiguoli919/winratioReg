############################################################
## C9633 FINAL PH--BB1 ANALYSIS: WR(age), OVERALL
##
## Primary event: death
## Secondary event: recurrence/progression before death
## Horizon: tau = 5 years
##
## This script:
##   1. Reads the C9633 data.
##   2. Constructs the two right-censored margins at tau = 5.
##   3. Uses a natural-spline PH model (df = 3) for age.
##   4. Fits arm-specific PH--BB1 sieve models.
##   5. Evaluates WR(age) on 21 points from the empirical
##      5th to 95th percentiles.
##   6. Uses an arm-specific Exp(1) weighted multiplier bootstrap.
##   7. Reports pointwise CIs, a constant projection with CI,
##      a linear projection, and tests of constancy/log-linearity.
##
## IMPORTANT:
##   Keep PH_BB1_recovered_functions.R in the same folder as this
##   script (or edit functions_file below).
############################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(survival)
  library(splines)
})

############################################################
## 1. SETTINGS
############################################################

file_path <- Sys.getenv("C9633_DATA", unset = "C9633.csv.xlsx")

functions_file <- Sys.getenv("WINRATIOREG_BB1_HELPERS", unset = "PH_BB1_recovered_functions.R")

tau <- 5
analysis_var <- "age"
spline_df <- 3
m1 <- 6
m2 <- 6
grid_length <- 21
grid_lower_q <- 0.05
grid_upper_q <- 0.95
B_boot <- 1000

## Reproducible seed for the age analysis.
seed_boot <- 20261012

alpha <- 0.05

############################################################
## 2. LOAD RECOVERED PH--BB1 FUNCTIONS
############################################################

if (!file.exists(functions_file)) {
  stop(
    "Cannot find PH_BB1_recovered_functions.R at:\n",
    functions_file
  )
}

source(functions_file)

required_functions <- c(
  "make_ns_spec",
  "eval_ns_spec",
  "fit_bb1_ph_sieve",
  "logwr_grid_bb1",
  "weighted_bootstrap_bb1"
)

missing_functions <- required_functions[
  !vapply(required_functions, exists, logical(1), mode = "function")
]

if (length(missing_functions) > 0L) {
  stop(
    "Missing required functions after source(): ",
    paste(missing_functions, collapse = ", ")
  )
}

############################################################
## 3. READ C9633 AND CONSTRUCT ENDPOINTS
############################################################

dat_raw <- read_excel(file_path)

dat <- dat_raw %>%
  rename_with(tolower) %>%
  mutate(
    treat = as.numeric(treat),
    age = as.numeric(age),
    survtime = as.numeric(survtime),
    survcens = as.numeric(survcens),
    failtime = as.numeric(failtime),
    failcens = as.numeric(failcens)
  )

required_vars <- c(
  "treat", "age",
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

analysis_dat <- dat %>%
  mutate(
    death_time = survtime,
    death_event = as.integer(survcens == 1),

    recurrence_time_raw = case_when(
      failcens == 1 & !is.na(failtime) ~ failtime,
      failcens %in% c(0, 2) ~ survtime,
      TRUE ~ NA_real_
    ),

    recurrence_event = as.integer(
      failcens == 1 & !is.na(failtime)
    ),

    U1 = pmin(death_time, tau),
    Delta1 = as.integer(
      death_event == 1 &
        death_time <= tau
    ),

    U2 = pmin(recurrence_time_raw, tau),
    Delta2 = as.integer(
      recurrence_event == 1 &
        recurrence_time_raw <= tau
    ),

    z = as.numeric(age)
  ) %>%
  filter(
    treat %in% c(0, 1),
    is.finite(z),
    is.finite(U1),
    is.finite(U2),
    !is.na(Delta1),
    !is.na(Delta2)
  ) %>%
  select(
    treat,
    z,
    U1,
    Delta1,
    U2,
    Delta2
  )

if (nrow(analysis_dat) == 0L) {
  stop("No complete observations are available.")
}

############################################################
## 4. NATURAL-SPLINE COVARIATE SPECIFICATION
############################################################

z_all <- analysis_dat$z

spline_spec <- make_ns_spec(
  z_all,
  df = spline_df
)

X_all <- eval_ns_spec(
  z_all,
  spline_spec
)

dat0 <- analysis_dat %>%
  filter(treat == 0) %>%
  select(U1, Delta1, U2, Delta2)

dat1 <- analysis_dat %>%
  filter(treat == 1) %>%
  select(U1, Delta1, U2, Delta2)

idx0 <- which(analysis_dat$treat == 0)
idx1 <- which(analysis_dat$treat == 1)

attr(dat0, "X_mat") <- X_all[idx0, , drop = FALSE]
attr(dat1, "X_mat") <- X_all[idx1, , drop = FALSE]

attr(dat0, "tmax") <- tau
attr(dat1, "tmax") <- tau

if (nrow(dat0) == 0L || nrow(dat1) == 0L) {
  stop("Both treatment arms must contain observations.")
}

############################################################
## 5. COVARIATE GRID
############################################################

z_grid <- seq(
  as.numeric(
    quantile(z_all, probs = grid_lower_q, na.rm = TRUE)
  ),
  as.numeric(
    quantile(z_all, probs = grid_upper_q, na.rm = TRUE)
  ),
  length.out = grid_length
)

X_grid <- eval_ns_spec(
  z_grid,
  spline_spec
)

cat("\n============================================================\n")
cat("C9633 PH--BB1 ANALYSIS: WR(age), OVERALL\n")
cat("============================================================\n")
cat("Total N:", nrow(analysis_dat), "\n")
cat("Control N:", nrow(dat0), "\n")
cat("Treatment N:", nrow(dat1), "\n")
cat("tau =", tau, "\n")
cat("Spline df =", spline_df, "\n")
cat("Spline knots:", spline_spec$knots, "\n")
cat("Boundary knots:", spline_spec$Boundary.knots, "\n")
cat("Reference:", spline_spec$reference, "\n")
cat("Grid range:", range(z_grid), "\n")
cat("Bootstrap B =", B_boot, "\n")
cat("============================================================\n\n")

############################################################
## 6. FIT ARM-SPECIFIC PH--BB1 MODELS
############################################################

cat("Fitting control-arm PH--BB1 model...\n")

fit0 <- fit_bb1_ph_sieve(
  dat0,
  m1 = m1,
  m2 = m2,
  multistart = TRUE,
  maxit = 2500
)

cat("Control fit complete.\n\n")

cat("Fitting treatment-arm PH--BB1 model...\n")

fit1 <- fit_bb1_ph_sieve(
  dat1,
  m1 = m1,
  m2 = m2,
  multistart = TRUE,
  maxit = 2500
)

cat("Treatment fit complete.\n\n")

dependence_results <- data.frame(
  Arm = c("Control", "Treatment"),
  theta = c(fit0$theta, fit1$theta),
  delta = c(fit0$delta, fit1$delta),
  Kendall_tau = c(fit0$kendall_tau, fit1$kendall_tau),
  Lower_tail = c(fit0$lower_tail, fit1$lower_tail),
  Upper_tail = c(fit0$upper_tail, fit1$upper_tail),
  logLik = c(fit0$loglik, fit1$loglik),
  AIC = c(fit0$AIC, fit1$AIC)
)

print(dependence_results, row.names = FALSE)

bb1_total_loglik <- fit0$loglik + fit1$loglik
bb1_total_AIC <- fit0$AIC + fit1$AIC

cat("\nBB1 total log-likelihood =", bb1_total_loglik, "\n")
cat("BB1 total AIC =", bb1_total_AIC, "\n\n")

############################################################
## 7. ORIGINAL-SAMPLE WR CURVE
############################################################

logWR_hat <- logwr_grid_bb1(
  fit1 = fit1,
  fit0 = fit0,
  X_grid = X_grid,
  tau = tau
)

if (any(!is.finite(logWR_hat))) {
  stop("Original PH--BB1 log-WR curve contains invalid values.")
}

WR_hat <- exp(logWR_hat)

############################################################
## 8. FIXED-GRID PROJECTION FUNCTIONS
############################################################

trapz_weights <- function(x) {

  J <- length(x)

  if (J < 2L) {
    stop("At least two grid points are required.")
  }

  w <- numeric(J)
  w[1] <- (x[2] - x[1]) / 2
  w[J] <- (x[J] - x[J - 1]) / 2

  if (J > 2L) {
    for (j in 2:(J - 1L)) {
      w[j] <- (x[j + 1L] - x[j - 1L]) / 2
    }
  }

  w / sum(w)
}

project_constant <- function(logcurve, qweight) {

  alpha_hat <- sum(qweight * logcurve)
  fitted <- rep(alpha_hat, length(logcurve))

  list(
    alpha = alpha_hat,
    fitted = fitted,
    resid = logcurve - fitted
  )
}

project_linear <- function(logcurve, x_grid, qweight) {

  X <- cbind(1, x_grid)
  WX <- qweight * X

  beta <- solve(
    t(X) %*% WX,
    t(X) %*% (qweight * logcurve)
  )

  fitted <- as.vector(X %*% beta)

  list(
    beta = as.vector(beta),
    fitted = fitted,
    resid = logcurve - fitted
  )
}

qweight <- trapz_weights(z_grid)

const_obs <- project_constant(
  logWR_hat,
  qweight
)

linear_obs <- project_linear(
  logWR_hat,
  z_grid,
  qweight
)

############################################################
## 9. WEIGHTED MULTIPLIER BOOTSTRAP
############################################################

cat("\nStarting weighted PH--BB1 bootstrap...\n\n")

boot <- weighted_bootstrap_bb1(
  dat1 = dat1,
  dat0 = dat0,
  fit1 = fit1,
  fit0 = fit0,
  X_grid = X_grid,
  tau = tau,
  B = B_boot,
  m1 = m1,
  m2 = m2,
  seed = seed_boot,
  maxit = 1000
)

cat(
  "\nSuccessful bootstrap fits:",
  boot$n_success,
  "of",
  B_boot,
  "\n\n"
)

if (boot$n_success < 20L) {
  stop("Too few successful bootstrap fits.")
}

############################################################
## 10. CONSTANCY AND LOG-LINEARITY TESTS
############################################################

D_const_obs <- sum(
  qweight *
    const_obs$resid^2
)

D_linear_obs <- sum(
  qweight *
    linear_obs$resid^2
)

D_const_boot <- rep(NA_real_, B_boot)
D_linear_boot <- rep(NA_real_, B_boot)
alpha_boot <- rep(NA_real_, B_boot)

beta_boot <- matrix(
  NA_real_,
  nrow = B_boot,
  ncol = 2
)

for (b in which(boot$success)) {

  log_b <- boot$boot_logWR[b, ]

  diff_b <- log_b - logWR_hat

  const_diff <- project_constant(
    diff_b,
    qweight
  )

  linear_diff <- project_linear(
    diff_b,
    z_grid,
    qweight
  )

  D_const_boot[b] <- sum(
    qweight *
      const_diff$resid^2
  )

  D_linear_boot[b] <- sum(
    qweight *
      linear_diff$resid^2
  )

  alpha_boot[b] <-
    project_constant(
      log_b,
      qweight
    )$alpha

  beta_boot[b, ] <-
    project_linear(
      log_b,
      z_grid,
      qweight
    )$beta
}

p_constancy <- mean(
  D_const_boot[is.finite(D_const_boot)] >=
    D_const_obs
)

p_linearity <- mean(
  D_linear_boot[is.finite(D_linear_boot)] >=
    D_linear_obs
)

############################################################
## 11. CENTERED/BASIC CI FOR CONSTANT PROJECTION
############################################################

alpha_hat <- const_obs$alpha

alpha_diff <- alpha_boot - alpha_hat
alpha_diff <- alpha_diff[is.finite(alpha_diff)]

alpha_q <- quantile(
  alpha_diff,
  probs = c(alpha / 2, 1 - alpha / 2),
  na.rm = TRUE,
  names = FALSE
)

constant_logWR_CI <- c(
  alpha_hat - alpha_q[2],
  alpha_hat - alpha_q[1]
)

constant_WR_hat <- exp(alpha_hat)
constant_WR_CI <- exp(constant_logWR_CI)

############################################################
## 12. FINAL RESULTS
############################################################

results <- data.frame(
  Age = z_grid,
  logWR = logWR_hat,
  SE_logWR = boot$se,
  WR = WR_hat,
  Normal_lower = exp(boot$normal_lower),
  Normal_upper = exp(boot$normal_upper),
  Basic_lower = exp(boot$basic_lower),
  Basic_upper = exp(boot$basic_upper)
)

summary_results <- data.frame(
  Analysis = "PH-BB1 WR(age)",
  N = nrow(analysis_dat),
  N_control = nrow(dat0),
  N_treatment = nrow(dat1),
  tau = tau,
  spline_df = spline_df,
  Constant_logWR = alpha_hat,
  Constant_WR = constant_WR_hat,
  Constant_WR_lower = constant_WR_CI[1],
  Constant_WR_upper = constant_WR_CI[2],
  Linear_intercept = linear_obs$beta[1],
  Linear_slope = linear_obs$beta[2],
  P_constancy = p_constancy,
  P_linearity = p_linearity,
  Bootstrap_success = boot$n_success,
  BB1_total_logLik = bb1_total_loglik,
  BB1_total_AIC = bb1_total_AIC
)

cat("\n============================================================\n")
cat("FINAL PH--BB1 RESULTS: WR(age), OVERALL\n")
cat("============================================================\n\n")

print(results, row.names = FALSE)

cat("\nSummary:\n\n")

print(summary_results, row.names = FALSE)

############################################################
## 13. VALIDATION AGAINST EARLIER AGE BB1 POINT ESTIMATES
############################################################

validation_targets <- data.frame(
  Age = c(44.000, 59.925, 75.850),
  Expected_WR = c(
    2.1199598,
    0.8136311,
    3.5368080
  )
)

validation_results <- merge(
  validation_targets,
  results[, c("Age", "WR")],
  by = "Age",
  all.x = TRUE
)

validation_results$Difference <-
  validation_results$WR -
  validation_results$Expected_WR

cat(
  "\nValidation against the earlier correct BB1 point estimates:\n\n"
)

print(validation_results, row.names = FALSE)

############################################################
## 14. PLOT
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
  results$Age,
  results$WR,
  type = "l",
  lwd = 2,
  ylim = ylim_plot,
  xlab = "Age",
  ylab = "Win ratio",
  main = "PH--BB1 model-based WR(age) at 5 years"
)

lines(
  results$Age,
  results$Basic_lower,
  lty = 2
)

lines(
  results$Age,
  results$Basic_upper,
  lty = 2
)

abline(h = 1, lty = 3)

rug(
  z_all,
  ticksize = 0.02
)

legend(
  "topleft",
  legend = c(
    "PH--BB1",
    "Pointwise 95% CI",
    "WR = 1"
  ),
  lty = c(1, 2, 3),
  lwd = c(2, 1, 1),
  bty = "n"
)

############################################################
## 15. SAVE RESULTS
############################################################

out_dir <- dirname(file_path)

prefix <- "C9633_PH_BB1_age_overall_B1000"

write.csv(
  results,
  file.path(
    out_dir,
    paste0(prefix, "_curve.csv")
  ),
  row.names = FALSE
)

write.csv(
  summary_results,
  file.path(
    out_dir,
    paste0(prefix, "_summary.csv")
  ),
  row.names = FALSE
)

write.csv(
  dependence_results,
  file.path(
    out_dir,
    paste0(prefix, "_dependence.csv")
  ),
  row.names = FALSE
)

saveRDS(
  list(
    analysis = "PH-BB1 WR(age), overall",
    tau = tau,
    spline_df = spline_df,
    spline_spec = spline_spec,
    z_grid = z_grid,
    X_grid = X_grid,
    dat0 = dat0,
    dat1 = dat1,
    fit0 = fit0,
    fit1 = fit1,
    dependence_results = dependence_results,
    bb1_total_loglik = bb1_total_loglik,
    bb1_total_AIC = bb1_total_AIC,
    results = results,
    summary = summary_results,
    bootstrap = boot,
    D_const_boot = D_const_boot,
    D_linear_boot = D_linear_boot,
    alpha_boot = alpha_boot,
    beta_boot = beta_boot,
    validation = validation_results
  ),
  file.path(
    out_dir,
    paste0(prefix, "_full_results.rds")
  )
)

cat("\nAnalysis complete.\n")

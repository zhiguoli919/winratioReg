############################################################
## PH--COPULA STAGE-2 VALIDATION
##
## Compares winratioReg 0.2.0 against the original
## BB1 / Clayton / Gumbel implementation used during
## manuscript development.
##
## Run from the package root after devtools::load_all():
##   source("inst/validation/phcopula_validation.R")
############################################################

cat("\n============================================================\n")
cat("PH--COPULA PACKAGE VALIDATION\n")
cat("============================================================\n")

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("Package 'survival' is required for this validation.")
}
suppressPackageStartupMessages(library(survival))

ref_file <- file.path(
  "inst", "validation", "reference", "PH_BB1_recovered_functions.R"
)
if (!file.exists(ref_file)) stop("Reference implementation not found: ", ref_file)

ref <- new.env(parent = globalenv())
sys.source(ref_file, envir = ref)

############################################################
## FIXED DATASET
############################################################

set.seed(20261020)
n0 <- 70L
n1 <- 70L
n <- n0 + n1
A <- sample(c(rep(0L, n0), rep(1L, n1)))
x <- rnorm(n)

## Dependence is generated only to create a realistic fixed dataset;
## the same dataset is fitted under all three working copulas.
z_common <- rnorm(n)
z1 <- 0.65 * z_common + sqrt(1 - 0.65^2) * rnorm(n)
z2 <- 0.65 * z_common + sqrt(1 - 0.65^2) * rnorm(n)
u1 <- pmin(pmax(pnorm(z1), 1e-8), 1 - 1e-8)
u2 <- pmin(pmax(pnorm(z2), 1e-8), 1 - 1e-8)

eta1 <- -0.25 * A + 0.25 * x - 0.20 * A * x
eta2 <- -0.18 * A + 0.15 * x - 0.12 * A * x
T1 <- qexp(u1, rate = 0.10 * exp(eta1))
T2 <- qexp(u2, rate = 0.18 * exp(eta2))
C <- pmin(rexp(n, rate = 0.06 * exp(0.10 * x)), 8)
U1 <- pmin(T1, C)
Delta1 <- as.integer(T1 <= C)
U2 <- pmin(T2, T1, C)
Delta2 <- as.integer(T2 <= pmin(T1, C))

dat <- data.frame(
  treat = A, x = x,
  U1 = U1, Delta1 = Delta1,
  U2 = U2, Delta2 = Delta2
)

grid <- seq(-0.75, 0.75, length.out = 5)
tau <- 5
upper <- 8
m1 <- 3L
m2 <- 3L
maxit <- 1200L

make_ref_arm <- function(d) {
  out <- d[, c("U1", "Delta1", "U2", "Delta2"), drop = FALSE]
  attr(out, "X_mat") <- matrix(d$x, ncol = 1L)
  attr(out, "tmax") <- upper
  out
}
ref0 <- make_ref_arm(dat[dat$treat == 0L, , drop = FALSE])
ref1 <- make_ref_arm(dat[dat$treat == 1L, , drop = FALSE])
X_grid_ref <- matrix(grid, ncol = 1L)

############################################################
## FIXED-FORMULA COMPONENT CHECKS
############################################################

u_chk <- c(0.20, 0.45, 0.80)
v_chk <- c(0.35, 0.60, 0.75)

component_diff <- function(a, b) {
  max(c(
    abs(a$C - b$C),
    abs(a$Cu - b$Cu),
    abs(a$Cv - b$Cv),
    abs(a$density - b$density)
  ))
}

comp_bb1 <- component_diff(
  winratioReg:::.ph_bb1_components(u_chk, v_chk, 0.7, 1.8),
  ref$bb1_components(u_chk, v_chk, 0.7, 1.8)
)
comp_clayton <- component_diff(
  winratioReg:::.ph_clayton_components(u_chk, v_chk, 1.2),
  ref$clayton_components(u_chk, v_chk, 1.2)
)
comp_gumbel <- component_diff(
  winratioReg:::.ph_gumbel_components(u_chk, v_chk, 2.0),
  ref$gumbel_components(u_chk, v_chk, 2.0)
)

############################################################
## FIT AND WR COMPARISONS
############################################################

fit_one_family <- function(family) {
  cat("\nValidating", toupper(family), "...\n")

  pkg <- winratio_phcopula(
    data = dat,
    arm = "treat",
    time1 = "U1", status1 = "Delta1",
    time2 = "U2", status2 = "Delta2",
    formula = ~ x,
    tau = tau,
    copula = family,
    degree = c(m1, m2),
    sieve_upper = upper,
    grid = grid,
    B = 0L,
    maxit = maxit,
    integration_reltol = 1e-7,
    integration_subdivisions = 300L,
    verbose = FALSE
  )

  if (family == "bb1") {
    r1 <- ref$fit_bb1_ph_sieve(ref1, m1 = m1, m2 = m2, maxit = maxit)
    r0 <- ref$fit_bb1_ph_sieve(ref0, m1 = m1, m2 = m2, maxit = maxit)
    rlog <- ref$logwr_grid_bb1(r1, r0, X_grid_ref, tau)
    delta_diff <- max(abs(c(pkg$fit1$delta - r1$delta, pkg$fit0$delta - r0$delta)))
  } else {
    r1 <- ref$fit_phcopula_sieve(ref1, family = family, m1 = m1, m2 = m2, maxit = maxit)
    r0 <- ref$fit_phcopula_sieve(ref0, family = family, m1 = m1, m2 = m2, maxit = maxit)
    rlog <- ref$logwr_grid(r1, r0, X_grid_ref, tau)
    delta_diff <- NA_real_
  }

  list(
    family = family,
    loglik_diff = abs(pkg$loglik - (r1$loglik + r0$loglik)),
    theta_diff = max(abs(c(pkg$fit1$theta - r1$theta, pkg$fit0$theta - r0$theta))),
    delta_diff = delta_diff,
    logwr_diff = max(abs(pkg$estimates$logWR - rlog)),
    pkg = pkg,
    ref1 = r1,
    ref0 = r0
  )
}

fits <- lapply(c("bb1", "clayton", "gumbel"), fit_one_family)

rows <- list(
  data.frame(check = "BB1 copula components", max_abs_difference = comp_bb1, tolerance = 1e-10),
  data.frame(check = "Clayton copula components", max_abs_difference = comp_clayton, tolerance = 1e-10),
  data.frame(check = "Gumbel copula components", max_abs_difference = comp_gumbel, tolerance = 1e-10)
)

for (z in fits) {
  fam <- toupper(z$family)
  rows[[length(rows) + 1L]] <- data.frame(
    check = paste(fam, "total logLik"),
    max_abs_difference = z$loglik_diff,
    tolerance = 1e-6
  )
  rows[[length(rows) + 1L]] <- data.frame(
    check = paste(fam, "theta"),
    max_abs_difference = z$theta_diff,
    tolerance = 1e-5
  )
  if (z$family == "bb1") {
    rows[[length(rows) + 1L]] <- data.frame(
      check = "BB1 delta",
      max_abs_difference = z$delta_diff,
      tolerance = 1e-5
    )
  }
  rows[[length(rows) + 1L]] <- data.frame(
    check = paste(fam, "log-WR grid"),
    max_abs_difference = z$logwr_diff,
    tolerance = 1e-5
  )
}

res <- do.call(rbind, rows)
res$pass <- is.finite(res$max_abs_difference) &
  res$max_abs_difference <= res$tolerance
rownames(res) <- NULL

cat("\n============================================================\n")
cat("VALIDATION RESULTS\n")
cat("============================================================\n")
print(res, row.names = FALSE)

out_file <- file.path("inst", "validation", "phcopula_validation_results.csv")
write.csv(res, out_file, row.names = FALSE)
cat("\nSaved comparison table to:\n", normalizePath(out_file, winslash = "/", mustWork = FALSE), "\n", sep = "")

cat("\n============================================================\n")
if (all(res$pass)) {
  cat("PASS: PH--copula package reproduces the original BB1/Clayton/Gumbel\n")
  cat("      implementation within the stated numerical tolerances.\n")
} else {
  cat("FAIL: At least one PH--copula validation comparison exceeded tolerance.\n")
  cat("      Please send the table above before changing any code.\n")
}
cat("============================================================\n")

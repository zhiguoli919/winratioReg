############################################################
## winratioReg: Frank / Gaussian PH--copula validation
##
## Purpose
## -------
## Independently validate the two Stage-3 copula additions against
## the CRAN 'copula' package, then run a fixed-data smoke fit through
## the full winratio_phcopula() interface.
##
## Run from the package root after installing the suggested package:
##   install.packages(c("copula", "pbivnorm"))
##   devtools::load_all()
##   source("inst/validation/frank_gaussian_validation.R")
############################################################

cat("\n============================================================\n")
cat("FRANK / GAUSSIAN PH--COPULA VALIDATION\n")
cat("============================================================\n")

needed <- c("copula", "pbivnorm", "survival")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) {
  stop("Please install the following package(s) first: ",
       paste(missing, collapse = ", "), call. = FALSE)
}

############################################################
## 1. INDEPENDENT COPULA-COMPONENT CHECKS
############################################################

u <- c(0.18, 0.31, 0.52, 0.77)
v <- c(0.27, 0.63, 0.44, 0.81)
h <- 1e-6

ref_components <- function(family, par) {
  cp <- if (family == "frank") {
    copula::frankCopula(param = par, dim = 2)
  } else {
    copula::normalCopula(param = par, dim = 2, dispstr = "un")
  }

  C <- copula::pCopula(cbind(u, v), cp)
  density <- copula::dCopula(cbind(u, v), cp)

  ## Independent numerical derivatives of the reference copula CDF.
  Cu <- (
    copula::pCopula(cbind(u + h, v), cp) -
      copula::pCopula(cbind(u - h, v), cp)
  ) / (2 * h)
  Cv <- (
    copula::pCopula(cbind(u, v + h), cp) -
      copula::pCopula(cbind(u, v - h), cp)
  ) / (2 * h)

  list(C = C, Cu = Cu, Cv = Cv, density = density, copula = cp)
}

pkg_frank <- winratioReg:::.ph_frank_components(u, v, theta = 4.0)
ref_frank <- ref_components("frank", 4.0)

pkg_gaussian <- winratioReg:::.ph_gaussian_components(u, v, rho = 0.55)
ref_gaussian <- ref_components("gaussian", 0.55)

rows <- list(
  data.frame(
    check = "Frank C",
    max_abs_difference = max(abs(pkg_frank$C - ref_frank$C)),
    tolerance = 1e-10
  ),
  data.frame(
    check = "Frank partial derivatives",
    max_abs_difference = max(abs(c(pkg_frank$Cu - ref_frank$Cu,
                                   pkg_frank$Cv - ref_frank$Cv))),
    tolerance = 5e-6
  ),
  data.frame(
    check = "Frank density",
    max_abs_difference = max(abs(pkg_frank$density - ref_frank$density)),
    tolerance = 1e-8
  ),
  data.frame(
    check = "Gaussian C",
    max_abs_difference = max(abs(pkg_gaussian$C - ref_gaussian$C)),
    tolerance = 2e-8
  ),
  data.frame(
    check = "Gaussian partial derivatives",
    max_abs_difference = max(abs(c(pkg_gaussian$Cu - ref_gaussian$Cu,
                                   pkg_gaussian$Cv - ref_gaussian$Cv))),
    tolerance = 5e-6
  ),
  data.frame(
    check = "Gaussian density",
    max_abs_difference = max(abs(pkg_gaussian$density - ref_gaussian$density)),
    tolerance = 1e-8
  )
)

## Independence limits.
ind_f <- winratioReg:::.ph_frank_components(u, v, theta = 0)
ind_g <- winratioReg:::.ph_gaussian_components(u, v, rho = 0)
ind_truth <- list(C = u * v, Cu = v, Cv = u, density = rep(1, length(u)))
ind_diff <- function(z) max(abs(c(
  z$C - ind_truth$C,
  z$Cu - ind_truth$Cu,
  z$Cv - ind_truth$Cv,
  z$density - ind_truth$density
)))
rows[[length(rows) + 1L]] <- data.frame(
  check = "Frank independence limit",
  max_abs_difference = ind_diff(ind_f),
  tolerance = 1e-10
)
rows[[length(rows) + 1L]] <- data.frame(
  check = "Gaussian independence limit",
  max_abs_difference = ind_diff(ind_g),
  tolerance = 1e-10
)

## Kendall's tau against the independent copula package implementation.
frank_tau_pkg <- winratioReg:::.ph_dependence_summary(
  "frank", list(theta = 4.0)
)$kendall_tau
frank_tau_ref <- as.numeric(copula::tau(copula::frankCopula(4.0, dim = 2)))

gauss_tau_pkg <- winratioReg:::.ph_dependence_summary(
  "gaussian", list(theta = 0.55)
)$kendall_tau
gauss_tau_ref <- as.numeric(copula::tau(copula::normalCopula(0.55, dim = 2)))

rows[[length(rows) + 1L]] <- data.frame(
  check = "Frank Kendall tau",
  max_abs_difference = abs(frank_tau_pkg - frank_tau_ref),
  tolerance = 1e-8
)
rows[[length(rows) + 1L]] <- data.frame(
  check = "Gaussian Kendall tau",
  max_abs_difference = abs(gauss_tau_pkg - gauss_tau_ref),
  tolerance = 1e-10
)

############################################################
## 2. FULL-INTERFACE FIXED-DATA SMOKE FITS
############################################################

simulate_survival_copula_data <- function(family, par, seed) {
  set.seed(seed)
  n0 <- 80L
  n1 <- 80L
  n <- n0 + n1
  A <- sample(c(rep(0L, n0), rep(1L, n1)))
  x <- stats::rnorm(n)

  cp <- if (family == "frank") {
    copula::frankCopula(param = par, dim = 2)
  } else {
    copula::normalCopula(param = par, dim = 2, dispstr = "un")
  }
  UV <- copula::rCopula(n, cp)

  eta1 <- -0.25 * A + 0.20 * x - 0.15 * A * x
  eta2 <- -0.18 * A + 0.12 * x - 0.10 * A * x
  rate1 <- 0.10 * exp(eta1)
  rate2 <- 0.18 * exp(eta2)

  ## U_j = S_j(T_j), so -log(U_j)/rate_j generates the selected
  ## copula directly on the survival scale used by the fitted model.
  T1 <- -log(UV[, 1L]) / rate1
  T2 <- -log(UV[, 2L]) / rate2
  C <- pmin(stats::rexp(n, rate = 0.05), 8)

  data.frame(
    treat = A,
    x = x,
    U1 = pmin(T1, C),
    Delta1 = as.integer(T1 <= C),
    U2 = pmin(T2, T1, C),
    Delta2 = as.integer(T2 <= pmin(T1, C))
  )
}

fit_new_family <- function(family, par, seed) {
  cat("\nFitting", toupper(family), "smoke-validation model ...\n")
  dat <- simulate_survival_copula_data(family, par, seed)
  fit <- winratio_phcopula(
    data = dat,
    arm = "treat",
    time1 = "U1", status1 = "Delta1",
    time2 = "U2", status2 = "Delta2",
    formula = ~ x,
    tau = 5,
    copula = family,
    degree = 3,
    sieve_upper = 8,
    grid = seq(-0.75, 0.75, length.out = 5),
    B = 0L,
    multistart = TRUE,
    maxit = 1200L,
    integration_reltol = 1e-7,
    integration_subdivisions = 300L,
    verbose = FALSE
  )
  list(
    family = family,
    fit = fit,
    ok = fit$fit0$convergence == 0L && fit$fit1$convergence == 0L &&
      is.finite(fit$loglik) && all(is.finite(fit$estimates$logWR)) &&
      all(is.finite(fit$estimates$WR)) && all(fit$estimates$WR > 0)
  )
}

smoke_frank <- fit_new_family("frank", 4.0, 20261031L)
smoke_gaussian <- fit_new_family("gaussian", 0.55, 20261101L)

rows[[length(rows) + 1L]] <- data.frame(
  check = "Frank full fit smoke test",
  max_abs_difference = if (smoke_frank$ok) 0 else Inf,
  tolerance = 0
)
rows[[length(rows) + 1L]] <- data.frame(
  check = "Gaussian full fit smoke test",
  max_abs_difference = if (smoke_gaussian$ok) 0 else Inf,
  tolerance = 0
)

############################################################
## 3. REPORT
############################################################

res <- do.call(rbind, rows)
res$pass <- is.finite(res$max_abs_difference) &
  res$max_abs_difference <= res$tolerance
rownames(res) <- NULL

cat("\n============================================================\n")
cat("VALIDATION RESULTS\n")
cat("============================================================\n")
print(res, row.names = FALSE)

cat("\nFitted dependence parameters from smoke fits:\n")
cat("  Frank control/treatment theta:",
    signif(smoke_frank$fit$fit0$theta, 6), "/",
    signif(smoke_frank$fit$fit1$theta, 6), "\n")
cat("  Gaussian control/treatment rho:",
    signif(smoke_gaussian$fit$fit0$theta, 6), "/",
    signif(smoke_gaussian$fit$fit1$theta, 6), "\n")

out_file <- file.path(
  "inst", "validation", "frank_gaussian_validation_results.csv"
)
write.csv(res, out_file, row.names = FALSE)
cat("\nSaved comparison table to:\n",
    normalizePath(out_file, winslash = "/", mustWork = FALSE), "\n", sep = "")

cat("\n============================================================\n")
if (all(res$pass)) {
  cat("PASS: Frank and Gaussian copula components agree with the independent\n")
  cat("      'copula' reference and both full PH--copula smoke fits succeeded.\n")
} else {
  cat("FAIL: At least one Frank/Gaussian validation check exceeded tolerance.\n")
  cat("      Please send the table above before changing any code.\n")
}
cat("============================================================\n")

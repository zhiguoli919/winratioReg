############################################################
## winratioReg: Bernstein-degree selection validation
##
## Purpose
## -------
## Validate degree = "auto" against explicit fits over the same
## candidate set. Automatic selection uses one common degree for both
## margins and minimizes total BIC under the specified copula family.
## The selected degree is fixed for any subsequent multiplier bootstrap.
##
## Run from the package root:
##   devtools::load_all()
##   source("inst/validation/degree_selection_validation.R")
############################################################

cat("\n============================================================\n")
cat("BERNSTEIN-DEGREE SELECTION VALIDATION\n")
cat("============================================================\n")

if (!requireNamespace("copula", quietly = TRUE)) {
  stop("Please install the suggested package 'copula' first.", call. = FALSE)
}

simulate_degree_data <- function(seed = 20261102L) {
  set.seed(seed)
  n0 <- 70L
  n1 <- 70L
  n <- n0 + n1
  A <- sample(c(rep(0L, n0), rep(1L, n1)))
  x <- stats::rnorm(n)

  cp <- copula::claytonCopula(param = 1.2, dim = 2)
  UV <- copula::rCopula(n, cp)

  eta1 <- -0.28 * A + 0.22 * x - 0.12 * A * x
  eta2 <- -0.20 * A + 0.14 * x - 0.08 * A * x
  rate1 <- 0.11 * exp(eta1)
  rate2 <- 0.19 * exp(eta2)

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

dat <- simulate_degree_data()
candidates <- 3:4
common_args <- list(
  data = dat,
  arm = "treat",
  time1 = "U1", status1 = "Delta1",
  time2 = "U2", status2 = "Delta2",
  formula = ~ x,
  tau = 5,
  copula = "clayton",
  sieve_upper = 8,
  grid = seq(-0.75, 0.75, length.out = 5),
  B = 0L,
  multistart = TRUE,
  maxit = 1000L,
  integration_reltol = 1e-7,
  integration_subdivisions = 250L,
  verbose = FALSE
)

cat("\nFitting automatic degree selector ...\n")
fit_auto <- do.call(
  winratio_phcopula,
  c(common_args, list(degree = "auto", degree_candidates = candidates))
)

auto_tab <- fit_auto$degree_selection$table

cat("\nFitting the same candidate degrees explicitly ...\n")
manual <- lapply(candidates, function(m) {
  z <- do.call(winratio_phcopula, c(common_args, list(degree = m)))
  data.frame(
    degree = m,
    logLik = z$loglik,
    npar = z$npar,
    AIC = z$AIC,
    BIC = z$BIC,
    stringsAsFactors = FALSE
  )
})
manual_tab <- do.call(rbind, manual)

idx <- match(candidates, auto_tab$degree)
rows <- list(
  data.frame(
    check = "candidate logLik values",
    max_abs_difference = max(abs(auto_tab$logLik[idx] - manual_tab$logLik)),
    tolerance = 1e-6
  ),
  data.frame(
    check = "candidate AIC values",
    max_abs_difference = max(abs(auto_tab$AIC[idx] - manual_tab$AIC)),
    tolerance = 1e-6
  ),
  data.frame(
    check = "candidate BIC values",
    max_abs_difference = max(abs(auto_tab$BIC[idx] - manual_tab$BIC)),
    tolerance = 1e-6
  ),
  data.frame(
    check = "selected degree equals minimum-BIC manual fit",
    max_abs_difference = abs(
      fit_auto$degree[1L] - manual_tab$degree[which.min(manual_tab$BIC)]
    ),
    tolerance = 0
  ),
  data.frame(
    check = "selected degree common to both margins",
    max_abs_difference = abs(fit_auto$degree[1L] - fit_auto$degree[2L]),
    tolerance = 0
  )
)

res <- do.call(rbind, rows)
res$pass <- is.finite(res$max_abs_difference) &
  res$max_abs_difference <= res$tolerance
rownames(res) <- NULL

cat("\nAutomatic-selection table:\n")
print(auto_tab, row.names = FALSE)
cat("\nManual-fit table:\n")
print(manual_tab, row.names = FALSE)
cat("\n============================================================\n")
cat("VALIDATION RESULTS\n")
cat("============================================================\n")
print(res, row.names = FALSE)

out_file <- file.path(
  "inst", "validation", "degree_selection_validation_results.csv"
)
write.csv(res, out_file, row.names = FALSE)
cat("\nSaved comparison table to:\n",
    normalizePath(out_file, winslash = "/", mustWork = FALSE), "\n", sep = "")

cat("\n============================================================\n")
if (all(res$pass)) {
  cat("PASS: Automatic Bernstein-degree selection reproduces explicit\n")
  cat("      candidate fits and selects the minimum-BIC common degree.\n")
} else {
  cat("FAIL: At least one degree-selection validation check exceeded tolerance.\n")
}
cat("============================================================\n")

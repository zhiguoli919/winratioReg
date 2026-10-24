test_that("selective Gaussian likelihood equals full-component likelihood", {
  data("winratio_example", package = "winratioReg", envir = environment())
  prep <- winratioReg:::.ph_prepare_data(
    winratio_example, "treat", "U1", "Delta1", "U2", "Delta2",
    ~ x1, 1, 0, 4, NULL, 2L, 2L
  )
  dat <- prep$dat1
  par <- winratioReg:::.ph_make_start(
    dat, 2L, 2L, "gaussian", list(theta = 0.2)
  )

  ll_fast <- winratioReg:::.ph_individual_loglik(par, dat, 2L, 2L, "gaussian")

  X <- as.matrix(attr(dat, "X_mat"))
  pars <- winratioReg:::.ph_unpack(par, ncol(X), 2L, 2L, "gaussian")
  lp1 <- as.vector(X %*% pars$beta1)
  lp2 <- as.vector(X %*% pars$beta2)
  b1 <- winratioReg:::.ph_baseline_eval_observed(dat, pars$phi1, 1L)
  b2 <- winratioReg:::.ph_baseline_eval_observed(dat, pars$phi2, 2L)
  logS1 <- -exp(lp1) * b1$Lambda
  logS2 <- -exp(lp2) * b2$Lambda
  S1 <- exp(pmax(logS1, -700))
  S2 <- exp(pmax(logS2, -700))
  logf1 <- log(b1$lambda) + lp1 + logS1
  logf2 <- log(b2$lambda) + lp2 + logS2
  cc <- winratioReg:::.ph_gaussian_components(S1, S2, pars$dependence$theta)

  d1 <- dat$Delta1
  d2 <- dat$Delta2
  both_events <- d1 == 1L & d2 == 1L
  only1 <- d1 == 1L & d2 == 0L
  only2 <- d1 == 0L & d2 == 1L
  both_cens <- d1 == 0L & d2 == 0L
  ll_full <- numeric(nrow(dat))
  ll_full[both_events] <- logf1[both_events] + logf2[both_events] + log(cc$density[both_events])
  ll_full[only1] <- logf1[only1] + log(cc$Cu[only1])
  ll_full[only2] <- logf2[only2] + log(cc$Cv[only2])
  ll_full[both_cens] <- log(cc$C[both_cens])

  expect_equal(ll_fast, ll_full, tolerance = 1e-12)
})

test_that("vectorized Simpson WR agrees with adaptive integration", {
  make_fit <- function(beta_shift, dep) {
    list(
      family = "gaussian", p = 1L,
      beta1 = 0.20 + beta_shift,
      beta2 = -0.10 + beta_shift,
      phi1 = c(0, cumsum(rep(0.18, 6L))),
      phi2 = c(0, cumsum(rep(0.22, 6L))),
      lower = 0, upper = 5,
      dependence = list(theta = dep)
    )
  }
  fit1 <- make_fit(-0.10, 0.35)
  fit0 <- make_fit(0.05, 0.25)
  Xg <- matrix(seq(-0.8, 0.8, length.out = 7L), ncol = 1L)

  a <- winratioReg:::.ph_logwr_grid(
    fit1, fit0, Xg, 4, method = "adaptive",
    rel.tol = 1e-9, subdivisions = 500L
  )
  s <- winratioReg:::.ph_logwr_grid(
    fit1, fit0, Xg, 4, method = "simpson", n_grid = 401L
  )
  expect_lt(max(abs(s - a)), 1e-4)
})

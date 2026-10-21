test_that("kernel weights are normalized", {
  X <- cbind(seq(-1, 1, length.out = 10))
  w <- winratioReg:::.local_kernel_weights(X, 0, 0.5)
  expect_equal(sum(w), 1, tolerance = 1e-12)
  expect_true(all(w >= 0))
})

test_that("Local IPCW produces finite one-dimensional estimates", {
  dat <- make_local_test_data()
  fit <- winratio_local(
    data = dat, arm = "arm", time1 = "U1", status1 = "Delta1",
    time2 = "U2", status2 = "Delta2", covariates = "x",
    tau = 2.2, grid = seq(-0.75, 0.75, length.out = 5),
    bandwidth = 0.65, B = 0
  )
  expect_s3_class(fit, "winratio_local")
  expect_true(all(is.finite(fit$estimates$WR)))
  expect_true(all(fit$estimates$WR > 0))
  expect_equal(length(fit$bandwidth), 1L)
})


test_that("Local IPCW supports two continuous covariates", {
  dat <- make_local_test_data(n_arm = 80L)
  dat$x2 <- rep(sin(seq(0, 2 * pi, length.out = 80L)), 2L)
  g <- expand.grid(x = c(-0.5, 0, 0.5), x2 = c(-0.5, 0, 0.5))
  fit <- winratio_local(
    data = dat, arm = "arm", time1 = "U1", status1 = "Delta1",
    time2 = "U2", status2 = "Delta2", covariates = c("x", "x2"),
    tau = 2.2, grid = g, bandwidth = c(0.7, 0.7), B = 0
  )
  expect_s3_class(fit, "winratio_local")
  expect_equal(fit$dimension, 2L)
  expect_true(all(is.finite(fit$estimates$WR)))
})

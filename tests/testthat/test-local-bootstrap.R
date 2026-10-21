test_that("bootstrap projections and tests are returned", {
  dat <- make_local_test_data()
  fit <- winratio_local(
    data = dat, arm = "arm", time1 = "U1", status1 = "Delta1",
    time2 = "U2", status2 = "Delta2", covariates = "x",
    tau = 2.2, grid = seq(-0.75, 0.75, length.out = 5),
    bandwidth = 0.65, B = 10, seed = 123
  )
  z <- blp(fit)
  expect_equal(nrow(z), 2L)
  expect_true(all(is.finite(z$estimate)))
  tt <- test_loglinearity(fit)
  expect_s3_class(tt, "htest")
  expect_true(is.na(tt$p.value) || (tt$p.value >= 0 && tt$p.value <= 1))
})

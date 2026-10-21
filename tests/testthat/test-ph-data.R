test_that("PH basis construction handles one-sided formulas", {
  d <- data.frame(x = seq(-1, 1, length.out = 20), z = seq(0, 2, length.out = 20))
  s1 <- winratioReg:::.ph_build_basis_spec(d, ~ x)
  expect_equal(ncol(s1$X), 1L)
  expect_equal(as.numeric(s1$X[, 1]), d$x)

  s2 <- winratioReg:::.ph_build_basis_spec(d, ~ x + z + I(x^2))
  expect_true(ncol(s2$X) >= 3L)
  nd <- data.frame(x = c(-0.5, 0.5), z = c(0.5, 1.5))
  Xn <- winratioReg:::.ph_basis_new(s2, nd)
  expect_equal(nrow(Xn), 2L)
  expect_equal(ncol(Xn), ncol(s2$X))
})

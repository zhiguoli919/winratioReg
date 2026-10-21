test_that("PH copula components are finite and valid", {
  u <- c(0.2, 0.5, 0.8)
  v <- c(0.3, 0.6, 0.7)

  b <- winratioReg:::.ph_bb1_components(u, v, 0.7, 1.8)
  c <- winratioReg:::.ph_clayton_components(u, v, 1.2)
  f <- winratioReg:::.ph_frank_components(u, v, 4.0)
  n <- winratioReg:::.ph_gaussian_components(u, v, 0.5)
  g <- winratioReg:::.ph_gumbel_components(u, v, 2.0)

  for (z in list(b, c, f, n, g)) {
    expect_true(all(is.finite(z$C)))
    expect_true(all(is.finite(z$Cu)))
    expect_true(all(is.finite(z$Cv)))
    expect_true(all(is.finite(z$density)))
    expect_true(all(z$C > 0 & z$C <= 1))
    expect_true(all(z$Cu > 0 & z$Cu <= 1))
    expect_true(all(z$Cv > 0 & z$Cv <= 1))
    expect_true(all(z$density > 0))
  }
})

test_that("Bernstein baseline is nonnegative", {
  phi <- c(0, 0.2, 0.5, 1.0)
  z <- winratioReg:::.ph_baseline_eval(seq(0, 5, length.out = 11), phi, 0, 5)
  expect_true(all(z$Lambda >= 0))
  expect_true(all(z$lambda > 0))
  expect_true(all(diff(z$Lambda) >= -1e-12))
})


test_that("Frank and Gaussian independence limits are correct", {
  u <- c(0.2, 0.5, 0.8)
  v <- c(0.3, 0.6, 0.7)

  f <- winratioReg:::.ph_frank_components(u, v, 0)
  g <- winratioReg:::.ph_gaussian_components(u, v, 0)

  for (z in list(f, g)) {
    expect_equal(z$C, u * v, tolerance = 1e-10)
    expect_equal(z$Cu, v, tolerance = 1e-10)
    expect_equal(z$Cv, u, tolerance = 1e-10)
    expect_equal(z$density, rep(1, length(u)), tolerance = 1e-10)
  }
})

test_that("new dependence transforms round trip", {
  for (th in c(-5, -1, 1, 5)) {
    z <- winratioReg:::.ph_dep_unconstrained("frank", list(theta = th))
    back <- winratioReg:::.ph_dep_from_unconstrained("frank", z)
    expect_equal(back$theta, th, tolerance = 1e-12)
  }
  for (rho in c(-0.8, -0.2, 0, 0.2, 0.8)) {
    z <- winratioReg:::.ph_dep_unconstrained("gaussian", list(theta = rho))
    back <- winratioReg:::.ph_dep_from_unconstrained("gaussian", z)
    expect_equal(back$theta, rho, tolerance = 1e-12)
  }
})

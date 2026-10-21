test_that("Bernstein degree specifications are parsed correctly", {
  expect_true(winratioReg:::.ph_is_auto_degree("auto"))
  expect_true(winratioReg:::.ph_is_auto_degree("AUTO"))
  expect_false(winratioReg:::.ph_is_auto_degree(6))

  expect_equal(winratioReg:::.ph_parse_fixed_degree(6), c(6L, 6L))
  expect_equal(winratioReg:::.ph_parse_fixed_degree(c(4, 6)), c(4L, 6L))
  expect_equal(winratioReg:::.ph_validate_degree_candidates(c(6, 4, 5, 4)),
               c(4L, 5L, 6L))

  expect_error(winratioReg:::.ph_parse_fixed_degree(c(0, 4)),
               "positive integers")
  expect_error(winratioReg:::.ph_validate_degree_candidates(c(3, 4.5)),
               "positive integers")
})

test_that("copula comparison requires fixed Bernstein degrees", {
  expect_error(
    compare_copulas(
      data = data.frame(), arm = "arm", time1 = "U1", status1 = "D1",
      time2 = "U2", status2 = "D2", formula = ~ x, tau = 5,
      degree = "auto"
    ),
    "requires an explicit Bernstein degree"
  )
})

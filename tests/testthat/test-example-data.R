test_that("synthetic example data obey the terminal-event observation scheme", {
  data("winratio_example", package = "winratioReg", envir = environment())
  expect_s3_class(winratio_example, "data.frame")
  expect_equal(nrow(winratio_example), 180L)
  expect_true(all(c("treat", "U1", "Delta1", "U2", "Delta2", "x1", "x2") %in% names(winratio_example)))
  expect_true(all(winratio_example$U2 <= winratio_example$U1 + 1e-12))
  expect_true(all(winratio_example$Delta1 %in% c(0, 1)))
  expect_true(all(winratio_example$Delta2 %in% c(0, 1)))
})

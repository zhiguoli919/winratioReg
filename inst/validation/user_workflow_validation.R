############################################################
## STAGE 5 USER-WORKFLOW VALIDATION
############################################################

cat("Loading synthetic example data ...\n")
data("winratio_example", package = "winratioReg", envir = environment())
stopifnot(is.data.frame(winratio_example))
stopifnot(nrow(winratio_example) == 180L)
stopifnot(all(winratio_example$U2 <= winratio_example$U1 + 1e-12))

cat("Fitting Local IPCW example ...\n")
fit_local <- winratio_local(
  data = winratio_example,
  arm = "treat",
  time1 = "U1", status1 = "Delta1",
  time2 = "U2", status2 = "Delta2",
  covariates = "x1",
  tau = 4,
  grid_length = 7,
  bandwidth = "rot",
  B = 0
)
stopifnot(inherits(fit_local, "winratio_local"))
stopifnot(all(is.finite(fit_local$estimates$WR)))
stopifnot(all(fit_local$estimates$WR > 0))
loc_pred <- predict(fit_local, data.frame(x1 = c(-0.5, 0, 0.5)))
stopifnot(length(loc_pred) == 3L, all(is.finite(loc_pred)), all(loc_pred > 0))
loc_blp <- blp(fit_local)
stopifnot(nrow(loc_blp) == 2L, all(is.finite(loc_blp$estimate)))

cat("Fitting Gaussian PH--copula example ...\n")
fit_ph <- winratio_phcopula(
  data = winratio_example,
  arm = "treat",
  time1 = "U1", status1 = "Delta1",
  time2 = "U2", status2 = "Delta2",
  formula = ~ x1,
  tau = 4,
  copula = "gaussian",
  degree = 2,
  grid_length = 7,
  B = 0,
  maxit = 800
)
stopifnot(inherits(fit_ph, "winratio_phcopula"))
stopifnot(all(is.finite(fit_ph$estimates$WR)))
stopifnot(all(fit_ph$estimates$WR > 0))
stopifnot(is.finite(as.numeric(logLik(fit_ph))))
ph_pred <- predict(fit_ph, data.frame(x1 = c(-0.5, 0, 0.5)))
stopifnot(length(ph_pred) == 3L, all(is.finite(ph_pred)), all(ph_pred > 0))
ph_blp <- blp(fit_ph)
stopifnot(nrow(ph_blp) == 2L, all(is.finite(ph_blp$estimate)))

cat("Comparing two working copulas ...\n")
cmp <- compare_copulas(
  data = winratio_example,
  arm = "treat",
  time1 = "U1", status1 = "Delta1",
  time2 = "U2", status2 = "Delta2",
  formula = ~ x1,
  tau = 4,
  copulas = c("frank", "gaussian"),
  degree = 2,
  grid_length = 7,
  maxit = 800
)
stopifnot(inherits(cmp, "winratio_copula_comparison"))
stopifnot(nrow(cmp) == 2L)
stopifnot(all(c("frank", "gaussian") %in% cmp$family))
stopifnot(all(is.finite(cmp$logLik)))

cat("\n============================================================\n")
cat("PASS: Example data and the primary user-facing workflows completed successfully.\n")
cat("============================================================\n")

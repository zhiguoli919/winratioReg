#' PH--copula estimation of a covariate-specific win ratio
#'
#' Fits arm-specific proportional hazards marginal models with Bernstein-sieve
#' baseline cumulative hazards and a copula model for dependence between the
#' prioritized event times. The fitted model is used to estimate the fixed-
#' horizon covariate-specific win ratio.
#'
#' @param data A data frame with one row per subject.
#' @param arm Column name identifying treatment arm.
#' @param time1,status1 Column names for the observed primary terminal-event
#'   follow-up time and its 0/1 event indicator.
#' @param time2,status2 Column names for the observed secondary-event time and
#'   its 0/1 event indicator. The observation scheme requires `time2 <= time1`.
#' @param formula One-sided formula defining the PH covariate basis, for example
#'   `~ age`, `~ x1 + x2`, or `~ splines::ns(age, df = 3)`.
#' @param tau Positive fixed analysis horizon.
#' @param treatment,control Values of `arm` identifying treatment and control.
#' @param copula Working survival-copula family. Supported values are
#'   `"bb1"`, `"clayton"`, `"frank"`, `"gaussian"`, and `"gumbel"`.
#' @param degree Bernstein polynomial degree. A single positive integer uses
#'   the same degree for both margins; a length-two vector specifies primary
#'   and secondary degrees separately. Use `"auto"` to select one common
#'   degree for both margins by BIC under the specified copula family.
#' @param degree_candidates Candidate common Bernstein degrees considered when
#'   `degree = "auto"`. Ignored for explicitly specified degrees.
#' @param sieve_upper Upper endpoint for the Bernstein sieve. By default the
#'   largest of the supplied follow-up times and `tau` is used.
#' @param grid Optional raw-covariate evaluation grid. For one underlying
#'   covariate, a numeric vector is accepted. For two, use a two-column matrix
#'   or data frame.
#' @param grid_length Number of grid points per raw covariate when `grid` is
#'   NULL. Defaults to 21 in one dimension and 5 per dimension in two.
#' @param grid_probs Quantile range used to construct an automatic grid.
#' @param B Number of arm-specific Exp(1) multiplier-bootstrap replicates.
#'   Set to zero for point estimation only.
#' @param conf.level Confidence level for pointwise intervals and bootstrap tests.
#' @param projection_weights Optional fixed weights over the raw-covariate grid.
#' @param seed Optional random seed for multiplier bootstrap.
#' @param multistart Logical; use multiple dependence-parameter starts for the
#'   original likelihood fit.
#' @param maxit Maximum optimizer iterations for the original arm fits.
#' @param bootstrap_maxit Maximum optimizer iterations for bootstrap refits.
#' @param integration_reltol Relative tolerance used in numerical integration
#'   of win/loss probabilities.
#' @param integration_subdivisions Maximum integration subdivisions.
#' @param verbose Logical; print fitting/bootstrap progress.
#'
#' @return An object of class `winratio_phcopula`.
#' @export
winratio_phcopula <- function(data, arm, time1, status1, time2, status2,
                              formula, tau, treatment = 1, control = 0,
                              copula = "bb1",
                              degree = 6L, degree_candidates = 3:7,
                              sieve_upper = NULL,
                              grid = NULL, grid_length = NULL,
                              grid_probs = c(0.05, 0.95),
                              B = 200L, conf.level = 0.95,
                              projection_weights = NULL,
                              seed = NULL, multistart = TRUE,
                              maxit = 2000L, bootstrap_maxit = 1200L,
                              integration_reltol = 1e-7,
                              integration_subdivisions = 300L,
                              verbose = FALSE) {
  call <- match.call()
  copula <- tolower(as.character(copula))
  if (length(copula) != 1L || !copula %in% c("bb1", "clayton", "frank", "gaussian", "gumbel")) {
    stop("copula must be one of 'bb1', 'clayton', 'frank', 'gaussian', or 'gumbel'.", call. = FALSE)
  }
  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau <= 0) {
    stop("tau must be a positive finite number.", call. = FALSE)
  }
  degree_auto <- .ph_is_auto_degree(degree)
  if (degree_auto) {
    degree_candidates <- .ph_validate_degree_candidates(degree_candidates)
    m1 <- m2 <- NA_integer_
  } else {
    degree_fixed <- .ph_parse_fixed_degree(degree)
    m1 <- degree_fixed[1L]
    m2 <- degree_fixed[2L]
  }
  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B < 0 || B != as.integer(B)) {
    stop("B must be a nonnegative integer.", call. = FALSE)
  }
  B <- as.integer(B)
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("conf.level must lie strictly between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(maxit) || maxit < 1 || maxit != as.integer(maxit) ||
      !is.numeric(bootstrap_maxit) || bootstrap_maxit < 1 || bootstrap_maxit != as.integer(bootstrap_maxit)) {
    stop("maxit and bootstrap_maxit must be positive integers.", call. = FALSE)
  }
  if (!is.numeric(integration_reltol) || length(integration_reltol) != 1L ||
      !is.finite(integration_reltol) || integration_reltol <= 0) {
    stop("integration_reltol must be positive.", call. = FALSE)
  }
  if (!is.numeric(integration_subdivisions) || length(integration_subdivisions) != 1L ||
      integration_subdivisions < 50 || integration_subdivisions != as.integer(integration_subdivisions)) {
    stop("integration_subdivisions must be an integer of at least 50.", call. = FALSE)
  }

  degree_selection <- NULL
  if (degree_auto) {
    sel <- .ph_degree_selection_table(
      data = data, arm = arm, time1 = time1, status1 = status1,
      time2 = time2, status2 = status2, formula = formula, tau = tau,
      treatment = treatment, control = control, copula = copula,
      sieve_upper = sieve_upper, degree_candidates = degree_candidates,
      multistart = multistart, maxit = as.integer(maxit), verbose = verbose
    )
    m1 <- m2 <- sel$selected
    prep <- sel$object$prep
    fit1 <- sel$object$fit1
    fit0 <- sel$object$fit0
    degree_selection <- list(
      method = "BIC",
      common_degree = TRUE,
      candidates = sel$candidates,
      selected = c(primary = m1, secondary = m2),
      table = sel$table
    )
  } else {
    prep <- .ph_prepare_data(
      data, arm, time1, status1, time2, status2, formula,
      treatment, control, tau, sieve_upper, m1, m2
    )
    degree_selection <- list(
      method = "fixed",
      common_degree = identical(m1, m2),
      candidates = NULL,
      selected = c(primary = m1, secondary = m2),
      table = NULL
    )
  }
  if (prep$n_omitted > 0L) {
    warning(prep$n_omitted, " incomplete observations were omitted.", call. = FALSE)
  }

  raw_grid <- .wr_prepare_grid(prep$raw_X, prep$basis$raw_vars,
                               grid, grid_length, grid_probs)
  grid_df <- as.data.frame(raw_grid, check.names = FALSE)
  names(grid_df) <- prep$basis$raw_vars
  X_grid <- .ph_basis_new(prep$basis, grid_df)
  qweight <- .local_projection_weights(raw_grid, projection_weights)

  if (!degree_auto) {
    if (isTRUE(verbose)) message("Fitting treatment-arm PH--", toupper(copula), " model...")
    fit1 <- .ph_fit_arm(prep$dat1, copula, m1, m2,
                        multistart = multistart, maxit = as.integer(maxit))
    if (isTRUE(verbose)) message("Fitting control-arm PH--", toupper(copula), " model...")
    fit0 <- .ph_fit_arm(prep$dat0, copula, m1, m2,
                        multistart = multistart, maxit = as.integer(maxit))
  } else if (isTRUE(verbose)) {
    message("Selected common Bernstein degree ", m1, " by BIC; reusing the selected arm fits.")
  }

  logWR_hat <- .ph_logwr_grid(
    fit1, fit0, X_grid, tau,
    rel.tol = integration_reltol,
    subdivisions = as.integer(integration_subdivisions)
  )
  if (any(!is.finite(logWR_hat))) {
    stop("The fitted PH--copula win ratio is invalid at one or more grid points.", call. = FALSE)
  }
  WR_hat <- exp(logWR_hat)
  proj_const <- .local_project(logWR_hat, raw_grid, qweight, "constant")
  proj_linear <- .local_project(logWR_hat, raw_grid, qweight, "linear")

  boot <- NULL
  if (B > 0L) {
    boot <- .ph_bootstrap(
      prep$dat1, prep$dat0, fit1, fit0, X_grid, raw_grid, tau,
      qweight, B, conf.level, seed, prep$n, m1, m2,
      as.integer(bootstrap_maxit), integration_reltol,
      as.integer(integration_subdivisions), verbose
    )
  }

  est <- as.data.frame(raw_grid, check.names = FALSE)
  est$logWR <- logWR_hat
  est$WR <- WR_hat
  if (is.null(boot)) {
    est$SE_logWR <- NA_real_
    est$Normal_lower <- NA_real_
    est$Normal_upper <- NA_real_
    est$Basic_lower <- NA_real_
    est$Basic_upper <- NA_real_
  } else {
    est$SE_logWR <- boot$se
    est$Normal_lower <- exp(boot$normal_lower)
    est$Normal_upper <- exp(boot$normal_upper)
    est$Basic_lower <- exp(boot$basic_lower)
    est$Basic_upper <- exp(boot$basic_upper)
  }

  total_loglik <- fit1$loglik + fit0$loglik
  total_npar <- fit1$npar + fit0$npar
  total_AIC <- 2 * total_npar - 2 * total_loglik
  total_BIC <- log(prep$n) * total_npar - 2 * total_loglik

  structure(list(
    call = call,
    method = "PH--copula",
    family = copula,
    tau = tau,
    formula = prep$basis$formula,
    covariates = prep$basis$raw_vars,
    basis_names = prep$basis$basis_names,
    dimension = length(prep$basis$raw_vars),
    treatment = treatment,
    control = control,
    degree = c(primary = m1, secondary = m2),
    degree_selection = degree_selection,
    sieve_upper = prep$sieve_upper,
    data_info = list(n = prep$n, n_treatment = prep$n1, n_control = prep$n0,
                     n_omitted = prep$n_omitted),
    grid = raw_grid,
    X_grid = X_grid,
    estimates = est,
    fit1 = fit1,
    fit0 = fit0,
    loglik = total_loglik,
    npar = total_npar,
    AIC = total_AIC,
    BIC = total_BIC,
    projection_weights = qweight,
    projection = list(constant = proj_const, linear = proj_linear),
    bootstrap = boot,
    B = B,
    conf.level = conf.level,
    basis_spec = prep$basis,
    analysis_data = list(dat1 = prep$dat1, dat0 = prep$dat0),
    integration = list(reltol = integration_reltol,
                       subdivisions = as.integer(integration_subdivisions))
  ), class = "winratio_phcopula")
}

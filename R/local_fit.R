#' Local IPCW estimation of a covariate-specific win ratio
#'
#' Estimates a fixed-horizon covariate-specific win ratio for prioritized
#' time-to-event outcomes using local kernel smoothing and inverse probability
#' of censoring weighting. One or two continuous baseline covariates are
#' currently supported.
#'
#' @param data A data frame with one row per subject.
#' @param arm Column name identifying treatment arm.
#' @param time1,status1 Column names for the observed primary terminal-event
#'   follow-up time and its 0/1 event indicator.
#' @param time2,status2 Column names for the observed secondary-event time and
#'   its 0/1 event indicator. The observation scheme requires `time2 <= time1`.
#' @param covariates Character vector containing one or two continuous covariate
#'   column names.
#' @param tau Positive fixed analysis horizon.
#' @param treatment,control Values of `arm` identifying treatment and control.
#' @param grid Optional evaluation grid. For one covariate, a numeric vector is
#'   accepted. For two covariates, supply a two-column matrix/data frame.
#' @param grid_length Number of grid points per covariate when `grid` is NULL.
#'   Defaults to 21 in one dimension and 5 per dimension in two dimensions.
#' @param grid_probs Quantile range used to construct an automatic grid.
#' @param bandwidth Either `"stability"`, `"rot"`, or a positive numeric
#'   bandwidth (one value per covariate). The stability selector compares each
#'   candidate with an inflated reference bandwidth and holds the selected
#'   bandwidth fixed throughout bootstrap resampling.
#' @param h_mult Candidate multipliers for stability bandwidth selection.
#' @param c_ref Inflation factor for the stability reference bandwidth.
#' @param bandwidth_tol Multiplicative tolerance above the minimum stability
#'   score; the first eligible (smallest) candidate is selected.
#' @param B Number of within-arm nonparametric bootstrap resamples. Set to zero
#'   for point estimation only.
#' @param conf.level Confidence level for pointwise intervals and bootstrap tests.
#' @param projection_weights Optional fixed weights over the evaluation grid.
#' @param seed Optional random seed for bootstrap resampling.
#' @param verbose Logical; print bandwidth/bootstrap progress messages.
#'
#' @return An object of class `winratio_local` containing the estimated win-ratio
#'   curve/surface, bandwidth, bootstrap inference, projection summaries, and
#'   constancy/log-linearity tests.
#' @export
winratio_local <- function(data, arm, time1, status1, time2, status2,
                           covariates, tau, treatment = 1, control = 0,
                           grid = NULL, grid_length = NULL,
                           grid_probs = c(0.05, 0.95),
                           bandwidth = c("stability", "rot"),
                           h_mult = c(0.7, 1.0, 1.3),
                           c_ref = 1.25, bandwidth_tol = 1.10,
                           B = 200L, conf.level = 0.95,
                           projection_weights = NULL,
                           seed = NULL, verbose = FALSE) {
  call <- match.call()
  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau <= 0) {
    stop("tau must be a positive finite number.", call. = FALSE)
  }
  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B < 0 || B != as.integer(B)) {
    stop("B must be a nonnegative integer.", call. = FALSE)
  }
  B <- as.integer(B)
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("conf.level must lie strictly between 0 and 1.", call. = FALSE)
  }

  prep <- .wr_prepare_local_data(data, arm, time1, status1, time2, status2,
                                 covariates, treatment, control)
  if (prep$n_omitted > 0L) {
    warning(prep$n_omitted, " incomplete observations were omitted.", call. = FALSE)
  }
  if (!any(prep$dat1$U1 > tau) || !any(prep$dat0$U1 > tau)) {
    warning("At least one arm has no primary-event follow-up beyond tau; secondary-event information at the horizon may be weak.", call. = FALSE)
  }

  G <- .wr_prepare_grid(prep$X, prep$covariate_names, grid, grid_length, grid_probs)
  bw <- .local_resolve_bandwidth(bandwidth, prep$dat1, prep$dat0, prep$X,
                                 G, tau, h_mult, c_ref, bandwidth_tol, verbose)
  h <- bw$h
  WR_hat <- .local_estimate_grid(prep$dat1, prep$dat0, G, h, tau)
  if (any(!is.finite(WR_hat)) || any(WR_hat <= 0)) {
    stop("The Local IPCW estimate is invalid at one or more grid points. Try a larger bandwidth or a narrower evaluation grid.", call. = FALSE)
  }
  logWR_hat <- log(WR_hat)
  qweight <- .local_projection_weights(G, projection_weights)
  proj_const <- .local_project(logWR_hat, G, qweight, "constant")
  proj_linear <- .local_project(logWR_hat, G, qweight, "linear")

  boot <- NULL
  if (B > 0L) {
    boot <- .local_bootstrap(prep$dat1, prep$dat0, G, h, tau, logWR_hat,
                             qweight, B, conf.level, seed, prep$n, verbose)
  }

  est <- as.data.frame(G, check.names = FALSE)
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

  structure(list(
    call = call,
    method = "Local IPCW",
    tau = tau,
    dimension = prep$dimension,
    covariates = prep$covariate_names,
    treatment = treatment,
    control = control,
    data_info = list(n = prep$n, n_treatment = prep$n1, n_control = prep$n0,
                     n_omitted = prep$n_omitted),
    grid = G,
    estimates = est,
    bandwidth = h,
    bandwidth_selection = bw,
    projection_weights = qweight,
    projection = list(constant = proj_const, linear = proj_linear),
    bootstrap = boot,
    B = B,
    conf.level = conf.level,
    analysis_data = list(dat1 = prep$dat1, dat0 = prep$dat0)
  ), class = "winratio_local")
}

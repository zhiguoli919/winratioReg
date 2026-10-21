.ph_match_formula <- function(formula) {
  if (inherits(formula, "formula")) {
    f <- formula
  } else if (is.character(formula) && length(formula) == 1L && nzchar(formula)) {
    f <- stats::as.formula(formula)
  } else {
    stop("formula must be a one-sided formula such as ~ age or ~ x1 + x2.", call. = FALSE)
  }
  if (length(f) != 2L) {
    stop("formula must be one-sided, for example ~ age or ~ x1 + x2.", call. = FALSE)
  }
  f
}

.ph_build_basis_spec <- function(data, formula) {
  formula <- .ph_match_formula(formula)
  raw_vars <- all.vars(formula)
  if (length(raw_vars) < 1L || length(raw_vars) > 2L) {
    stop("The PH--copula interface currently supports one or two underlying continuous covariates.", call. = FALSE)
  }
  missing_vars <- setdiff(raw_vars, names(data))
  if (length(missing_vars) > 0L) {
    stop("Variables in formula not found in data: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  raw <- as.matrix(data[, raw_vars, drop = FALSE])
  storage.mode(raw) <- "double"
  if (any(!is.finite(raw))) {
    stop("Covariates appearing in formula must be finite numeric values.", call. = FALSE)
  }
  sds <- apply(raw, 2L, stats::sd)
  if (any(!is.finite(sds)) || any(sds <= 0)) {
    stop("Each underlying covariate must have positive sample standard deviation.", call. = FALSE)
  }

  mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  trm <- stats::terms(mf)
  X <- stats::model.matrix(trm, data = mf)
  if ("(Intercept)" %in% colnames(X)) {
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  }
  storage.mode(X) <- "double"
  if (ncol(X) < 1L || any(!is.finite(X))) {
    stop("The covariate formula produced an empty or non-finite design matrix.", call. = FALSE)
  }

  list(
    formula = formula,
    terms = trm,
    raw_vars = raw_vars,
    raw_X = raw,
    X = X,
    basis_names = colnames(X),
    contrasts = attr(X, "contrasts")
  )
}

.ph_basis_new <- function(spec, newdata) {
  if (!is.data.frame(newdata)) newdata <- as.data.frame(newdata)
  missing_vars <- setdiff(spec$raw_vars, names(newdata))
  if (length(missing_vars) > 0L) {
    stop("newdata is missing: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  mf <- stats::model.frame(spec$terms, data = newdata, na.action = stats::na.pass)
  X <- stats::model.matrix(spec$terms, data = mf, contrasts.arg = spec$contrasts)
  if ("(Intercept)" %in% colnames(X)) {
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  }
  storage.mode(X) <- "double"
  if (ncol(X) != length(spec$basis_names)) {
    stop("newdata produced a design matrix incompatible with the fitted model.", call. = FALSE)
  }
  colnames(X) <- spec$basis_names
  if (any(!is.finite(X))) stop("newdata produced non-finite basis values.", call. = FALSE)
  X
}

.ph_prepare_data <- function(data, arm, time1, status1, time2, status2,
                             formula, treatment, control, tau, sieve_upper,
                             m1, m2) {
  if (!is.data.frame(data)) stop("data must be a data.frame.", call. = FALSE)
  arm <- .wr_col_name(arm, "arm")
  time1 <- .wr_col_name(time1, "time1")
  status1 <- .wr_col_name(status1, "status1")
  time2 <- .wr_col_name(time2, "time2")
  status2 <- .wr_col_name(status2, "status2")
  formula <- .ph_match_formula(formula)
  raw_vars <- all.vars(formula)

  vars <- unique(c(arm, time1, status1, time2, status2, raw_vars))
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars) > 0L) {
    stop("Missing variables: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }

  d <- data[, vars, drop = FALSE]
  keep <- stats::complete.cases(d)
  n_omitted <- sum(!keep)
  d <- d[keep, , drop = FALSE]
  if (nrow(d) == 0L) stop("No complete observations are available.", call. = FALSE)

  if (identical(treatment, control)) {
    stop("treatment and control must identify different arm values.", call. = FALSE)
  }
  arm_values <- d[[arm]]
  allowed <- arm_values == treatment | arm_values == control
  if (any(!allowed)) {
    stop("The arm variable contains values other than treatment and control.", call. = FALSE)
  }
  if (!any(arm_values == treatment) || !any(arm_values == control)) {
    stop("Both treatment and control arms must contain observations.", call. = FALSE)
  }

  U1 <- suppressWarnings(as.numeric(d[[time1]]))
  U2 <- suppressWarnings(as.numeric(d[[time2]]))
  D1 <- suppressWarnings(as.numeric(d[[status1]]))
  D2 <- suppressWarnings(as.numeric(d[[status2]]))
  if (any(!is.finite(U1)) || any(!is.finite(U2)) || any(U1 < 0) || any(U2 < 0)) {
    stop("Event/follow-up times must be finite and nonnegative.", call. = FALSE)
  }
  if (any(!is.finite(D1)) || any(!is.finite(D2)) ||
      any(!D1 %in% c(0, 1)) || any(!D2 %in% c(0, 1))) {
    stop("status1 and status2 must be coded 0/1.", call. = FALSE)
  }
  if (any(U2 > U1 + 1e-10)) {
    stop("The observation scheme requires time2 <= time1 for every subject.", call. = FALSE)
  }

  basis <- .ph_build_basis_spec(d, formula)
  if (nrow(basis$X) != nrow(d)) stop("Could not construct the covariate design matrix.", call. = FALSE)

  if (is.null(sieve_upper)) {
    sieve_upper <- max(c(U1, U2, tau), na.rm = TRUE)
  }
  if (!is.numeric(sieve_upper) || length(sieve_upper) != 1L ||
      !is.finite(sieve_upper) || sieve_upper <= 0) {
    stop("sieve_upper must be a positive finite number.", call. = FALSE)
  }
  if (sieve_upper + 1e-10 < max(c(U1, U2), na.rm = TRUE)) {
    stop("sieve_upper must be at least as large as all supplied follow-up times.", call. = FALSE)
  }
  if (sieve_upper + 1e-10 < tau) {
    stop("sieve_upper must be at least tau.", call. = FALSE)
  }

  internal <- data.frame(
    U1 = U1,
    Delta1 = as.integer(D1),
    U2 = U2,
    Delta2 = as.integer(D2),
    check.names = FALSE
  )
  idx1 <- arm_values == treatment
  idx0 <- arm_values == control

  prep_arm <- function(idx) {
    out <- internal[idx, , drop = FALSE]
    X <- basis$X[idx, , drop = FALSE]
    attr(out, "X_mat") <- X
    attr(out, "tmax") <- sieve_upper
    attr(out, "B1_m") <- .ph_bernstein_basis(out$U1, m = m1, lower = 0, upper = sieve_upper)
    attr(out, "B1_m1") <- .ph_bernstein_basis(out$U1, m = m1 - 1L, lower = 0, upper = sieve_upper)
    attr(out, "B2_m") <- .ph_bernstein_basis(out$U2, m = m2, lower = 0, upper = sieve_upper)
    attr(out, "B2_m1") <- .ph_bernstein_basis(out$U2, m = m2 - 1L, lower = 0, upper = sieve_upper)
    out
  }

  list(
    dat1 = prep_arm(idx1),
    dat0 = prep_arm(idx0),
    raw_X = basis$raw_X,
    basis = basis,
    arm = arm_values,
    n = nrow(d),
    n1 = sum(idx1),
    n0 = sum(idx0),
    n_omitted = n_omitted,
    sieve_upper = sieve_upper
  )
}

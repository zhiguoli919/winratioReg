.wr_col_name <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || !nzchar(x)) {
    stop(arg, " must be a single column name.", call. = FALSE)
  }
  x
}

.wr_prepare_local_data <- function(data, arm, time1, status1, time2, status2,
                                   covariates, treatment, control) {
  if (!is.data.frame(data)) {
    stop("data must be a data.frame.", call. = FALSE)
  }

  arm <- .wr_col_name(arm, "arm")
  time1 <- .wr_col_name(time1, "time1")
  status1 <- .wr_col_name(status1, "status1")
  time2 <- .wr_col_name(time2, "time2")
  status2 <- .wr_col_name(status2, "status2")

  if (!is.character(covariates) || length(covariates) < 1L ||
      length(covariates) > 2L || any(!nzchar(covariates))) {
    stop("covariates must contain one or two column names.", call. = FALSE)
  }

  vars <- unique(c(arm, time1, status1, time2, status2, covariates))
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
  X <- as.matrix(d[, covariates, drop = FALSE])
  storage.mode(X) <- "double"

  if (any(!is.finite(U1)) || any(!is.finite(U2)) || any(U1 < 0) || any(U2 < 0)) {
    stop("Event/follow-up times must be finite and nonnegative.", call. = FALSE)
  }
  if (any(!is.finite(D1)) || any(!is.finite(D2)) ||
      any(!D1 %in% c(0, 1)) || any(!D2 %in% c(0, 1))) {
    stop("status1 and status2 must be coded 0/1.", call. = FALSE)
  }
  if (any(!is.finite(X))) {
    stop("Covariates must be finite numeric values.", call. = FALSE)
  }
  if (any(U2 > U1 + 1e-10)) {
    stop("The observation scheme requires time2 <= time1 for every subject.", call. = FALSE)
  }

  sds <- apply(X, 2L, stats::sd)
  if (any(!is.finite(sds)) || any(sds <= 0)) {
    stop("Each covariate must have positive sample standard deviation.", call. = FALSE)
  }

  internal <- data.frame(
    X,
    U1 = U1,
    Delta1 = as.integer(D1),
    U2 = U2,
    Delta2 = as.integer(D2),
    check.names = FALSE
  )
  names(internal)[seq_along(covariates)] <- paste0("X", seq_along(covariates))

  idx1 <- arm_values == treatment
  idx0 <- arm_values == control

  list(
    dat1 = internal[idx1, , drop = FALSE],
    dat0 = internal[idx0, , drop = FALSE],
    X = X,
    arm = arm_values,
    covariate_names = covariates,
    dimension = length(covariates),
    n = nrow(internal),
    n1 = sum(idx1),
    n0 = sum(idx0),
    n_omitted = n_omitted
  )
}

.wr_prepare_grid <- function(X, covariate_names, grid = NULL,
                             grid_length = NULL, grid_probs = c(0.05, 0.95)) {
  p <- ncol(X)
  if (length(grid_probs) != 2L || any(!is.finite(grid_probs)) ||
      grid_probs[1L] < 0 || grid_probs[2L] > 1 || grid_probs[1L] >= grid_probs[2L]) {
    stop("grid_probs must be two increasing probabilities in [0, 1].", call. = FALSE)
  }

  if (is.null(grid)) {
    if (is.null(grid_length)) grid_length <- if (p == 1L) 21L else 5L
    if (length(grid_length) != 1L || !is.finite(grid_length) || grid_length < 2L) {
      stop("grid_length must be at least 2.", call. = FALSE)
    }
    grid_length <- as.integer(grid_length)

    seqs <- lapply(seq_len(p), function(j) {
      qq <- stats::quantile(X[, j], probs = grid_probs, na.rm = TRUE, names = FALSE)
      seq(qq[1L], qq[2L], length.out = grid_length)
    })
    if (p == 1L) {
      G <- matrix(seqs[[1L]], ncol = 1L)
    } else {
      G <- as.matrix(expand.grid(seqs[[1L]], seqs[[2L]]))
    }
  } else {
    if (p == 1L && is.atomic(grid) && is.null(dim(grid))) {
      G <- matrix(as.numeric(grid), ncol = 1L)
    } else {
      G <- as.matrix(grid)
      storage.mode(G) <- "double"
      if (ncol(G) != p) {
        stop("The supplied grid must have one column per covariate.", call. = FALSE)
      }
    }
  }

  storage.mode(G) <- "double"
  if (nrow(G) < p + 1L || any(!is.finite(G))) {
    stop("The evaluation grid is too small or contains non-finite values.", call. = FALSE)
  }
  if (p == 1L) G <- G[order(G[, 1L]), , drop = FALSE]
  colnames(G) <- covariate_names
  G
}

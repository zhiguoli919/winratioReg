.local_rot_bandwidth <- function(X, test_rule = TRUE) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  mult <- if (isTRUE(test_rule) && p == 1L) 1.06 else 1
  h <- apply(X, 2L, stats::sd) * n^(-1 / 5) * mult
  as.numeric(h)
}

.local_select_bandwidth <- function(dat1, dat0, X_all, grid, tau,
                                    h_mult, c_ref, tol, verbose) {
  if (length(h_mult) < 1L || any(!is.finite(h_mult)) || any(h_mult <= 0)) {
    stop("h_mult must contain positive finite multipliers.", call. = FALSE)
  }
  h_mult <- sort(unique(as.numeric(h_mult)))
  if (!is.finite(c_ref) || c_ref <= 1) stop("c_ref must exceed 1.", call. = FALSE)
  if (!is.finite(tol) || tol < 1) stop("bandwidth_tol must be at least 1.", call. = FALSE)

  h0 <- .local_rot_bandwidth(X_all, test_rule = FALSE)
  scores <- rep(Inf, length(h_mult))
  candidates <- lapply(h_mult, function(a) a * h0)

  for (k in seq_along(candidates)) {
    h <- candidates[[k]]
    wr_h <- .local_estimate_grid(dat1, dat0, grid, h, tau)
    wr_ref <- .local_estimate_grid(dat1, dat0, grid, c_ref * h, tau)
    good <- is.finite(wr_h) & is.finite(wr_ref) & wr_h > 0 & wr_ref > 0
    if (any(good)) scores[k] <- mean((log(wr_h[good]) - log(wr_ref[good]))^2)
    if (isTRUE(verbose)) {
      message("Bandwidth candidate ", k, "/", length(candidates),
              ": multiplier=", signif(h_mult[k], 4),
              ", score=", signif(scores[k], 6))
    }
  }

  if (!any(is.finite(scores))) {
    stop("All stability-bandwidth candidates failed. Try a wider/narrower grid or specify bandwidth manually.", call. = FALSE)
  }
  qmin <- min(scores, na.rm = TRUE)
  eligible <- which(is.finite(scores) & scores <= tol * qmin)
  idx <- eligible[1L]

  list(
    method = "stability",
    h0 = h0,
    h = candidates[[idx]],
    multipliers = h_mult,
    scores = scores,
    selected_index = idx,
    c_ref = c_ref,
    tol = tol
  )
}

.local_resolve_bandwidth <- function(bandwidth, dat1, dat0, X_all, grid, tau,
                                     h_mult, c_ref, tol, verbose) {
  p <- ncol(X_all)
  if (is.numeric(bandwidth)) {
    h <- as.numeric(bandwidth)
    if (length(h) == 1L && p > 1L) h <- rep(h, p)
    if (length(h) != p || any(!is.finite(h)) || any(h <= 0)) {
      stop("A numeric bandwidth must contain one positive value per covariate.", call. = FALSE)
    }
    return(list(method = "manual", h = h))
  }

  method <- match.arg(bandwidth, c("stability", "rot"))
  if (method == "rot") {
    return(list(method = "rot", h = .local_rot_bandwidth(X_all, test_rule = TRUE)))
  }
  .local_select_bandwidth(dat1, dat0, X_all, grid, tau, h_mult, c_ref, tol, verbose)
}

.local_X <- function(dat) {
  as.matrix(dat[, grep("^X[0-9]+$", names(dat)), drop = FALSE])
}

.local_kernel_weights <- function(X, x0, h, eps = 1e-12) {
  X <- as.matrix(X)
  x0 <- as.numeric(x0)
  h <- as.numeric(h)
  if (ncol(X) != length(x0) || length(h) != ncol(X) || any(h <= 0)) {
    stop("Incompatible covariate point or bandwidth.", call. = FALSE)
  }
  U <- sweep(X, 2L, x0, FUN = "-")
  U <- sweep(U, 2L, h, FUN = "/")
  K <- matrix(stats::dnorm(as.vector(U)), nrow = nrow(U), ncol = ncol(U))
  w <- apply(K, 1L, prod)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= eps) return(rep(1 / nrow(X), nrow(X)))
  w / sw
}

.local_weighted_G <- function(times, dat, w, eps = 1e-8) {
  times <- as.numeric(times)
  if (length(times) == 0L) return(numeric(0))

  censor <- as.integer(dat$Delta1 == 0L)
  censor_times <- sort(unique(dat$U1[censor == 1L]))
  if (length(censor_times) == 0L) return(rep(1, length(times)))

  Gvals <- numeric(length(censor_times))
  G <- 1
  for (j in seq_along(censor_times)) {
    tj <- censor_times[j]
    risk <- sum(w[dat$U1 >= tj])
    dcen <- sum(w[dat$U1 == tj & censor == 1L])
    if (is.finite(risk) && risk > eps) {
      G <- G * (1 - dcen / risk)
    }
    Gvals[j] <- max(G, eps)
  }

  out <- vapply(times, function(tt) {
    idx <- which(censor_times <= tt)
    if (length(idx) == 0L) 1 else Gvals[max(idx)]
  }, numeric(1))
  pmax(out, eps)
}

.local_S1 <- function(times, dat, w) {
  if (length(times) == 0L) return(numeric(0))
  G <- .local_weighted_G(times, dat, w)
  vapply(times, function(tt) sum(w * (dat$U1 > tt)), numeric(1)) / G
}

.local_dF1 <- function(event_times, dat, w) {
  if (length(event_times) == 0L) return(numeric(0))
  G <- .local_weighted_G(event_times, dat, w)
  mass <- vapply(event_times, function(tt) {
    sum(w * (dat$Delta1 == 1L) * (dat$U1 == tt))
  }, numeric(1))
  mass / G
}

.local_J <- function(tau, u, dat, w) {
  if (length(u) == 0L) return(numeric(0))
  Gtau <- .local_weighted_G(tau, dat, w)[1L]
  vapply(u, function(uu) {
    sum(w * (dat$U1 > tau) * (dat$U2 > uu))
  }, numeric(1)) / Gtau
}

.local_minus_dJ <- function(tau, event_times, dat, w) {
  if (length(event_times) == 0L) return(numeric(0))
  Gtau <- .local_weighted_G(tau, dat, w)[1L]
  mass <- vapply(event_times, function(uu) {
    sum(w * (dat$U1 > tau) * (dat$Delta2 == 1L) * (dat$U2 == uu))
  }, numeric(1))
  mass / Gtau
}

.local_wr_one <- function(x0, dat1, dat0, h, tau) {
  X1 <- .local_X(dat1)
  X0 <- .local_X(dat0)
  w1 <- .local_kernel_weights(X1, x0, h)
  w0 <- .local_kernel_weights(X0, x0, h)

  e10 <- sort(unique(dat0$U1[dat0$Delta1 == 1L & dat0$U1 <= tau]))
  e11 <- sort(unique(dat1$U1[dat1$Delta1 == 1L & dat1$U1 <= tau]))

  W_primary <- if (length(e10) > 0L) {
    sum(.local_S1(e10, dat1, w1) * .local_dF1(e10, dat0, w0))
  } else 0
  L_primary <- if (length(e11) > 0L) {
    sum(.local_S1(e11, dat0, w0) * .local_dF1(e11, dat1, w1))
  } else 0

  e20 <- sort(unique(dat0$U2[
    dat0$Delta2 == 1L & dat0$U2 <= tau & dat0$U1 > tau
  ]))
  e21 <- sort(unique(dat1$U2[
    dat1$Delta2 == 1L & dat1$U2 <= tau & dat1$U1 > tau
  ]))

  W_secondary <- if (length(e20) > 0L) {
    sum(.local_J(tau, e20, dat1, w1) * .local_minus_dJ(tau, e20, dat0, w0))
  } else 0
  L_secondary <- if (length(e21) > 0L) {
    sum(.local_J(tau, e21, dat0, w0) * .local_minus_dJ(tau, e21, dat1, w1))
  } else 0

  W <- W_primary + W_secondary
  L <- L_primary + L_secondary
  if (!is.finite(W) || !is.finite(L) || W <= 0 || L <= 0) return(NA_real_)
  W / L
}

.local_estimate_grid <- function(dat1, dat0, grid, h, tau) {
  apply(grid, 1L, function(x0) {
    .local_wr_one(as.numeric(x0), dat1, dat0, h, tau)
  })
}

make_local_test_data <- function(n_arm = 60L) {
  x0 <- seq(-1, 1, length.out = n_arm)
  x1 <- seq(-1, 1, length.out = n_arm)
  k <- seq_len(n_arm)
  T10 <- 1.45 + 0.45 * (x0 + 1) + 0.12 * sin(k)
  T11 <- 1.52 + 0.42 * (x1 + 1) + 0.12 * cos(k)
  T20 <- pmax(0.15, 0.65 + 0.20 * (x0 + 1) + 0.08 * cos(k))
  T21 <- pmax(0.15, 0.70 + 0.18 * (x1 + 1) + 0.08 * sin(k))
  C0 <- rep(4, n_arm)
  C1 <- rep(4, n_arm)
  data.frame(
    arm = rep(c(0, 1), each = n_arm),
    x = c(x0, x1),
    U1 = c(pmin(T10, C0), pmin(T11, C1)),
    Delta1 = 1L,
    U2 = c(pmin(T20, T10, C0), pmin(T21, T11, C1)),
    Delta2 = 1L
  )
}

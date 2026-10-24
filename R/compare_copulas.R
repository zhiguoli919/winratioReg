#' Compare candidate PH--copula families
#'
#' Fits the same PH marginal specification and Bernstein degrees under several
#' copula families and returns likelihood-based fit summaries. This function
#' performs point estimation only; use [winratio_phcopula()] for bootstrap
#' inference after choosing a working family.
#'
#' @inheritParams winratio_phcopula
#' @param copulas Character vector of candidate families selected from
#'   `"bb1"`, `"clayton"`, `"frank"`, `"gaussian"`, and `"gumbel"`.
#' @details `degree` must be explicitly specified for copula comparison.
#'   Automatic degree selection is intentionally not performed separately for
#'   each copula family, because doing so would confound copula-family
#'   comparison with changes in marginal sieve flexibility. If desired, first
#'   fit one working family with `winratio_phcopula(..., degree = "auto")`,
#'   then pass its selected `degree` to `compare_copulas()`.
#' @return A data frame with total log likelihood, AIC, BIC, convergence, and
#'   arm-specific dependence summaries. The fitted objects are stored in the
#'   `fits` attribute.
#' @export
compare_copulas <- function(data, arm, time1, status1, time2, status2,
                            formula, tau, treatment = 1, control = 0,
                            copulas = c("bb1", "clayton", "frank", "gaussian", "gumbel"),
                            degree = 6L, sieve_upper = NULL,
                            grid = NULL, grid_length = NULL,
                            grid_probs = c(0.05, 0.95),
                            projection_weights = NULL,
                            multistart = TRUE, maxit = 2000L,
                            integration_method = c("simpson", "adaptive"),
                            integration_grid = 401L,
                            integration_reltol = 1e-7,
                            integration_subdivisions = 300L,
                            verbose = FALSE) {
  integration_method <- match.arg(integration_method)
  if (.ph_is_auto_degree(degree)) {
    stop(
      "compare_copulas() requires an explicit Bernstein degree so every copula is compared with the same marginal sieve specification. Select a degree first with winratio_phcopula(..., degree = 'auto', B = 0), then pass fit$degree here.",
      call. = FALSE
    )
  }
  degree <- .ph_parse_fixed_degree(degree)

  copulas <- unique(tolower(as.character(copulas)))
  bad <- setdiff(copulas, c("bb1", "clayton", "frank", "gaussian", "gumbel"))
  if (length(copulas) < 1L || length(bad) > 0L) {
    stop("copulas must be selected from 'bb1', 'clayton', 'frank', 'gaussian', and 'gumbel'.", call. = FALSE)
  }

  fits <- vector("list", length(copulas))
  names(fits) <- copulas
  rows <- vector("list", length(copulas))
  for (k in seq_along(copulas)) {
    fam <- copulas[k]
    if (isTRUE(verbose)) message("Comparing copula family: ", fam)
    z <- try(winratio_phcopula(
      data = data, arm = arm, time1 = time1, status1 = status1,
      time2 = time2, status2 = status2, formula = formula, tau = tau,
      treatment = treatment, control = control, copula = fam,
      degree = degree, sieve_upper = sieve_upper,
      grid = grid, grid_length = grid_length, grid_probs = grid_probs,
      B = 0L, projection_weights = projection_weights,
      multistart = multistart, maxit = maxit,
      integration_method = integration_method,
      integration_grid = integration_grid,
      integration_reltol = integration_reltol,
      integration_subdivisions = integration_subdivisions,
      verbose = verbose
    ), silent = TRUE)
    fits[[k]] <- z
    if (inherits(z, "try-error")) {
      rows[[k]] <- data.frame(
        family = fam, logLik = NA_real_, AIC = Inf, BIC = Inf,
        control_theta = NA_real_, control_rho = NA_real_, control_delta = NA_real_, control_tau = NA_real_,
        treatment_theta = NA_real_, treatment_rho = NA_real_, treatment_delta = NA_real_, treatment_tau = NA_real_,
        control_convergence = NA_integer_, treatment_convergence = NA_integer_,
        stringsAsFactors = FALSE
      )
    } else {
      rows[[k]] <- data.frame(
        family = fam, logLik = z$loglik, AIC = z$AIC, BIC = z$BIC,
        control_theta = if (fam == "gaussian") NA_real_ else z$fit0$theta,
        control_rho = if (fam == "gaussian") z$fit0$theta else NA_real_,
        control_delta = if (is.null(z$fit0$delta)) NA_real_ else z$fit0$delta,
        control_tau = z$fit0$kendall_tau,
        treatment_theta = if (fam == "gaussian") NA_real_ else z$fit1$theta,
        treatment_rho = if (fam == "gaussian") z$fit1$theta else NA_real_,
        treatment_delta = if (is.null(z$fit1$delta)) NA_real_ else z$fit1$delta,
        treatment_tau = z$fit1$kendall_tau,
        control_convergence = z$fit0$convergence,
        treatment_convergence = z$fit1$convergence,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "fits") <- fits
  class(out) <- c("winratio_copula_comparison", class(out))
  out
}

#' @export
print.winratio_copula_comparison <- function(x, ...) {
  y <- x
  class(y) <- "data.frame"
  print(y, row.names = FALSE, ...)
  invisible(x)
}

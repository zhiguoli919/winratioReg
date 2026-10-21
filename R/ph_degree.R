.ph_is_auto_degree <- function(degree) {
  is.character(degree) && length(degree) == 1L &&
    identical(tolower(trimws(degree)), "auto")
}

.ph_parse_fixed_degree <- function(degree) {
  z <- suppressWarnings(as.numeric(degree))
  if (length(z) == 1L) z <- rep(z, 2L)
  if (length(z) != 2L || any(!is.finite(z)) ||
      any(z < 1) || any(z != as.integer(z))) {
    stop("degree must be 'auto', one positive integer, or two positive integers.",
         call. = FALSE)
  }
  as.integer(z)
}

.ph_validate_degree_candidates <- function(degree_candidates) {
  z <- suppressWarnings(as.numeric(degree_candidates))
  if (length(z) < 1L || any(!is.finite(z)) ||
      any(z < 1) || any(z != as.integer(z))) {
    stop("degree_candidates must contain one or more positive integers.",
         call. = FALSE)
  }
  sort(unique(as.integer(z)))
}

.ph_degree_selection_table <- function(data, arm, time1, status1, time2, status2,
                                       formula, tau, treatment, control, copula,
                                       sieve_upper, degree_candidates,
                                       multistart, maxit, verbose) {
  candidates <- .ph_validate_degree_candidates(degree_candidates)
  rows <- vector("list", length(candidates))
  objects <- vector("list", length(candidates))

  for (k in seq_along(candidates)) {
    m <- candidates[k]
    if (isTRUE(verbose)) {
      message("Evaluating Bernstein degree ", m, " (common to both margins) by BIC...")
    }

    z <- try({
      prep <- .ph_prepare_data(
        data, arm, time1, status1, time2, status2, formula,
        treatment, control, tau, sieve_upper, m, m
      )
      fit1 <- .ph_fit_arm(
        prep$dat1, copula, m, m,
        multistart = multistart, maxit = as.integer(maxit)
      )
      fit0 <- .ph_fit_arm(
        prep$dat0, copula, m, m,
        multistart = multistart, maxit = as.integer(maxit)
      )
      total_loglik <- fit1$loglik + fit0$loglik
      total_npar <- fit1$npar + fit0$npar
      total_AIC <- 2 * total_npar - 2 * total_loglik
      total_BIC <- log(prep$n) * total_npar - 2 * total_loglik
      list(
        prep = prep,
        fit1 = fit1,
        fit0 = fit0,
        logLik = total_loglik,
        npar = total_npar,
        AIC = total_AIC,
        BIC = total_BIC,
        converged = fit1$convergence == 0L && fit0$convergence == 0L
      )
    }, silent = TRUE)

    if (inherits(z, "try-error")) {
      rows[[k]] <- data.frame(
        degree = m,
        logLik = NA_real_,
        npar = NA_integer_,
        AIC = Inf,
        BIC = Inf,
        control_convergence = NA_integer_,
        treatment_convergence = NA_integer_,
        converged = FALSE,
        stringsAsFactors = FALSE
      )
      objects[k] <- list(NULL)
    } else {
      rows[[k]] <- data.frame(
        degree = m,
        logLik = z$logLik,
        npar = z$npar,
        AIC = z$AIC,
        BIC = z$BIC,
        control_convergence = z$fit0$convergence,
        treatment_convergence = z$fit1$convergence,
        converged = z$converged,
        stringsAsFactors = FALSE
      )
      objects[[k]] <- z
    }
  }

  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL
  good <- which(tab$converged & is.finite(tab$BIC))
  if (length(good) == 0L) {
    stop(
      "Automatic Bernstein-degree selection failed because no candidate degree converged in both treatment arms.",
      call. = FALSE
    )
  }
  best <- good[which.min(tab$BIC[good])]

  list(
    selected = as.integer(tab$degree[best]),
    table = tab,
    object = objects[[best]],
    criterion = "BIC",
    candidates = candidates
  )
}

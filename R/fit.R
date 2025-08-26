#' Fit parametric power model with robust start/jitter loop
#'
#' Tries optimization from reasonable starts, jitters on failure,
#' and only accepts solutions that are converged, finite, and non-degenerate.
#'
#' @param dat Data frame with columns y, m1, m2.
#' @param alpha Significance level.
#' @param test 'onesided' or 'twosided'.
#' @param parameterization Passed to power_fun.
#' @param start Optional numeric vector of length 3 (phi). If NULL,
#'   defaults depend on \code{transform}.
#' @param transform If TRUE, \code{phi} are log-parameters (mapped via exp()).
#' @param max_tries Maximum optimization attempts.
#' @param jitter_scale Jitter amplitude for new starts.
#' @param eps_flat_sd Threshold for sd(P) to flag flat power curves.
#' @param eps_edges Threshold near 0/1 to flag degenerate power.
#' @return An object of class \code{power_model}.
#' @export
fit_model <- function(dat, alpha = 0.05, test = "onesided",
                      parameterization = "natural",
                      start = NULL,
                      transform = TRUE,
                      max_tries = 25,
                      jitter_scale = 0.5,
                      eps_flat_sd = 1e-6,
                      eps_edges = 1e-10) {

  stopifnot(is.data.frame(dat), all(c("y","m1","m2") %in% names(dat)))
  y  <- dat$y
  m1 <- dat$m1
  m2 <- dat$m2

  # choose starting point
  if (is.null(start)) {
    cur_start <- .default_start(transform)
  } else {
    cur_start <- as.numeric(start)
    if (length(cur_start) != 3L) stop("'start' must be length 3")
  }

  tries <- 0L
  success <- FALSE
  opt_phi <- NULL; V_mod <- NULL; V_rob <- NULL; conv <- NA_integer_; ll <- NA_real_
  start_used <- cur_start

  repeat {
    tries <- tries + 1L

    opt_try <- try(
      nlminb(
        start = cur_start,
        objective = negLL,
        y = y, m1 = m1, m2 = m2,
        alpha = alpha, test = test,
        parameterization = parameterization,
        transform = transform,
        lower = if (isTRUE(transform)) rep(-Inf, 3) else c(-Inf, 1e-8, 1e-8)
      ),
      silent = TRUE
    )

    bad <- FALSE
    if (inherits(opt_try, "try-error")) {
      bad <- TRUE
    } else if (!is.finite(opt_try$objective)) {
      bad <- TRUE
    } else if (opt_try$convergence != 0) {
      bad <- TRUE
    }

    if (!bad) {
      opt_phi_tmp <- opt_try$par
      pars <- stabilityTransformation(opt_phi_tmp, transform = transform)
      P_hat <- power_fun(m1 = m1, m2 = m2,
                         beta0 = pars$beta0, beta1 = pars$beta1, beta2 = pars$beta2,
                         alpha = alpha, test = test,
                         parameterization = parameterization)
      if (.is_flat_power(P_hat, eps_flat_sd = eps_flat_sd, eps_edges = eps_edges)) {
        bad <- TRUE
      }
    }

    if (!bad) {
      H <- hess_fun(opt_phi = opt_try$par, y = y, m1 = m1, m2 = m2,
                    alpha = alpha, test = test,
                    parameterization = parameterization, transform = transform)
      V_mod_try <- try(solve(H), silent = TRUE)
      if (inherits(V_mod_try, "try-error") || any(!is.finite(V_mod_try))) {
        bad <- TRUE
      } else {
        S <- score_mat(opt_phi = opt_try$par, y = y, m1 = m1, m2 = m2,
                       alpha = alpha, test = test,
                       parameterization = parameterization, transform = transform)
        B <- t(S) %*% S
        V_rob_try <- V_mod_try %*% B %*% V_mod_try

        opt_phi <- opt_try$par
        V_mod   <- V_mod_try
        V_rob   <- V_rob_try
        conv    <- opt_try$convergence
        ll      <- -opt_try$objective
        start_used <- cur_start
        success <- TRUE
      }
    }

    if (success) break
    if (tries >= max_tries) {
      stop(sprintf("fit_model: failed after %d attempts (non-converged or degenerate).", max_tries))
    }

    base <- if (!inherits(opt_try, "try-error")) as.numeric(opt_try$par) else cur_start
    cur_start <- .jitter(base, scale = jitter_scale)
  }

  out <- list(
    opt_phi = opt_phi,
    V_mod   = V_mod,
    V_rob   = V_rob,
    conv    = conv,
    ll      = ll,
    alpha   = alpha,
    test    = test,
    parameterization = parameterization,
    transform = transform,
    tries   = tries,
    start_used = start_used,
    # store data
    y  = y,
    m1 = m1,
    m2 = m2
  )
  class(out) <- "power_model"
  out
}

#' @importFrom stats AIC logLik nobs
NULL

#' AIC method for power_model objects
#'
#' @param object A 'power_model' object.
#' @param ... Further arguments passed to methods.
#' @return A numeric value with class "AIC".
#' @export
AIC.power_model <- function(object, ...) {
  ll  <- logLik(object)
  k   <- attr(ll, "df")
  aic <- -2 * as.numeric(ll) + 2 * k
  attr(aic, "df")   <- k
  attr(aic, "nobs") <- nobs(object)
  class(aic) <- "AIC"
  aic
}

#' logLik method for power_model objects
#'
#' @param object A 'power_model' object.
#' @param ... Further arguments (ignored).
#' @return An object of class "logLik".
#' @export
logLik.power_model <- function(object, ...) {
  val <- object$ll
  attr(val, "df")   <- length(object$opt_phi)
  attr(val, "nobs") <- nobs(object)
  class(val) <- "logLik"
  val
}

#' nobs method for power_model objects
#'
#' @param object A 'power_model' object.
#' @param ... Unused.
#' @return Number of observations used in the fit.
#' @export
nobs.power_model <- function(object, ...) {
  length(object$y)
}

#' Extract model data from a power_model
#'
#' @param object A 'power_model' object.
#' @return A data.frame with columns y, m1, m2.
#' @export
model_data.power_model <- function(object) {
  data.frame(y = object$y, m1 = object$m1, m2 = object$m2)
}


# --- Helpers ---------------------------------------------------------------

# check if predicted power curve is flat/degenerate
.is_flat_power <- function(P, eps_flat_sd = 1e-6, eps_edges = 1e-10) {
  if (!all(is.finite(P))) return(TRUE)
  all(P > 1 - eps_edges) || all(P < eps_edges) || (stats::sd(P) < eps_flat_sd)
}

# default start depending on transform
.default_start <- function(transform) {
  if (isTRUE(transform)) {
    # log-phi = (0,0,0) → beta0 = -1, beta1 = 1, beta2 = 1
    c(0, 0, 0)
  } else {
    # natural phi = (1,1,1) → beta0 = -1, beta1 = 1, beta2 = 1
    c(1, 1, 1)
  }
}

# jitter helper
.jitter <- function(x, scale = 0.5) {
  x + stats::runif(length(x), min = -scale, max = scale)
}


#' coef method for power_model objects
#'
#' @param object A 'power_model' object.
#' @param transform Logical, whether to back-transform from log-phi.
#'   Defaults to object$transform.
#' @param ... Unused.
#' @return Named numeric vector of coefficients (beta0, beta1, beta2).
#' @export
coef.power_model <- function(object, transform = object$transform, ...) {
  pars <- stabilityTransformation(object$opt_phi, transform = transform)
  c(beta0 = pars$beta0,
    beta1 = pars$beta1,
    beta2 = pars$beta2)
}


#' Print method for power_model objects
#'
#' @param x A 'power_model' object.
#' @param ... Unused.
#' @export
print.power_model <- function(x, ...) {
  cat("Power model fit\n")
  cat(paste0("  logLik:", round(x$ll, 3), "  (nobs =", nobs(x), ")\n"))
  cat("  Converged:", x$conv == 0, " after", x$tries, "tries\n\n")
  cat(paste0("Coefficients (",x$parameterization," parameterization):\n"))
  print(coef(x), digits = 4)
  invisible(x)  # return full object invisibly
}

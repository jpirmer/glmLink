#' Fit parametric power model
#'
#' @param dat Data frame with y, m1, m2.
#' @param alpha Significance level.
#' @param test 'onesided' or 'twosided'.
#' @param start Named numeric vector of log-parameters.
#' @return An object of class 'power_model'.
#' @import numDeriv
#' @export
fit_model <- function(dat, alpha = 0.05, test = "onesided",
                      start = c(logbeta0 = 0, logtheta = 0, logr = 0))
{
  opt <- nlminb(start = start, objective = negLL,
                y = dat$y, m1 = dat$m1, m2 = dat$m2,
                alpha = alpha, test = test)
  opt_phi <- opt$par

  H <- hess_fun(opt_phi = opt_phi, dat$y, dat$m1, dat$m2, alpha, test)
  V_mod <- solve(H)

  S <- score_mat(opt_phi = opt_phi, dat$y, dat$m1, dat$m2, alpha, test)
  B <-   t(S) %*% S
  V_rob <- V_mod %*% B %*% V_mod # sandwich estimator for SEs

  out <-  list(opt_phi = opt_phi,
         V_mod   = V_mod,
         V_rob   = V_rob,
         conv    = opt$convergence,
         ll      = -opt$objective,
         alpha = alpha,
         test = test)
  class(out) <- "power_model"
  return(out)
}

#' AIC method for power_model objects
#'
#' @param object A 'power_model' object
#' @param ... Additional arguments (ignored)
#' @return AIC value
#' @export
AIC.power_model <- function(object, ...) {
  k <- length(object$opt_phi)
  -2 * object$ll + 2 * k
}


#' @export
negLL <- function(phi, y, m1, m2, alpha = 0.05, test = "onesided",
                  parameterization = "natural",
                  transform = TRUE) {
  pars <- stabilityTransformation(phi, transform = transform)
  P <- power_fun(m1 = m1, m2 = m2,
                 beta0 = pars$beta0, beta1 = pars$beta1, beta2 = pars$beta2,
                 alpha = alpha, test = test, parameterization = parameterization)
  # numeric stability
  P <- pmin(pmax(P, 1e-12), 1 - 1e-12)
  -sum(dbinom(y, 1, P, log = TRUE))
}

#' Stability transformation for optimization parameters
#'
#' @param phi numeric vector (length 3). If transform=TRUE, phi are log-params.
#' @param transform logical. If TRUE: beta0=-exp(phi1), beta1=exp(phi2), beta2=exp(phi3).
#'                  If FALSE: beta0=-phi1, beta1=phi2, beta2=phi3 (natural scale).
#' @return list(beta0=..., beta1=..., beta2=...)
#' @export
stabilityTransformation <- function(phi, transform = TRUE) {
  stopifnot(length(phi) >= 3)
  if (!transform) {
    return(list(beta0 = -phi[1], beta1 = phi[2], beta2 = phi[3]))
  } else {
    return(list(beta0 = -exp(phi[1]), beta1 = exp(phi[2]), beta2 = exp(phi[3])))
  }
}

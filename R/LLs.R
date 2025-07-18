#' @export
negLL <- function(phi, y, m1, m2, alpha = 0.05, test = "onesided") {
  beta0 <- -exp(phi[1]); theta <- exp(phi[2]);   r <- exp(phi[3])
  P <- power_fun(m1, m2, beta0, theta, r, alpha, test)
  ## for numeric stability
  P  <- pmin(pmax(P, 1e-12), 1 - 1e-12)
  -sum(dbinom(y, 1, P, log = TRUE))
}

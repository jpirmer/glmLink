#' Compute power for cross-classified LMM
#'
#' Calculates analytic power based on freely estimated intercept, adapted from
#' Jiang et al. (2024).
#'
#' @param m1 Number of level-1 units per cluster.
#' @param m2 Number of level-2 units per cluster.
#' @param beta0 beta0. This maps to the critical z-value (e.g., -1.96), negative for left tail.
#' @param beta1 beta1. This maps to the first random variance if \code{parameterization} is \code{"natural"}.
#' @param beta2 beta2. This maps to the second random variance if \code{parameterization} is \code{"natural"}.
#' @param alpha Significance level.
#' @param test 'onesided' or 'twosided'.
#' @param parameterization If 'natural', then beta1 and beta2 correspond to first and second random variance. Might be numerically unstable.
#' @return Power (probability of rejecting null).
#' @export
power_fun <- function(m1, m2, beta0, beta1, beta2, alpha = 0.05,
                      test = "onesided",
                      parameterization = "natural") {

  if(parameterization != "natural")
  {
    sqrtM1M2 <- beta1 * sqrt(m1 * m2 / (beta2 * m2 + m1))
  }else{
    sqrtM1M2 <- sqrt(1 / (beta1/m1 + beta2/m2))
  }
  zcrit <- beta0
  if(test == "onesided")
  {
    pnorm(zcrit + sqrtM1M2)
  }else if(test == "twosided"){
    1 - pnorm(-zcrit - sqrtM1M2) + pnorm(zcrit - sqrtM1M2)
  }else{
    stop("Wrong argument to in 'test'!")
  }
}

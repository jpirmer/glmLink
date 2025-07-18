#' Compute power for cross-classified LMM
#'
#' Calculates analytic power based on freely estimated intercept, adapted from
#' Jiang et al. (2024).
#'
#' @param m1 Number of level-1 units per cluster.
#' @param m2 Number of level-2 units per cluster.
#' @param beta0 Critical z-value (e.g., -1.96), negative for left tail.
#' @param theta True effect size parameter.
#' @param r Intraclass correlation.
#' @param alpha Significance level.
#' @param test 'onesided' or 'twosided'.
#' @return Power (probability of rejecting null).
#' @export
power_fun <- function(m1, m2, beta0, theta, r, alpha = 0.05,
                      test = "onesided") {
  delta <- theta * sqrt(m1 * m2 / (r * m2 + m1))
  zcrit <- beta0
  if(test == "onesided")
  {
    pnorm(zcrit + delta)
  }else if(test == "twosided"){
    1 - pnorm(-zcrit - delta) + pnorm(zcrit - delta)
  }else{
    stop("Wrong argument to in 'test'!")
  }
}

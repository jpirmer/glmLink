#' Delta-method confidence intervals on power predictions
#'
#' @param object A 'power_model' fit.
#' @param newdata Data frame with m1, m2.
#' @param level Confidence level.
#' @param use.robust Use robust variance.
#' @return Data frame with estimates and CIs.
#' @export

predict.power_model <- function(object, newdata,
                                level = 0.95,
                           use.robust = FALSE)
{
  # start
  M1M2_grid_new <- newdata
  opt_phi <- object$opt_phi
  if(use.robust){
    cov_logphi <-  object$V_rob
  }else{
    cov_logphi <-  object$V_mod
  }
  alpha <- object$alpha
  test <- object$test


  beta0 <- -exp(opt_phi[1]); theta <- exp(opt_phi[2]);   r <- exp(opt_phi[3])
  m1_new <- M1M2_grid_new$m1;   m2_new <- M1M2_grid_new$m2; N <- nrow(M1M2_grid_new)
  Power_hat <- power_fun(m1_new, m2_new, beta0, theta, r, alpha, test)

  ## Gradient dP/dphi (Log-Parameter)
  g <- t(sapply(seq_along(m1_new), FUN = function(i){
    numDeriv::grad(func = function(ph) {
      bet <- -exp(ph[1]); th <- exp(ph[2]); rr <- exp(ph[3])
      power_fun(m1 = m1_new[i], m2 = m2_new[i], beta0 = bet, theta = th, r = rr,
                alpha = alpha, test = test)
    }, x = opt_phi)}))

  var_p <- sapply(seq_along(m1_new), FUN = function(i){
    g[i,,drop=FALSE] %*% cov_logphi %*% t(g[i,,drop=FALSE])})
  se_p  <- sqrt(var_p)

  z <- qnorm(1 - (1 - level)/2)
  ci <- cbind(Power_hat, Power_hat) + cbind(rep(-1, N), rep(1,N)) * z * se_p
  colnames(ci) <- c("cilb", "ciub")
  return(list(Power_hat = Power_hat, se = se_p, ci = ci))
}

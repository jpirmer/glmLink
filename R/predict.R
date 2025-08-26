#' Predict power with Delta-method CIs for power_model
#'
#' @param object A 'power_model' fit.
#' @param newdata Optional data.frame with columns m1, m2. If NULL, uses the
#'   training design stored in the fitted object.
#' @param level Confidence level for two-sided CIs (default 0.95).
#' @param use.robust If TRUE, use sandwich covariance; otherwise model-based.
#' @param transform Logical, whether to interpret \code{opt_phi} on the log-phi
#'   scale (as in the fit). Defaults to \code{object$transform}.
#' @return A list with components:
#'   \item{Power_hat}{Vector of predicted power at rows of \code{newdata}.}
#'   \item{se}{Vector of Delta-method standard errors.}
#'   \item{ci}{N x 2 matrix with columns \code{cilb}, \code{ciub}.}
#' @export
predict.power_model <- function(object, newdata = NULL,
                                level = 0.95,
                                use.robust = FALSE,
                                transform = object$transform) {

  # 1) Choose design (m1, m2)
  if (is.null(newdata)) {
    M1M2_grid_new <- data.frame(m1 = object$m1, m2 = object$m2)
  } else {
    stopifnot(is.data.frame(newdata), all(c("m1","m2") %in% names(newdata)))
    M1M2_grid_new <- newdata
  }

  # 2) Extract parameters/covariance from fit
  opt_phi <- object$opt_phi
  cov_phi <- if (isTRUE(use.robust)) object$V_rob else object$V_mod
  alpha <- object$alpha
  test  <- object$test
  parameterization <- object$parameterization

  # 3) Map phi -> (beta0, beta1, beta2) on the natural scale
  pars <- stabilityTransformation(opt_phi, transform = transform)
  beta0 <- pars$beta0; beta1 <- pars$beta1; beta2 <- pars$beta2

  m1_new <- M1M2_grid_new$m1
  m2_new <- M1M2_grid_new$m2
  N <- length(m1_new)

  # 4) Point predictions
  Power_hat <- power_fun(m1 = m1_new, m2 = m2_new,
                         beta0 = beta0, beta1 = beta1, beta2 = beta2,
                         alpha = alpha, test = test,
                         parameterization = parameterization)

  # 5) Gradient dP/dphi at opt_phi for each row (size N x p)
  g <- t(vapply(seq_len(N), function(i) {
    numDeriv::grad(
      func = function(ph) {
        p <- stabilityTransformation(ph, transform = transform)
        power_fun(m1 = m1_new[i], m2 = m2_new[i],
                  beta0 = p$beta0, beta1 = p$beta1, beta2 = p$beta2,
                  alpha = alpha, test = test,
                  parameterization = parameterization)
      },
      x = opt_phi
    )
  }, numeric(length(opt_phi))))

  # 6) Delta-method variance per row: g_i %*% Cov(phi) %*% g_i^T
  var_p <- vapply(seq_len(N), function(i) {
    gi <- matrix(g[i, ], nrow = 1)
    as.numeric(gi %*% cov_phi %*% t(gi))
  }, numeric(1L))
  se_p <- sqrt(var_p)

  # 7) Two-sided normal CIs
  z <- stats::qnorm(1 - (1 - level) / 2)
  ci <- cbind(Power_hat - z * se_p,
              Power_hat + z * se_p)
  colnames(ci) <- c("cilb", "ciub")

  list(Power_hat = Power_hat, se = se_p, ci = ci)
}

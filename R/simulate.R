#' Simulate binary outcomes under power model
#'
#' @inheritParams power_fun
#' @param M1M2_grid Data frame with columns m1 and m2.
#' @param seed Random seed for reproducibility.
#' @return Data frame with simulated y and true power.
#' @export
sim_data <- function(M1M2_grid, beta0_true, theta_true, r_true,
                     alpha = 0.05, seed = NULL, test = "onesided") {
  if (!is.null(seed)) set.seed(seed)
  m1 <- M1M2_grid$m1; m2 <- M1M2_grid$m2; N <- nrow(M1M2_grid)
  P <- power_fun(m1, m2, beta0_true, theta_true, r_true, alpha, test)
  y <- rbinom(N, size = 1, prob = P)
  data.frame(y = y, m1 = m1, m2 = m2, power_true = P)
}

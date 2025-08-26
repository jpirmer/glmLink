## Hessian at optimum_phi (dimension = length(opt_phi) × length(opt_phi))
hess_fun <- function(opt_phi, y, m1, m2, alpha = 0.05, test = "onesided",
                     parameterization = "natural",
                     transform = TRUE) {
  numDeriv::hessian(negLL, x = opt_phi,
                    y = y, m1 = m1, m2 = m2,
                    alpha = alpha, test = test,
                    parameterization = parameterization,
                    transform = transform)
}

## Score matrix (one row per observation)
score_mat <- function(opt_phi, y, m1, m2, alpha = 0.05,
                      test = "onesided", parameterization = "natural",
                      transform = TRUE) {
  rows <- lapply(seq_along(y), function(i)
    numDeriv::grad(negLL, x = opt_phi,
                   y = y[i], m1 = m1[i], m2 = m2[i],
                   alpha = alpha, test = test,
                   parameterization = parameterization,
                   transform = transform))
  do.call(rbind, rows)
}

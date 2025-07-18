## Hessian an optimum_phi (2 × 2)
hess_fun <- function(opt_phi, y, m1, m2, alpha = 0.05, test = "onesided") {
  numDeriv::hessian(negLL, x = opt_phi, y = y, m1 = m1, m2 = m2,
                    alpha = alpha, test = test)
}

## Score function for sandwich estimator of SEs
score_mat <- function(opt_phi, y, m1, m2, alpha = 0.05, test = "onesided") {
  t(sapply(seq_along(y), function(i)
    numDeriv::grad(negLL, x = opt_phi,
                   y = y[i], m1 = m1[i], m2 = m2[i],
                   alpha = alpha, test = test)))
}

# glmLink

## Introduction

This package collects link functions needed for simulation based power analysis.

---

## Installation

You can install `glmLink` from GitHub using the `remotes` package:

```r
# Install remotes package if not already installed
install.packages("remotes") 

# Install glmLink package
remotes::install_github("jpirmer/glmLink")

# Load glmLink package
library(glmLink)
```

---

## Example

```r
Reps <- 10
M1M2_grid_new <- M1M2_grid <- expand.grid("m1" = rep(round(seq(20, 300, length.out = 12)), Reps), 
                         "m2" = rep(round(seq(20, 300, length.out = 12)), 1))
data <- sim_data(M1M2_grid = M1M2_grid, beta0_true = -1.96,
                 theta_true = .3, r_true = .5, seed = 2025)
head(data)
fit_new <- fit_model(dat = data)
fit_glm <- glm(y ~ sqrt((m1*m2) / (m1 + m2)), data = data)
AIC(fit_new) # better
AIC(fit_glm)


pred_new <- predict(fit_new, newdata = M1M2_grid_new, use.robust = FALSE)$Power_hat
pred_glm <- predict(fit_glm, newdata = M1M2_grid_new, type = "response")

power_true <- power_fun(m1 = M1M2_grid_new$m1, 
                        m2 = M1M2_grid_new$m2, beta0 = -1.96, 
                        theta = .3, r = .5)

# RMSE:
sqrt(mean((pred_new - power_true)^2)) #better
sqrt(mean((pred_glm - power_true)^2))
```

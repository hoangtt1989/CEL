rm(list = ls())

pacman::p_load(tidyverse, MASS)

source('mirror_descent.R')
source('lm_admm.R')
source('lm_admm_outlier.R')
# source('lm_admm_fixed.R')
source('outlier_loop.R')
source('owenEL.R')
source('lm_profiler.R')
source('lm_nested_concord.R')
source('owenEL_old.R')
source('damped_newton.R')
source('lbfgs_step.R')
source('grad_step.R')
source('wts_beta_fun.R')



set.seed(123)

n <- 50
p <- 5

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
beta <- as.vector(rnorm(p))
noise <- as.vector(rnorm(n))

y <- X %*% beta + noise
mean_init <- solve(t(X) %*% X, t(X) %*% y)

fix_index <- 1
alpha_add <- .2

beta_test <- mean_init
beta_test[fix_index] <- mean_init[fix_index] + alpha_add

R <- as.vector(y - X %*% beta_test)
Z_mat <- R * X

owen_soln <- emplik(Z_mat)


########add outliers
set.seed(123)
outlier_len <- 5

# noise[1:outlier_len] <- noise[1:outlier_len] + sign(noise[1:outlier_len]) * runif(outlier_len, min = 5, max = 7)
y2 <- y
y2[1:outlier_len] <- sign(y2[1:outlier_len]) * runif(outlier_len, min = 12, max = 14)

# y <- X %*% beta + noise
mean_init <- solve(t(X) %*% X, t(X) %*% y2)

head(y - X %*% mean_init)

fix_index <- 1
alpha_add <- .2

beta_test <- mean_init
beta_test[fix_index] <- mean_init[fix_index] + alpha_add

head(y2 - X %*% beta_test)

test_fit <- CEL_lm_admm_outlier(X, y2, beta_test, q = outlier_len, outer_tol_type = 'primal_dual', mirror_arg = list(line_search = T), outlier_loop_arg = list(line_search = T, maxit = 2000))
test_out_idx <- which(test_fit$delta_opt != 0)


test_fit$wts[1:3]
test_fit$wts[test_fit$delta_opt != 0]

R <- as.vector(y2 - X %*% beta_test)
Z_mat <- R * X


owen_out_fit <- emplik(Z_mat)
owen_out_fit$wts[1:outlier_len]
owen_out_fit$wts[test_fit$delta_opt != 0]


owen_soln$wts[1:outlier_len]
owen_soln$wts[test_fit$delta_opt != 0]



df <- data_frame(idx = 1:n, pod = test_fit$wts, owen = as.vector(owen_soln$wts), owen_out = as.vector(owen_out_fit$wts)) %>% 
  gather(key, value, -idx)

ggplot(df, aes(x = idx, y = value)) +
  geom_point(aes(color = key), alpha = .7) +
  geom_vline(xintercept = c(1:outlier_len, which(test_fit$delta_opt != 0)), linetype = 2, alpha = .3)



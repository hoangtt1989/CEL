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
p <- 2

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
beta <- as.vector(rnorm(p))
noise <- as.vector(rnorm(n))

y <- X %*% beta + noise
mean_init1 <- solve(t(X) %*% X, t(X) %*% y)

fix_index <- 1
alpha_add <- 1

beta_test <- mean_init1
beta_test[fix_index] <- mean_init1[fix_index] + alpha_add

R <- as.vector(y - X %*% beta_test)
Z_mat <- R * X

owen_soln <- emplik(Z_mat)


########add outliers
set.seed(123)
outlier_len <- 1


# noise[1:outlier_len] <- noise[1:outlier_len] + sign(noise[1:outlier_len]) * runif(outlier_len, min = 5, max = 7)
y2 <- y
y2[1:outlier_len] <- y2[1:outlier_len] + runif(outlier_len, min = 18, max = 20)

# y <- X %*% beta + noise
mean_init2 <- solve(t(X) %*% X, t(X) %*% y2)


fix_index <- 1
alpha_add <- 1

beta_test <- mean_init1
beta_test[fix_index] <- mean_init1[fix_index] + alpha_add

head(y - X %*% beta_test)

head(y2 - X %*% beta_test)

# ord_fit <- CEL_lm_admm(beta_test, X, y2, F_reparam = diag(p), s_hat_reparam = rep(0, p), outer_tol_type = 'primal_dual', vary_penalty = 'none', dual_step = .5)

test_fit1 <- CEL_lm_admm_outlier(X, y2, beta_test, q = 0, outer_tol_type = 'primal_dual', dual_step = 5, 
                                mirror_arg = list(line_search = T), outlier_loop_arg = list(line_search = T, maxit = 2000))

test_fit <- CEL_lm_admm_outlier(X, y2, beta_test, q = 1, outer_tol_type = 'primal_dual', dual_step = 5, 
                                mirror_arg = list(line_search = T), outlier_loop_arg = list(line_search = T, maxit = 2000))
test_out_idx <- which(test_fit$delta_opt != 0)

hat_mat <- X %*% solve(t(X) %*% X, t(X))
diag(hat_mat)[test_out_idx]

crossprod(test_fit$wts * X, as.vector(test_fit$y_adj - X %*% test_fit$beta_opt))

test_fit$wts[1:outlier_len]
test_fit$wts[test_fit$delta_opt != 0]

R <- as.vector(y2 - X %*% beta_test)
Z_mat <- R * X

plot(Z_mat[, 1], Z_mat[, 2])


Z_hat <- Z_mat %*% solve(t(Z_mat) %*% Z_mat, t(Z_mat))
# cor(t(Z_mat))
# Z_svd <- svd(Z_mat)
# Z_svd$d


owen_out_fit <- emplik(Z_mat)
owen_out_fit$wts[1:outlier_len]
owen_out_fit$wts[test_fit$delta_opt != 0]


owen_soln$wts[1:outlier_len]
owen_soln$wts[test_fit$delta_opt != 0]



df <- data_frame(idx = 1:n, admm = test_fit1$wts, pod = test_fit$wts, owen_out = as.vector(owen_out_fit$wts)) %>% 
  gather(key, value, -idx)

ggplot(df, aes(x = idx, y = value)) +
  geom_point(aes(color = key), alpha = .7) +
  geom_vline(xintercept = 1:outlier_len, linetype = 2, alpha = .3, color = 'black') +
  geom_vline(xintercept = which(test_fit$delta_opt != 0), linetype = 2, alpha = .3, color = 'magenta') +
  geom_hline(yintercept = 1/n, linetype = 2, alpha = .3)



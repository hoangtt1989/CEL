rm(list = ls())

pacman::p_load(tidyverse, MASS)

source('mirror_descent.R')
source('lm_admm.R')
source('lm_admm_outlier.R')
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
alpha_add <- .3

##set up the reparam
T_reparam <- t(as.matrix(rep(0, p)))
T_reparam[fix_index] <- 1
alpha_test <- mean_init[fix_index] + alpha_add
s_hat_reparam <- as.vector(rep(0, p))
s_hat_reparam[fix_index] <- alpha_test
# s_hat_reparam <- pinv(T_reparam) * alpha_test
F_reparam <- orth(nullspace(T_reparam))
##
##set up the test
beta_test <- mean_init
beta_test[fix_index] <- alpha_test

# F_reparam <- NULL
# s_hat_reparam <- NULL


# T_reparam <- diag(p)
# alpha_test <- mean_init
# alpha_test[fix_index] <- mean_init[fix_index] + .1
# s_hat_reparam <- pinv(T_reparam) %*% alpha_test
# F_reparam <- orth(nullspace(T_reparam))

test_fun <- CEL_lm_admm_outlier(X, y, beta_test = beta_test, F_reparam = F_reparam, s_hat_reparam = s_hat_reparam, 
                                q = 1, outer_tol_type = 'primal_dual', outer_maxit = 200, 
                                mirror_arg = list(line_search = T), dual_step = 20, outlier_loop_arg = list(outlier_eps = 1e-4))


test_fun_noout <- CEL_lm_admm(X, y, beta_test = beta_test, F_reparam = F_reparam, s_hat_reparam = s_hat_reparam, 
                                 outer_tol_type = 'primal_dual',
                                mirror_arg = list(line_search = T), dual_step = 20)

test_fun$beta_opt - test_fun_noout$beta_opt

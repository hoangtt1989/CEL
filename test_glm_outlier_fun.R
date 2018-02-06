rm(list = ls())

pacman::p_load(tidyverse, MASS, doMC)

source('mirror_descent.R')
source('lm_admm.R')
source('glm_admm_outlier.R')
source('owenEL.R')
source('lm_profiler.R')
source('glm_profiler.R')
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
alpha_add <- .1

##set up the reparam
T_reparam <- t(as.matrix(rep(0, p)))
T_reparam[fix_index] <- 1
alpha_test <- mean_init[fix_index] + alpha_add
##set up the test
beta_test <- mean_init
beta_test[fix_index] <- alpha_test

source('proximal_gradient_outlier.R')
test_fit <- CEL_glm_admm_outlier(beta_test, X, y, outer_tol_type = 'primal_dual', 
                                 prox_fun = prox_affine, 
                                 prox_fun_arg = list(T_reparam = T_reparam, alpha_reparam = alpha_test), mirror_arg = list(line_search = T, mirror_eps = 1e-5), 
                                 dual_step = 20, RB_tau = 1.5, RB_mu = 10, 
                                 prox_optim_arg = list(prox_eps = 1e-4, tol_type = 'fval', maxit = 2000))

plot(test_fit$wts)

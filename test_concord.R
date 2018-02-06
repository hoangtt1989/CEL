rm(list = ls())

pacman::p_load(tidyverse, MASS)

source('mirror_descent.R')
source('lm_admm.R')
source('owenEL.R')
source('lm_profiler.R')
source('owenEL_old.R')
source('lbfgs_step.R')
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


###no concord
source('lm_nested.R')
nest <- CEL_lm_nested(beta_test, X, y, F_reparam, s_hat_reparam)

###concord
source('lm_nested_concord.R')
source('damped_newton.R')
source('lbfgs_step.R')
source('grad_step.R')
nest_concord <- CEL_lm_nested(beta_test, X, y, F_reparam, s_hat_reparam)

nest_concord_lbfgs <- CEL_lm_nested(beta_test, X, y, F_reparam, s_hat_reparam, beta_opt_type = 'LBFGS')

nest_concord_grad <- CEL_lm_nested(beta_test, X, y, F_reparam, s_hat_reparam, beta_opt_type = 'gradient')

nest$outer_nits
nest_concord$outer_nits

nest$inner_nits
nest_concord$inner_nits

nest$outer_fval[nest$outer_nits] - nest_concord$outer_fval[nest_concord$outer_nits]

###admm
admm <- CEL_lm_admm(beta_test, X, y, F_reparam, s_hat_reparam)
admm$outer_nits
admm$logelr

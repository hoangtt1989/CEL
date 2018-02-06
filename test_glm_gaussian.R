rm(list = ls())

pacman::p_load(tidyverse, MASS)

source('mirror_descent.R')
source('lm_admm.R')
source('glm_admm.R')
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


lm_res1 <- CEL_lm_admm(beta_test, X, y, F_reparam, s_hat_reparam, beta_opt_type = 'closed_form')
glm_res1 <- CEL_glm_admm(beta_test, X, y, F_reparam, s_hat_reparam, beta_opt_type = 'closed_form', family = 'gaussian')


lm_res2 <- CEL_lm_admm(beta_test, X, y, F_reparam, s_hat_reparam, beta_opt_type = 'LBFGS', verbose = T)
source('glm_admm.R')
source('glm_beta_wrapper.R')
glm_res2 <- CEL_glm_admm(beta_test, X, y, F_reparam, s_hat_reparam, beta_opt_type = 'LBFGS', family = 'gaussian', verbose = T)

rm(list = ls())

pacman::p_load(tidyverse, MASS, doMC)

source('mirror_descent.R')
source('lm_admm.R')
source('glm_admm.R')
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


reparam_res <- CEL_glm_admm(beta_test, X, y, F_reparam, s_hat_reparam, family = 'gaussian')
prox_res <- CEL_glm_admm(beta_test, X, y, family = 'gaussian', prox_fun = prox_affine, prox_fun_arg = list(T_reparam = T_reparam, alpha_reparam = alpha_test))

prox_res$beta_opt - reparam_res$beta_opt


prox_ineq <- CEL_glm_admm(beta_test, X, y, family = 'gaussian', prox_fun = prox_inequality, prox_fun_arg = list(test_idx = 1:p, test_val = 0, dir = 'gt'))
prox_interval_res <- CEL_glm_admm(beta_test, X, y, family = 'gaussian', prox_fun = prox_interval, prox_fun_arg = list(test_idx = 1:p, left_val = 0, right_val = 1))



###########binomial
set.seed(123)

n <- 50
p <- 5

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
beta <- as.vector(rnorm(p))
noise <- as.vector(rnorm(n))

## generate probabilities
noisy_prob <- 1 / (1 + exp(-(X %*% beta + noise)))
y <- rbinom(n, 1, noisy_prob)
#######

# y <- X %*% beta + noise
# mean_init <- solve(t(X) %*% X, t(X) %*% y)
glm_fit <- glm.fit(X, y, family = binomial(link = 'logit'))
mean_init <- glm_fit$coefficients

fix_index <- 1
alpha_add <- .1

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

reparam_res <- CEL_glm_admm(beta_test, X, y, F_reparam, s_hat_reparam, family = 'binomial')
prox_res <- CEL_glm_admm(beta_test, X, y, family = 'binomial', prox_fun = prox_affine, prox_fun_arg = list(T_reparam = T_reparam, alpha_reparam = alpha_test))

prox_res$beta_opt - reparam_res$beta_opt


beta_test <- mean_init
prox_ineq <- CEL_glm_admm(beta_test, X, y, family = 'binomial', prox_fun = prox_inequality, prox_fun_arg = list(test_idx = 1:p, test_val = 0, dir = 'gt'))
prox_interval_res <- CEL_glm_admm(beta_test, X, y, family = 'binomial', prox_fun = prox_interval, prox_fun_arg = list(test_idx = 1:p, left_val = 0, right_val = 1))

##test boot
registerDoMC(parallel::detectCores())
prox_boot_res <- CEL_glm_boot(X, y, prox_interval, family = 'binomial', prox_fun_arg = list(test_idx = 1:p, left_val = 0, right_val = 1), B = 100, outer_tol_type = 'primal_dual', mirror_arg = list(mirror_eps = 1e-6, maxit = 2000))
prox_boot_ineq_res <- CEL_glm_boot(X, y, prox_inequality, family = 'binomial', #####fails!!!
                                   prox_fun_arg = list(test_idx = 1:p, test_val = 0, dir = 'gt'), 
                                   B = 100, outer_tol_type = 'primal_dual', dual_step = 20, wts_beta_rep = 3,
                                   outer_rel_eps = 1e-3, prox_optim_arg = list(prox_eps = 1e-5), mirror_arg = list(mirror_eps = 1e-6, maxit = 2000))
prox_boot_ineq_res <- CEL_glm_boot(X, y, prox_inequality, family = 'binomial',
                                   prox_fun_arg = list(test_idx = 1, test_val = 0, dir = 'gt'), 
                                   B = 100, outer_tol_type = 'primal_dual', dual_step = 20, wts_beta_rep = 3,
                                   outer_rel_eps = 1e-3, prox_optim_arg = list(prox_eps = 1e-5), mirror_arg = list(mirror_eps = 1e-6, maxit = 2000, line_search = T))
prox_boot_ineq_res2 <- CEL_glm_boot(X, y, prox_inequality, family = 'binomial', prox_fun_arg = list(test_idx = 1, test_val = .5, dir = 'gt'), prox_optim_arg = list(prox_eps = 1e-5, maxit = 2000),
                                    B = 100, outer_tol_type = 'primal_dual', dual_step = 20, wts_beta_rep = 2,
                                    outer_rel_eps = 1e-3, mirror_arg = list(mirror_eps = 1e-5, maxit = 2000, line_search = T))
##
###########

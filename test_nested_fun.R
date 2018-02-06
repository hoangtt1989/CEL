rm(list = ls())

pacman::p_load(tidyverse, MASS)

source('mirror_descent.R')
source('lm_admm.R')
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

# nest2 <- CEL_lm_nested2(beta_test, X, y, F_reparam, s_hat_reparam)
# nest2$time
# nest_lbfgs <- CEL_lm_nested(beta_test, X, y, F_reparam, s_hat_reparam, beta_opt_type = 'LBFGS')
nest_grad <- CEL_lm_nested(beta_test, X, y, F_reparam, s_hat_reparam, beta_opt_type = 'gradient')
nest1 <- CEL_lm_nested(beta_test, X, y, F_reparam, s_hat_reparam, beta_newton_arg = list(mod_hess = F))
nest_mod <- CEL_lm_nested(beta_test, X, y, F_reparam, s_hat_reparam, beta_newton_arg = list(mod_hess = T))
nest_inner_lbfgs <- CEL_lm_nested(beta_test, X, y, F_reparam, s_hat_reparam, beta_newton_arg = list(mod_hess = F), inner_opt_type = 'LBFGS')

admm1 <- CEL_lm_admm(beta_test, X, y, F_reparam, s_hat_reparam, beta_opt_type = 'closed_form', mirror_arg = list(vary_step = T))
# admm1_ls <- CEL_lm_admm(beta_test, X, y, F_reparam, s_hat_reparam, beta_opt_type = 'closed_form', mirror_arg = list(line_search = T))
admm1_RB2 <- CEL_lm_admm(beta_test, X, y, F_reparam, s_hat_reparam, beta_opt_type = 'closed_form', vary_penalty = 'RB2', RB2_ksi = .5)
admm_pd <- CEL_lm_admm(beta_test, X, y, F_reparam, s_hat_reparam, beta_opt_type = 'closed_form', outer_tol_type = 'primal_dual')
admm1_nvary <- CEL_lm_admm(beta_test, X, y, F_reparam, s_hat_reparam, beta_opt_type = 'closed_form', mirror_arg = list(vary_step = F))
admm2 <- CEL_lm_admm(beta_test, X, y, F_reparam, s_hat_reparam, beta_opt_type = 'closed_form', wts_beta_rep = 1, mirror_arg = list(line_search = T))
admm3 <- CEL_lm_admm2(beta_test, X, y, F_reparam, s_hat_reparam, beta_opt_type = 'closed_form', wts_beta_rep = 4)
admm4 <- CEL_lm_admm2(beta_test, X, y, F_reparam, s_hat_reparam, beta_opt_type = 'closed_form', wts_beta_rep = 1, random_order = F)
admm6 <- CEL_lm_admm2(beta_test, X, y, F_reparam, s_hat_reparam, beta_opt_type = 'closed_form', wts_beta_rep = 1, random_order = T)
admm5 <- CEL_lm_admm2(beta_test, X, y, F_reparam, s_hat_reparam, beta_opt_type = 'closed_form', wts_beta_rep = 4, random_order = T)
# admm$time


nest1$time
admm$time

# nest2$logelr
nest1$logelr
admm$logelr

# nest2$outer_nits
nest1$outer_nits

# nest2$inner_nits
nest1$inner_nits

# sum(nest2$beta_opt - nest1$beta_opt)

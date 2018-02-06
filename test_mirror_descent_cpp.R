rm(list = ls())

set.seed(123)

pacman::p_load(tidyverse, benchmark, microbenchmark, Rcpp, RcppEigen)

source('lm_admm.R')
source('lm_nested_concord.R')
source('lm_profiler.R')
source('mirror_descent.R')


set.seed(123)

n <- 100
p <- 5

beta_true <- c(5, 4, 4, 3, 2)

# alpha_seq <- seq(1e-3, .7, length.out = 10)

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
noise <- rnorm(n)

y <- X %*% beta_true + noise

beta_OLS <- solve(crossprod(X), crossprod(X, y))
beta_test <- beta_OLS
beta_test[1] <- beta_OLS[1] + .1

gamma <- rnorm(p)

dual_step <- 2

fun_md <- fun_lm_AL2
grad_md <- grad_lm_AL2

R <- as.vector(y - X %*% beta_test)
Xtr_R <- t(R * X)
wts <- rep(1, n)

A_max <- dual_step * max(abs(crossprod(Xtr_R))) * n ##t(Xtr_R) %*% Xtr_R
mirror_step <- 1 / (A_max + max(abs(crossprod(gamma, Xtr_R))))


cpp_res <- mirror_descent_lm_AL_cpp(wts, gamma, X, R, dual_step, mirror_step = mirror_step, line_search = T, ls_beta = .5)
fun_arg <- list(gamma = gamma, X = X, R = R, dual_step = dual_step)
reg_res <- mirror_descent(fun_md, grad_md, wts, fun_arg = fun_arg, grad_arg = fun_arg, mirror_step = mirror_step, line_search = T)
max(abs(cpp_res$wts - reg_res$wts))

# tst1 <- runif(p)
# tst2 <- runif(p)
# bregman_entropy(tst1, tst2)
# bregman_entropy_cpp(tst1, tst2)
# identical(as.numeric(bregman_entropy(tst1, tst2)), bregman_entropy_cpp(tst1, tst2))


microbenchmark(mirror_descent_lm_AL_cpp(wts, gamma, X, R, dual_step, mirror_step = mirror_step),
               mirror_descent(fun_md, grad_md, wts, fun_arg = fun_arg, grad_arg = fun_arg, mirror_step = mirror_step))

# 
# set.seed(123)
# 
# n <- 100
# p <- 5
# 
# X <- matrix(rnorm(n * p), nrow = n, ncol = p)
# beta <- as.vector(rnorm(p))
# noise <- as.vector(rnorm(n))
# 
# y <- X %*% beta + noise
# mean_init <- solve(t(X) %*% X, t(X) %*% y)
# 
# fix_index <- 1
# alpha_add <- .3
# 
# ##set up the reparam
# T_reparam <- t(as.matrix(rep(0, p)))
# T_reparam[fix_index] <- 1
# alpha_test <- mean_init[fix_index] + alpha_add
# s_hat_reparam <- as.vector(rep(0, p))
# s_hat_reparam[fix_index] <- alpha_test
# # s_hat_reparam <- pinv(T_reparam) * alpha_test
# F_reparam <- orth(nullspace(T_reparam))
# ##
# ##set up the test
# beta_test <- mean_init
# beta_test[fix_index] <- alpha_test
# 
# admm_reg <- CEL_lm_admm(beta_test, X, y, F_reparam, s_hat_reparam, beta_opt_type = 'closed_form')
# admm_cpp <- CEL_lm_admm_cpp(beta_test, X, y, F_reparam, s_hat_reparam, beta_opt_type = 'closed_form')
# 
# max(abs(admm_reg$wts - admm_cpp$wts))
# max(abs(admm_reg$logelr - admm_cpp$logelr))
# max(abs(admm_reg$beta_opt - admm_cpp$beta_opt))
# 
# microbenchmark(CEL_lm_admm(beta_test, X, y, F_reparam, s_hat_reparam, beta_opt_type = 'closed_form'), CEL_lm_admm_cpp(beta_test, X, y, F_reparam, s_hat_reparam, beta_opt_type = 'closed_form'))

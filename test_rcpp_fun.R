rm(list = ls())

pacman::p_load(MASS, Rcpp, RcppEigen, microbenchmark, benchmark)

source('mirror_descent.R')
source('lm_admm.R')
source('owenEL.R')
source('lm_profiler.R')
source('lm_nested.R')
source('owenEL_old.R')
source('damped_newton.R')
source('lbfgs_step.R')
source('wts_beta_fun.R')

sourceCpp('rcpp_fun.cpp')

set.seed(123)

n <- 1000
p <- 50

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
beta <- as.vector(rnorm(p))
noise <- as.vector(rnorm(n))

y <- X %*% beta + noise

mean_init <- solve(t(X) %*% X, t(X) %*% y)

fix_index <- 1
alpha_add <- .9

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

R <- as.vector(y - X %*% beta_test)

# microbenchmark(crossprod(X), AtA(X))

wts <- rep(1, n)
gamma <- rep(0, p)
dual_step <- 1

tst <- fun_lm_AL_cpp(wts, gamma, X, R, 1)
all.equal(tst, as.numeric(fun_lm_AL2(wts, gamma, X, R, 1)))

tst_grad <- grad_lm_AL_cpp(wts, gamma, X, R, 1)
all.equal(tst_grad, as.numeric(grad_lm_AL2(wts, gamma, X, R, 1)))

microbenchmark(fun_lm_AL_cpp(wts, gamma, X, R, 1), fun_lm_AL2(wts, gamma, X, R, 1))

microbenchmark(grad_lm_AL_cpp(wts, gamma, X, R, 1), grad_lm_AL2(wts, gamma, X, R, 1))
# microbenchmark(fun_lm_AL_cpp(X, R), R * X)

tmp_fun <- function(wts, gamma, y, X, F_reparam, s_hat_reparam, dual_step) {
  Xtr_W_X <- crossprod(X, as.vector(wts) * X)
  z_hat_new <- solve(crossprod(F_reparam, crossprod(Xtr_W_X) %*% F_reparam), 
                     crossprod(F_reparam, Xtr_W_X) %*% (crossprod(X, wts * y) + gamma / dual_step - Xtr_W_X %*% s_hat_reparam))
}

microbenchmark(closed_beta_lm_AL_cpp(wts, gamma, y, X, F_reparam, s_hat_reparam, dual_step), tmp_fun(wts, gamma, y, X, F_reparam, s_hat_reparam, dual_step))
tst <- tmp_fun(wts, gamma, y, X, F_reparam, s_hat_reparam, dual_step)

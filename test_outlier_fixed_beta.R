rm(list = ls())

pacman::p_load(tidyverse, MASS)

source('mirror_descent.R')
source('lm_admm.R')
source('lm_admm_fixed.R')
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
alpha_add <- .1


beta_test <- mean_init
beta_test[fix_index] <- mean_init[fix_index] + alpha_add

res <- EL_lm_AL_outlier(beta_test, X, y, q = 3, dual_step = 3)

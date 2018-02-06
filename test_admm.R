rm(list = ls())

pacman::p_load(MASS)

source('mirror_descent.R')
source('lm_admm.R')
source('lm_nested.R')
source('lm_profiler.R')
source('owenEL.R')

set.seed(123)

n <- 1000
p <- 500

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
beta <- as.vector(rnorm(p))
noise <- as.vector(rnorm(n))

y <- X %*% beta + noise
beta_OLS <- solve(crossprod(X), crossprod(X, y))

prof_rand_lbfgs <- CEL_lm_beta_fixed(.1, 1, X, y, verbose = T, outer_eps = 1e-5, random_order = T, mirror_arg = list(mirror_eps = 1e-7, maxit = 2000), beta_opt_type = 'LBFGS', lbfgs_arg = list(invisible = 1))
prof_admm <- CEL_lm_beta_fixed(.1, 1, X, y, verbose = T, outer_eps = 1e-5, mirror_arg = list(mirror_eps = 1e-7, maxit = 2000), beta_opt_type = 'LBFGS', lbfgs_arg = list(invisible = 1))
prof_admm_closed <- CEL_lm_beta_fixed(.1, 1, X, y, verbose = T, outer_eps = 1e-5, mirror_arg = list(mirror_eps = 1e-7, maxit = 2000), beta_opt_type = 'closed_form')


prof_nested <- CEL_lm_beta_fixed(.1, 1, X, y, verbose = T, outer_eps = 1e-5, algorithm = 'nested')

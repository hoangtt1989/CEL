rm(list = ls())

pacman::p_load(MASS)

source('mirror_descent.R')
source('lm_admm.R')
source('owenEL.R')
source('lm_profiler.R')
source('lm_nested.R')
source('owenEL_old.R')

set.seed(123)

n <- 50
p <- 5

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
beta <- as.vector(rnorm(p))
noise <- as.vector(rnorm(n))

y <- X %*% beta + noise

nested_res <- emplik_lm_beta_fixed(1, 5, X, y, algorithm = 'nested')

admm_res <- emplik_lm_beta_fixed(.3, 5, X, y, algorithm = 'ADMM', verbose = T, mirror_arg = list(mirror_eps = 1e-8, max_iter = 1000, tol_type = 'logelr'))
admm_res$inner_nits[1:5]

admm_res$logelr - nested_res$logelr

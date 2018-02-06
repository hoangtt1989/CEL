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


ADMM_res <- emplik_lm_pred_fixed(5, 1, X, y, algorithm = 'ADMM')

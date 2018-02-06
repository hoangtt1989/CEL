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
p <- 10

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
beta <- as.vector(rnorm(p))
noise <- as.vector(rnorm(n))

y <- X %*% beta + noise

idx <- 2


nested_int <- emplik_lm_beta_profiler(idx, X, y, upper_divide = 4, algorithm = 'nested', outer_eps = 1e-6)
nested_int


ADMM_int <- emplik_lm_beta_profiler(idx, X, y, upper_divide = 4, algorithm = 'ADMM', outer_eps = 1e-8)
ADMM_int


rm(list = ls())

set.seed(123)

pacman::p_load(tidyverse, doMC)

source('lm_admm.R')
source('lm_nested_concord.R')
source('lm_profiler.R')


n <- 100
p <- 5

beta_true <- c(5, 4, 4, 3, 2)

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
noise <- rnorm(n)

y <- X %*% beta_true + noise

fix_index <- 1

registerDoMC(cores = parallel::detectCores())
res_ADMM <- CEL_lm_beta_boot(fix_index, X, y, B = 200, algorithm = 'ADMM', parallel = T, outer_tol_type = 'primal_dual', dual_step = 10, RB_tau = 1.5)
res_nested <- CEL_lm_beta_boot(fix_index, X, y, B = 1000, algorithm = 'nested', parallel = T)

test_ADMM <- CEL_lm_beta_fixed(5, 1, X, y, algorithm = 'ADMM')

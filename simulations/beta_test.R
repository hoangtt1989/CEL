rm(list = ls())

set.seed(123)

pacman::p_load(tidyverse)

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
beta_fix_val <- beta_true[fix_index]


res_ADMM <- CEL_lm_beta_fixed(beta_fix_val, fix_index, X, y, algorithm = 'ADMM')
res_nested <- CEL_lm_beta_fixed(beta_fix_val, fix_index, X, y, algorithm = 'nested')

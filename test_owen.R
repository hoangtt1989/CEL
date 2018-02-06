rm(list = ls())

pacman::p_load(MASS)

source('lm_profiler.R')
source('lm_nested.R')
source('owenEL.R')
source('owenEL_old.R')

set.seed(123)

n <- 50
p <- 5

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
beta <- as.vector(rnorm(p))
noise <- as.vector(rnorm(n))

y <- X %*% beta + noise
beta_OLS <- solve(crossprod(X), crossprod(X, y))


##initial values
beta_test <- beta_OLS
beta_test[1] <- beta_OLS[1] + .1
R <- as.vector(y - X %*% beta_test)
score_eq <- R * X

owen_old <- elm(score_eq)
# owen_old$lambda

owen_cc <- emplik(score_eq)
max(abs(owen_cc$lam - owen_old$lambda))

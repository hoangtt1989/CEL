rm(list = ls())

pacman::p_load(MASS)

source('lm_admm.R')
source('owenEL_old.R')
source('mirror_descent.R')

set.seed(123)

n <- 1000
p <- 200

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
beta <- as.vector(rnorm(p))
noise <- as.vector(rnorm(n))

y <- X %*% beta + noise
# beta_OLS <- ginv(t(X) %*% X) %*% (t(X) %*% y)
beta_OLS <- solve(t(X) %*% X, t(X) %*% y)

##initial values
beta_test <- beta_OLS
beta_test[1] <- beta_OLS[1] + .4

res_AL <- EL_lm_AL(beta_test, X, y, outer_tol_type = 'primal_dual', dual_step = .7, vary_penalty = T,
                   mirror_arg = list(tol_type = 'fval', maxit = 1000, line_search = F), verbose = T)
res_AL$time
res_AL$outer_nits
res_AL$logelr
median(res_AL$inner_nits)
plot(res_AL$outer_fval)

R_curr <- diag(as.vector(y - X %*% beta_test))
Z <- R_curr %*% X
system.time(res_owen <- elm(Z))
res_owen$logelr

res_AL$logelr - res_owen$logelr

rm(list = ls())

pacman::p_load(MASS)

source('mirror_descent.R')
source('owenEL.R')

set.seed(123)

n <- 50
p <- 5

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
beta <- as.vector(rnorm(p))
noise <- as.vector(rnorm(n))

y <- X %*% beta + noise
beta_OLS <- ginv(t(X) %*% X) %*% (t(X) %*% y)
beta_OLS2 <- solve(t(X) %*% X, t(X) %*% y)

##initial values
beta_test <- beta_OLS
beta_test[1] <- beta_OLS[1] + .1
R_curr <- diag(as.vector(y - X %*% beta_test))
wts_curr <- as.vector(rep(1, n))
gamma_curr <- as.vector(rep(0, p))

outer_tol_type <- 'fval'
outer_stop <- 5000
outer_converged <- F
outer_eps <- 1e-8
dual_step <- .03

Xtr_R <- t(X) %*% R_curr
A_mat <- t(Xtr_R) %*% Xtr_R
A_max <- n * dual_step * max(abs(A_mat))

gammatr_Xtr_R <- t(gamma_curr) %*% Xtr_R
R_X_gamma <- R_curr %*% X %*% gamma_curr
outer_fval_curr <- do.call(fun_lm_AL, list(wts = wts_curr, gammatr_Xtr_R, Xtr_R, dual_step))
outer_fval <- NULL
inner_nits <- NULL
inner_tol <- NULL

f_arg <- list(gammatr_Xtr_R = gammatr_Xtr_R, Xtr_R = Xtr_R, dual_step = dual_step)
g_arg <- list(R_X_gamma = R_X_gamma, Xtr_R = Xtr_R, dual_step = dual_step)

##compare with owen
Z <- R_curr %*% X
owen_res <- emplik(Z)
##

for (j in 1:outer_stop) {
  
  ##mirror descent
  mirror_step <- 1 / (A_max + max(abs(R_curr %*% (X %*% gamma_curr))))
  mirror_res <- mirror_descent(fun_lm_AL, grad_lm_AL, wts_curr, f_arg, g_arg, mirror_step = mirror_step, tol_type = 'logelr', max_iter = 1000, mirror_eps = 1e-7)
  wts_curr <- mirror_res$wts
  inner_nits[j] <- mirror_res$iter
  inner_tol[j] <- mirror_res$tol
  ##
  
  ##gamma update
  gamma_new <- gamma_curr + dual_step * Xtr_R %*% wts_curr
  gamma_diff <- max(abs(gamma_new - gamma_curr)) / max(abs(gamma_curr))
  gamma_curr <- gamma_new
  ##
  ##function value update
  g_arg$R_X_gamma <- t(Xtr_R) %*% gamma_curr
  gammatr_Xtr_R <- t(gamma_curr) %*% Xtr_R
  f_arg$gammatr_Xtr_R <- gammatr_Xtr_R
  outer_fval[j] <- do.call(fun_lm_AL, c(list(wts = wts_curr), f_arg))
  outer_fval_new <- outer_fval[j]
  outer_fval_diff <- abs(outer_fval_new - outer_fval_curr) / abs(outer_fval_curr)
  outer_fval_curr <- outer_fval_new
  ##
  ##primal residual
  primal_resid <- sum((Xtr_R %*% wts_curr)^2)
  ##
  ##stopping
  if (outer_tol_type == 'fval') {
    outer_tol <- outer_fval_diff
  }
  if (outer_tol < outer_eps) {
    outer_converged <- T
    break
  }
}

my_logelr <- sum(log(wts_curr))

my_logelr - owen_res$logelr
my_logelr > owen_res$logelr

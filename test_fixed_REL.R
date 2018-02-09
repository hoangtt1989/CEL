rm(list = ls())

pacman::p_load(tidyverse, MASS)

source('mirror_descent.R')
source('lm_admm.R')
source('lm_admm_REL.R')
# source('lm_admm_fixed.R')
source('REL_gradient.R')
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
p <- 2

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
beta <- as.vector(rnorm(p))
noise <- as.vector(rnorm(n))


#####clean data
y <- X %*% beta + noise
mean_init <- solve(t(X) %*% X, t(X) %*% y)

fix_index <- 1
alpha_add <- .5

beta_test <- mean_init
beta_test[fix_index] <- mean_init[fix_index] + alpha_add

R <- as.vector(y - X %*% beta_test)
Z_mat <- R * X

owen_soln <- emplik(Z_mat)

REL_soln <- CEL_lm_admm_outlier(X, y, beta_test = beta_test, q = 3, 
                                outer_tol_type = 'primal_dual', 
                                mirror_arg = list(line_search = T), 
                                outer_rel_eps = 1e-4, dual_step = 5, outer_maxit = 100)
#####

#####outliers
outlier_len <- 3
y2 <- y
y2[1:outlier_len] <- y2[1:outlier_len] + runif(outlier_len, min = 12, max = 14)

R_out <- as.vector(y2 - X %*% beta_test)
Z_mat_out <- R_out * X

set.seed(123)
REL_out_soln <- CEL_lm_admm_outlier(X, y2, beta_test = beta_test, q = outlier_len, outer_tol_type = 'primal_dual', 
                                    mirror_arg = list(line_search = T), outlier_loop_arg = list(line_search = T))
REL_out_soln$delta_opt
REL_out_soln$primal_resid
REL_out_soln$dual_resid
sqrt(sum(crossprod(as.vector(REL_out_soln$wts * n - REL_out_soln$delta_opt) * X, y2 - X %*% beta_test)^2))
crossprod(as.vector(1 - REL_out_soln$delta_opt) * X, y2 - X %*% beta_test)
REL_out_idx <- which(REL_out_soln$delta_opt != 0)
(REL_out_soln$wts * n - REL_out_soln$delta_opt)[REL_out_idx]

set.seed(123)
REL_out_soln_pos <- CEL_lm_admm_outlier(X, y2, beta_test = beta_test, q = outlier_len, outer_tol_type = 'primal_dual', outlier_loop_arg = list(delta_pos = T))
crossprod(as.vector(1 - REL_out_soln_pos$delta_opt) * X, y2 - X %*% beta_test)
REL_out_soln_pos$delta_opt
REL_out_pos_idx <- which(REL_out_soln_pos$delta_opt != 0)
(REL_out_soln_pos$wts * n - REL_out_soln_pos$delta_opt)[REL_out_pos_idx]

##delta plus
set.seed(123)
REL_out_soln_plus <- CEL_lm_admm_outlier(X, y2, beta_test = beta_test, q = outlier_len, outer_tol_type = 'primal_dual', delta_plus = T)
REL_out_plus_idx <- which(REL_out_soln_plus$delta_opt != 0)
(REL_out_soln_plus$wts * n - REL_out_soln_plus$delta_opt)[REL_out_plus_idx]
##


owen_out_soln <- emplik(Z_mat_out)
#####

#####

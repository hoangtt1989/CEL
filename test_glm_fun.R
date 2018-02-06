rm(list = ls())

pacman::p_load(tidyverse, MASS, doMC)

source('mirror_descent.R')
source('glm_admm.R')
source('glm_profiler.R')
source('owenEL.R')
source('lm_profiler.R')
source('lm_nested_concord.R')
source('owenEL_old.R')
source('damped_newton.R')
source('lbfgs_step.R')
source('grad_step.R')


#####binomial
set.seed(123)

n <- 50
p <- 5

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
beta <- as.vector(rnorm(p))
noise <- as.vector(rnorm(n))

## generate probabilities
noisy_prob <- 1 / (1 + exp(-(X %*% beta + noise)))
y <- rbinom(n, 1, noisy_prob)
#######

# y <- X %*% beta + noise
# mean_init <- solve(t(X) %*% X, t(X) %*% y)
glm_fit <- glm.fit(X, y, family = binomial(link = 'logit'))
mean_init <- glm_fit$coefficients

fix_index <- 1
alpha_add <- .1

##set up the reparam
T_reparam <- t(as.matrix(rep(0, p)))
T_reparam[fix_index] <- 1
alpha_test <- mean_init[fix_index] + alpha_add
s_hat_reparam <- as.vector(rep(0, p))
s_hat_reparam[fix_index] <- alpha_test
# s_hat_reparam <- pinv(T_reparam) * alpha_test
F_reparam <- orth(nullspace(T_reparam))
##
##set up the test
beta_test <- mean_init
beta_test[fix_index] <- alpha_test

res <- CEL_glm_admm(beta_test, X, y, F_reparam, s_hat_reparam, family = 'binomial', beta_opt_type = 'LBFGS')

###test bootstrap
registerDoMC(parallel::detectCores())
boot_res <- CEL_glm_beta_boot(1, X, y, family = 'binomial', B = 200, 
                              outer_tol_type = 'primal_dual', dual_step = 20, 
                              wts_beta_rep = 2, outer_rel_eps = 1e-3, mirror_arg = list(mirror_eps = 1e-6, maxit = 2000))
boot_res$boot_stat
boot_res$boot_vec
####
#########

#########poisson
set.seed(123)

n <- 50
p <- 5

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
beta <- as.vector(rnorm(p))
noise <- as.vector(rnorm(n))

## generate means
noisy_mean <- exp(X %*% beta + noise)
# noisy_prob <- 1 / (1 + exp(-(X %*% beta + noise)))
# y <- rbinom(n, 1, noisy_prob)
y <- rpois(n, noisy_mean)
#######

# y <- X %*% beta + noise
# mean_init <- solve(t(X) %*% X, t(X) %*% y)
glm_fit <- glm.fit(X, y, family = poisson(link = 'log'))
mean_init <- glm_fit$coefficients

fix_index <- 1
alpha_add <- .1

##set up the reparam
T_reparam <- t(as.matrix(rep(0, p)))
T_reparam[fix_index] <- 1
alpha_test <- mean_init[fix_index] + alpha_add
s_hat_reparam <- as.vector(rep(0, p))
s_hat_reparam[fix_index] <- alpha_test
# s_hat_reparam <- pinv(T_reparam) * alpha_test
F_reparam <- orth(nullspace(T_reparam))
##
##set up the test
beta_test <- mean_init
beta_test[fix_index] <- alpha_test

res <- CEL_glm_admm(beta_test, X, y, F_reparam, s_hat_reparam, family = 'poisson', beta_opt_type = 'LBFGS')

###test bootstrap
registerDoMC(parallel::detectCores())
boot_res <- CEL_glm_beta_boot(1, X, y, family = 'poisson', B = 200, 
                              outer_tol_type = 'primal_dual', dual_step = 20, 
                              wts_beta_rep = 2, outer_rel_eps = 1e-3, mirror_arg = list(mirror_eps = 1e-6, maxit = 2000))
boot_res$boot_stat
boot_res$boot_vec
#########

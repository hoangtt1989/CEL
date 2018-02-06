rm(list = ls())

set.seed(123)

pacman::p_load(tidyverse)

source('lm_admm.R')
source('lm_nested_concord.R')
source('lm_profiler.R')


n <- 100
p <- 40

# beta_true <- c(5, 4, 4, 3, 2)
beta_true <- runif(p, min = -2, max = 5)
beta_true[1] <- 5

fix_idx <- 1

# alpha_seq <- seq(1e-3, .7, length.out = 10)
alpha_add <- .5

nsim <- 50

mirror_maxit <- 2000
mirror_eps <- 1e-6
ADMM_maxit <- 2000

res_ADMM <- list()
res_nested <- list()

# mirror_converged <- NULL
# ADMM_converged <- NULL
# logelr_ADMM <- NULL
# sum_ADMM <- NULL
# time_ADMM <- NULL
# nits_ADMM <- NULL
# 
# nested_converged <- NULL
# logelr_nested <- NULL
# sum_nested <- NULL
# time_nested <- NULL
# nits_nested <- NULL
for (sim_it in 1:nsim) {
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  noise <- rnorm(n, mean = 0, sd = 1)
  # noise <- rt(n, df = 4)
  # noise <- rcauchy(n, location = 0, scale = .1) ##interesting
  # noise <- rlnorm(n, 0, 1.5) ##interesting results
  # noise <- rexp(n, rate = 1)
  
  y <- X %*% beta_true + noise
  
  
  res_ADMM[[sim_it]] <- CEL_lm_beta_fixed(beta_true[fix_idx], fix_idx, X, y, algorithm = 'ADMM', outer_maxit = ADMM_maxit, 
                                          vary_penalty = 'RB', outer_eps = 1e-8,
                                          outer_rel_eps = 1e-4, outer_tol_type = 'primal_dual', beta_opt_type = 'LBFGS', 
                                          wts_beta_rep = 2, dual_step = 8, RB_tau = 1.25, 
                                          mirror_arg = list(line_search = T, maxit = mirror_maxit, mirror_eps = mirror_eps))
  
  res_nested[[sim_it]] <- CEL_lm_beta_fixed(beta_true[fix_idx], fix_idx, X, y, algorithm = 'nested', inner_opt_type = 'LBFGS')
  
}


time_ADMM <- map_dbl(res_ADMM, 'time')
converged_ADMM <- map_lgl(res_ADMM, 'outer_converged')
mirror_ADMM <- map_dbl(res_ADMM, ~ max(.$mirror_nits))
logelr_ADMM <- map_dbl(res_ADMM, 'logelr')
primal_ADMM <- map_dbl(res_ADMM, 'primal_resid') / n

df_ADMM <- data_frame(logelr = logelr_ADMM, time = time_ADMM, primal_resid = primal_ADMM, idx = 1:nsim)

time_nested <- map_dbl(res_nested, 'time')
logelr_nested <- map_dbl(res_nested, 'logelr')
sum_wts <- map_dbl(res_nested, 'sum_wts')
primal_nested <- map_dbl(res_nested, 'primal_resid')

div_idx <- which(abs(sum_wts - 1) > 1e-5)

df_nested <- data_frame(logelr = logelr_nested, time = time_nested, primal_resid = primal_nested, idx = 1:nsim)

comb_df <- bind_rows(list(BICEP = df_ADMM, Nested = df_nested), .id = 'Algorithm')


ggplot(comb_df) +
  geom_boxplot(aes(x = '', y = time)) +
  labs(x = '', y = 'Time (s)') +
  facet_wrap(~ Algorithm, scales = 'free_y') +
  theme_bw()

ggplot(comb_df %>% filter(idx %in% div_idx, Algorithm == 'BICEP')) +
  geom_boxplot(aes(x = '', y = logelr)) +
  labs(x = '', y = 'Log EL ratio') +
  theme_bw()

rm(list = ls())

set.seed(123)

pacman::p_load(tidyverse)

source('lm_admm.R')
source('lm_nested_concord.R')
source('lm_profiler.R')


n <- 100
p <- 5

beta_true <- c(5, 4, 4, 3, 2)

# alpha_seq <- seq(1e-3, .7, length.out = 10)
alpha_add <- .5

nsim <- 100

mirror_maxit <- 3000
mirror_eps <- 1e-7
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
  # noise <- rnorm(n, mean = 0, sd = 2)
  # noise <- rt(n, df = 4)
  noise <- rcauchy(n, location = 0, scale = .1) ##interesting
  # noise <- rlnorm(n, 0, 1.5) ##interesting results
  # noise <- rexp(n, rate = 1)
  
  y <- X %*% beta_true + noise
  
  
  res_ADMM[[sim_it]] <- CEL_lm_beta_fixed(5, 1, X, y, algorithm = 'ADMM', outer_maxit = ADMM_maxit, vary_penalty = 'RB', outer_eps = 1e-8,
                                          outer_rel_eps = 1e-3, outer_tol_type = 'primal_dual',
                                        wts_beta_rep = 2, dual_step = 8, RB_tau = 1.25, mirror_arg = list(maxit = mirror_maxit, mirror_eps = mirror_eps))
  
  res_nested[[sim_it]] <- CEL_lm_beta_fixed(5, 1, X, y, algorithm = 'nested')
  
}


conv_ADMM <- map_lgl(res_ADMM, 'outer_converged')
constr_ADMM <- map_dbl(res_ADMM, 'primal_resid')
logelr_ADMM <- map_dbl(res_ADMM, 'logelr')
nits_ADMM <- map_dbl(res_ADMM, 'outer_nits')
time_ADMM <- map_dbl(res_ADMM, 'time')
mirror_nits <- map_dbl(res_ADMM, ~ max(.$mirror_nits))

conv_nested <- map_lgl(res_nested, 'outer_converged')
sum_nested <- map_dbl(res_nested, 'sum_wts')
constr_nested <- map_dbl(res_nested, 'primal_resid')
logelr_nested <- map_dbl(res_nested, 'logelr')
nits_nested <- map_dbl(res_nested, 'outer_nits')
time_nested <- map_dbl(res_nested, 'time')
inner_nits <- map_dbl(res_nested, ~ max(.$inner_nits))

div_idx <- which(abs(sum_nested - 1) > 1e-4)

median(abs(logelr_ADMM[-div_idx] - logelr_nested[-div_idx]))
median(logelr_ADMM[-div_idx] - logelr_nested[-div_idx])


div_ADMM <- data_frame(logelr = logelr_ADMM[div_idx])
ggplot(div_ADMM, aes(x = '', y = logelr)) +
  geom_boxplot(fatten = 2) +
  labs(x = '', y = 'Log EL ratio') +
  theme_bw()

ADMM_df <- data_frame(logelr = logelr_ADMM,
                      div_idx = abs(sum_nested - 1) > 1e-4,
                      type = 'BICEP')

nested_df <- data_frame(logelr = logelr_nested,
                        div_idx = abs(sum_nested - 1) > 1e-4,
                        type = 'Nested')

comb_df <- bind_rows(ADMM_df, nested_df)


ggplot(comb_df %>% filter((!div_idx & type == 'Nested') | type == 'BICEP'), aes(x = '', y = logelr)) +
  geom_boxplot(fatten = 3) +
  labs(x = '', y = 'Log EL ratio') +
  facet_wrap(~ type, scales = 'free_y') +
  theme_bw()

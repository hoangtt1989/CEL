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
  
  y <- X %*% beta_true + noise
  
  res_ADMM[[sim_it]] <- CEL_lm_beta_add(alpha_add, 1, X, y, algorithm = 'ADMM', outer_maxit = ADMM_maxit, vary_penalty = 'RB', outer_eps = 1e-7,
                                        wts_beta_rep = 4, dual_step = 8, mirror_arg = list(maxit = mirror_maxit, mirror_eps = mirror_eps))
  
  res_nested[[sim_it]] <- CEL_lm_beta_add(alpha_add, 1, X, y, algorithm = 'nested')
  
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

inner_nits[div_idx]


# trace_ADMM <- lapply(alpha_seq, function(x) {CEL_lm_beta_add(x, 1, X, y, algorithm = 'ADMM', mirror_arg = list(maxit = mirror_maxit))})
# map_lgl(trace_ADMM, 'outer_converged')
# nits_admm <- map_dbl(trace_ADMM, 'outer_nits')
# logelr_ADMM <- map_dbl(trace_ADMM, 'logelr')
# 
# trace_nested <- lapply(alpha_seq, function(x) {CEL_lm_beta_add(x, 1, X, y, algorithm = 'nested')})
# logelr_nested <- map_dbl(trace_nested, 'logelr')
# nits_nested <- map_dbl(trace_nested, 'outer_nits')
# map_lgl(trace_nested, 'outer_converged')
# 
# 
# logelr_df <- data_frame(alpha = alpha_seq,
#                         logelr_ADMM = logelr_ADMM,
#                         logelr_nested = logelr_nested) %>% 
#   mutate(idx = row_number())
# 
# plot(trace_nested[[7]]$outer_fval)

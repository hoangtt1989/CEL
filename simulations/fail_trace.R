rm(list = ls())

set.seed(123)

pacman::p_load(tidyverse)

source('lm_admm.R')
source('lm_nested_concord.R')
source('lm_profiler.R')


n <- 100
p <- 5

beta_true <- c(5, 4, 4, 3, 2)

alpha_seq <- seq(1e-3, .7, length.out = 10)


mirror_maxit <- 3000
mirror_eps <- 1e-7
ADMM_maxit <- 2000

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
noise <- rnorm(n)

y <- X %*% beta_true + noise


trace_ADMM <- lapply(alpha_seq, function(x) {CEL_lm_beta_add(x, 5, X, y, algorithm = 'ADMM', outer_maxit = ADMM_maxit, outer_eps = 1e-12, outer_rel_eps = 1e-4,
                                                             dual_step = 8, RB_tau = 1.5, 
                                                             outer_tol_type = 'primal_dual', wts_beta_rep = 4, mirror_arg = list(maxit = mirror_maxit))})
conv_ADMM <- map_lgl(trace_ADMM, 'outer_converged')
nits_ADMM <- map_dbl(trace_ADMM, 'outer_nits')
mirror_nits <- map_dbl(trace_ADMM, ~ max(.$mirror_nits))
logelr_ADMM <- map_dbl(trace_ADMM, 'logelr')
constr_ADMM <- map_dbl(trace_ADMM, 'primal_resid')

trace_nested <- lapply(alpha_seq, function(x) {CEL_lm_beta_add(x, 5, X, y, algorithm = 'nested')})
logelr_nested <- map_dbl(trace_nested, 'logelr')
nits_nested <- map_dbl(trace_nested, 'outer_nits')
conv_nested <- map_lgl(trace_nested, 'outer_converged')
sum_nested <- map_dbl(trace_nested, 'sum_wts')
inner_nits <- map_dbl(trace_nested, ~ max(.$inner_nits))
constr_nested <- map_dbl(trace_nested, 'primal_resid')


plot(trace_nested[[10]]$outer_fval)
trace_nested[[10]]$sum_wts
trace_nested[[10]]$logelr
trace_nested[[10]]$outer_tol

plot(trace_ADMM[[10]]$outer_fval)
trace_ADMM[[10]]$logelr
trace_ADMM[[10]]$outer_tol


ADMM_df <- data_frame(logelr = logelr_ADMM,
                      alpha = alpha_seq,
                      type = 'BICEP')
nested_df <- data_frame(logelr = logelr_nested,
                        alpha = alpha_seq,
                        type = 'Nested')

comb_df <- bind_rows(ADMM_df, nested_df)

ggplot(comb_df, aes(x = alpha, y = logelr)) +
  geom_point() +
  labs(y = 'Log EL ratio', x = expression(alpha)) +
  facet_wrap(~ type, scales = 'free_y') +
  theme_bw()


# tst <- CEL_lm_beta_add(alpha_seq[10], 5, X, y, algorithm = 'ADMM', outer_maxit = ADMM_maxit, outer_eps = 1e-12, outer_rel_eps = 1e-4, dual_step = 10,
#                        vary_penalty = 'RB', RB2_ksi = 2, outer_tol_type = 'primal_dual', wts_beta_rep = 4, mirror_arg = list(maxit = mirror_maxit))



alpha_seq2 <- seq(1e-2, .5, length.out = 10)

trace_ADMM2 <- lapply(alpha_seq2, function(x) {CEL_lm_beta_add(x, 5, X, y, algorithm = 'ADMM', outer_maxit = ADMM_maxit, outer_eps = 1e-12, outer_rel_eps = 1e-4,
                                                              dual_step = 10,
                                                              outer_tol_type = 'primal_dual', wts_beta_rep = 4, mirror_arg = list(maxit = mirror_maxit))})

conv_ADMM2 <- map_lgl(trace_ADMM2, 'outer_converged')
nits_ADMM2 <- map_dbl(trace_ADMM2, 'outer_nits')
mirror_nits2 <- map_dbl(trace_ADMM2, ~ max(.$mirror_nits))
logelr_ADMM2 <- map_dbl(trace_ADMM2, 'logelr')
constr_ADMM2 <- map_dbl(trace_ADMM2, 'primal_resid')


trace_nested2 <- lapply(alpha_seq2, function(x) {CEL_lm_beta_add(x, 5, X, y, algorithm = 'nested')})
logelr_nested2 <- map_dbl(trace_nested2, 'logelr')
nits_nested2 <- map_dbl(trace_nested2, 'outer_nits')
conv_nested2 <- map_lgl(trace_nested2, 'outer_converged')
sum_nested2 <- map_dbl(trace_nested2, 'sum_wts')
inner_nits2 <- map_dbl(trace_nested2, ~ max(.$inner_nits))
constr_nested2 <- map_dbl(trace_nested2, 'primal_resid')


tst <- CEL_lm_beta_add(alpha_seq2[5], 5, X, y, algorithm = 'ADMM', outer_maxit = ADMM_maxit, outer_rel_eps = 1e-4,
                       dual_step = 2,
                       outer_tol_type = 'primal_dual', wts_beta_rep = 4, mirror_arg = list(maxit = mirror_maxit))

rm(list = ls())

set.seed(123)

pacman::p_load(tidyverse, ggpubr)

source('lm_admm.R')
source('lm_nested_concord.R')
source('lm_profiler.R')


n <- 100
p <- 5

beta_true <- c(5, 4, 4, 3, 2)

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
noise <- rnorm(n)

y <- X %*% beta_true + noise

mean_init <- solve(t(X) %*% X, t(X) %*% y)

nsim <- 100

mirror_maxit <- 2000
mirror_eps <- 1e-6
ADMM_maxit <- 2000
ADMM_eps <- 1e-8

res_ADMM <- list()
res_nested <- list()


for (sim_it in 1:nsim) {
  
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
  # beta_test <- mean_init + rnorm(p)
  # beta_test <- mean_init + runif(p, min = -.2, max = .2) ##decrease these numbers and nested behaves better
  # beta_test <- mean_init + runif(p, min = -.3, max = .3)
  # beta_test <- mean_init + runif(p, min = -.4, max = .4)
  beta_test <- mean_init + runif(p, min = -.5, max = .5)
  beta_test[fix_index] <- alpha_test
  
  wts_init <- runif(n, min = 1, max = 10)
  wts_init <- wts_init / sum(wts_init) * n
  
  res_ADMM[[sim_it]] <- CEL_lm_admm(beta_test, X, y, F_reparam, s_hat_reparam, wts_init = NULL,
                                    mirror_arg = list(maxit = mirror_maxit, mirror_eps = mirror_eps, line_search = T), outer_rel_eps = 1e-4,
                                    wts_beta_rep = 4, RB_tau = 1.5, outer_tol_type = 'primal_dual',
                                    outer_eps = ADMM_eps, outer_maxit = ADMM_maxit, dual_step = 8)
  
  res_nested[[sim_it]] <- CEL_lm_nested(beta_test, X, y, F_reparam, s_hat_reparam, gamma_init = 'random')
  
}

###soln without random start
soln_ADMM <- CEL_lm_beta_add(.1, 1, X, y, algorithm = 'ADMM', mirror_arg = list(maxit = mirror_maxit, mirror_eps = mirror_eps, line_search = T), outer_rel_eps = 1e-4,
                             wts_beta_rep = 4, RB_tau = 1.5, outer_tol_type = 'primal_dual',
                             outer_eps = ADMM_eps, outer_maxit = ADMM_maxit, dual_step = 8)
soln_nested <- CEL_lm_beta_add(.1, 1, X, y, algorithm = 'nested')
###

mean(abs(soln_ADMM$mean_init - soln_ADMM$beta_opt))

time_ADMM <- map_dbl(res_ADMM, 'time')
nits_ADMM <- map_dbl(res_ADMM, 'outer_nits')
mirror_nits <- map_dbl(res_ADMM, ~ max(.$mirror_nits))
logelr_ADMM <- map_dbl(res_ADMM, 'logelr')
beta_ADMM <- map_dbl(res_ADMM, ~ median(.$beta_opt))
wts_ADMM <- map_dbl(res_ADMM, ~ min(.$wts))
constr_ADMM <- map_dbl(res_ADMM, ~ median(.$primal_resid))

ADMM_df <- data_frame(type = 'BICEP',
                      logelr_soln = soln_ADMM$logelr,
                      logelr = logelr_ADMM,
                      beta = beta_ADMM,
                      beta_soln = median(soln_ADMM$beta_opt),
                      wts_soln = min(soln_ADMM$wts),
                      wts = wts_ADMM,
                      constr = constr_ADMM / n,
                      constr_soln = soln_ADMM$primal_resid /n)

# ggplot(ADMM_df, aes(y = logelr, x = '')) +
#   geom_boxplot() +
#   geom_hline(yintercept = soln_ADMM$logelr, linetype = 2)

time_nested <- map_dbl(res_nested, 'time')
logelr_nested <- map_dbl(res_nested, 'logelr')
beta_nested <- map_dbl(res_nested, ~ median(.$beta_opt))
wts_nested <- map_dbl(res_nested, ~ min(.$wts))
constr_nested <- map_dbl(res_nested, ~ median(.$primal_resid))
sum_nested <- map_dbl(res_nested, 'sum_wts')

sum1_idx <- abs(sum_nested - 1) < 1e-6
# logelr_nested[sum1_idx]

# boxplot(logelr_nested[sum1_idx])
# boxplot(constr_nested[sum1_idx])

nested_df <- data_frame(type = 'Nested',
                      logelr_soln = soln_nested$logelr,
                      logelr = logelr_nested,
                      beta = beta_nested,
                      beta_soln = median(soln_nested$beta_opt),
                      wts = wts_nested,
                      wts_soln = min(soln_nested$wts),
                      constr = constr_nested,
                      constr_soln = soln_nested$primal_resid) %>% 
  mutate(sum1_idx = sum1_idx)

ADMM_df2 <- ADMM_df %>% 
  mutate(sum1_idx = sum1_idx)

# ggplot(nested_df, aes(y = logelr, x = '')) +
  # geom_boxplot()


comb_df <- bind_rows(ADMM_df2, nested_df)

logelr_pt <- ggplot(comb_df %>% filter((sum1_idx & type == 'Nested') | type == 'BICEP'), aes(y = logelr, x = '')) +
  geom_boxplot(fatten = 3) +
  # geom_violin() +
  geom_hline(aes(yintercept = logelr_soln), linetype = 4) +
  facet_wrap(~ type, scales = 'free_y', shrink = F) +
  labs(x = '', y = 'Log EL ratio') +
  theme_bw()

beta_pt <- ggplot(comb_df%>% filter((sum1_idx & type == 'Nested') | type == 'BICEP'), aes(y = abs(beta), x = '')) +
  geom_boxplot(fatten = 3) +
  geom_hline(aes(yintercept = beta_soln), linetype = 4) +
  facet_wrap(~ type, scales = 'free_y') +
  labs(x = '', y = expression('Median of'~paste('|', beta, '|', sep = ''))) +
  theme_bw()

wts_pt <- ggplot(comb_df%>% filter((sum1_idx & type == 'Nested') | type == 'BICEP'), aes(y = wts, x = '')) +
  geom_boxplot(fatten = 3) +
  geom_hline(aes(yintercept = wts_soln), linetype = 4) +
  facet_wrap(~ type, scales = 'free_y') +
  labs(x = '', y = 'Minimum of weights') +
  theme_bw()

primal_pt <- ggplot(comb_df%>% filter((sum1_idx & type == 'Nested') | type == 'BICEP'), aes(y = constr, x = '')) +
  geom_boxplot(fatten = 3) +
  # geom_violin() +
  geom_hline(aes(yintercept = constr_soln), linetype = 4) +
  facet_wrap(~ type, scales = 'free_y') +
  labs(x = '', y = 'Primal residual') +
  theme_bw()

# gath_df <- comb_df %>% 
#   gather(key, value, beta, wts)
# 
# ggplot(gath_df%>% filter((sum1_idx & type == 'Nested') | type == 'ADMM'), aes(y = value, x = '')) +
#   geom_boxplot(fatten = 3) +
#   facet_grid(key~ type, scales = 'free')

ggarrange(logelr_pt, primal_pt, ncol = 2, labels = c('A', 'B'))

ggarrange(beta_pt, wts_pt, ncol = 2, labels = c('A', 'B'))

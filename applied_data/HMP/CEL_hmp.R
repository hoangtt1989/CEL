rm(list = ls())

pacman::p_load(tidyverse, MASS, doMC, phyloseq, stringr)

source('mirror_descent.R')
source('glm_admm.R')
source('glm_profiler.R')
source('owenEL.R')
# source('lm_profiler.R')
# source('lm_nested_concord.R')
# source('owenEL_old.R')
# source('damped_newton.R')
# source('lbfgs_step.R')
# source('grad_step.R')


OTU <- readRDS('applied_data/HMP/OTU_top10.rds')
sex <- readRDS('applied_data/HMP/sample_data_sex.rds')
site <- readRDS('applied_data/HMP/sample_data_site.rds')

predictors <- model.matrix(~ site + sex)

###set up response and predictors
confint_trace <- function(idx, coef_idx = 2) {
  resp <- OTU[, idx]
  boot_res <- CEL_glm_beta_boot(coef_idx, predictors, resp, parallel = T, 
                                family = 'poisson', dual_step = 20, wts_beta_rep = 2, RB_tau = 1.5, 
                                outer_rel_eps = 1e-3, outer_tol_type = 'primal_dual', mirror_arg = list(mirror_eps = 1e-6, maxit = 2000))
  boot_confint <- CEL_glm_beta_profiler(coef_idx, predictors, resp, family = 'poisson', 
                                        test_thresh = boot_res$boot_stat, dual_step = 20, wts_beta_rep = 2, 
                                        RB_tau = 1.5, outer_rel_eps = 1e-3, outer_tol_type = 'primal_dual', mirror_arg = list(mirror_eps = 1e-6, maxit = 2000))
  chisq_confint <- CEL_glm_beta_profiler(coef_idx, predictors, resp, family = 'poisson',
                                         test_thresh = 'chisq', dual_step = 20, wts_beta_rep = 2, 
                                         RB_tau = 1.5, outer_rel_eps = 1e-3, outer_tol_type = 'primal_dual', mirror_arg = list(mirror_eps = 1e-6, maxit = 2000))
  glm_mod <- glm(resp ~ predictors - 1, family = poisson(link = 'log'))
  glm_confint <- confint(glm_mod, coef_idx)
  confint_df <- data_frame(OTU = colnames(OTU)[idx], lower_boot = boot_confint$lower_int, upper_boot = boot_confint$upper_int,
                           lower_chisq = chisq_confint$lower_int, upper_chisq = chisq_confint$upper_int,
                           lower_glm = glm_confint[1], upper_glm = glm_confint[2], mean = glm_mod$coefficients[coef_idx])
  return(confint_df)
}
registerDoMC(parallel::detectCores())
confint_df <- map_df(1:ncol(OTU), confint_trace)

confint_df2 <- confint_df %>% 
  gather(key, value, -OTU, -mean) %>% 
  mutate(confint_bound = str_split_fixed(key, '_', 2)[, 1],
         confint_type = str_split_fixed(key, '_', 2)[, 2])

confint_df3 <- confint_df2 %>% 
  group_by(confint_type) %>% 
  dplyr::select(-key) %>% 
  spread(confint_bound, value) %>% 
  ungroup() %>% 
  mutate(confint_type = factor(confint_type, labels = c('BS-CEL', 'CS-CEL', 'Prof')))

ggplot(confint_df3, aes(x = confint_type, y = mean)) +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = '', y = '') +
  facet_wrap(~ OTU) +
  theme_bw()


######repeat for sex
registerDoMC(parallel::detectCores())
confint_sex_df <- map_df(1:ncol(OTU), confint_trace, coef_idx = 3)

confint_sex_df2 <- confint_sex_df %>% 
  gather(key, value, -OTU, -mean) %>% 
  mutate(confint_bound = str_split_fixed(key, '_', 2)[, 1],
         confint_type = str_split_fixed(key, '_', 2)[, 2])

confint_sex_df3 <- confint_sex_df2 %>% 
  group_by(confint_type) %>% 
  dplyr::select(-key) %>% 
  spread(confint_bound, value) %>% 
  ungroup() %>% 
  mutate(confint_type = factor(confint_type, labels = c('BS-CEL', 'CS-CEL', 'Prof')))

ggplot(confint_sex_df3, aes(x = confint_type, y = mean)) +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = '', y = '') +
  facet_wrap(~ OTU) +
  theme_bw()
######


##some interval tests
resp <- OTU[, 'OTU_97.32083']
beta_test <- glm.fit(predictors, resp, family = poisson(link = 'log'), intercept = F)$coefficients
ineq_test <- CEL_glm_admm(beta_test, predictors, resp, family = 'poisson', prox_fun = prox_inequality, prox_fun_arg = list(test_idx = 2:3, test_val = 1), 
                         outer_tol_type = 'primal_dual', wts_beta_rep = 2, dual_step = 20, RB_tau = 1.5, mirror_arg = list(line_search = T, mirror_eps = 1e-6, maxit = 2000))
ineq_boot <- CEL_glm_boot(predictors, resp, B = 200, family = 'poisson', prox_fun = prox_inequality, prox_fun_arg = list(test_idx = 2:3, test_val = 1), 
                          outer_tol_type = 'primal_dual', wts_beta_rep = 2, dual_step = 20, RB_tau = 1.5, mirror_arg = list(line_search = T, mirror_eps = 1e-6, maxit = 2000))
boot2 <- CEL_glm_beta_boot(2:3, predictors, resp, family = 'poisson', outer_tol_type = 'primal_dual', wts_beta_rep = 2, dual_step = 20, RB_tau = 1.5, mirror_arg = list(line_search = T, mirror_eps = 1e-6, maxit = 2000))
length(which(boot2$boot_vec >= ineq_test$logelr_stat)) / length(boot2$boot_vec)
ineq_boot$boot_stat
ineq_boot$boot_pval
##

####debugging
# idx <- 3
# resp <- OTU[, idx]
# boot_res <- CEL_glm_beta_boot(2, predictors, resp, parallel = T, 
#                               family = 'poisson', dual_step = 20, wts_beta_rep = 2, RB_tau = 1.5, 
#                               outer_rel_eps = 1e-3, outer_tol_type = 'primal_dual', mirror_arg = list(mirror_eps = 1e-6, maxit = 2000))
# boot_confint <- CEL_glm_beta_profiler(2, predictors, resp, family = 'poisson', 
#                                       test_thresh = boot_res$boot_stat, dual_step = 20, wts_beta_rep = 2, 
#                                       RB_tau = 1.5, outer_rel_eps = 1e-3, outer_tol_type = 'primal_dual', mirror_arg = list(mirror_eps = 1e-6, maxit = 2000))
# chisq_confint <- CEL_glm_beta_profiler(2, predictors, resp, family = 'poisson',
#                                        test_thresh = 'chisq', dual_step = 20, wts_beta_rep = 2, 
#                                        RB_tau = 1.5, outer_rel_eps = 1e-3, outer_tol_type = 'primal_dual', mirror_arg = list(mirror_eps = 1e-6, maxit = 2000))
#####



# #######contours
# g <- expand.grid(site = seq(-.7, 1, by = .1), theta=seq(-.7, .7, by = .1))
# # elm_fixed_fun <- function(X, y, beta, family, ...) {
# #   link_fun <- switch(family,
# #                      gaussian = base::identity,
# #                      binomial = binomial_link,
# #                      poisson = exp)
# #   R <- as.vector(y - link_fun(X %*% beta))
# #   Z <- R * X
# #   elm_res <- emplik(Z, ...)
# #   return(-2 * elm_res$logelr)
# # }
# ll_grid <- apply(g, 1, function(inp){elm_fixed_fun(predictors, resp, c(glm_mod$coefficients[1], inp[1], inp[2]), family = 'poisson')})
# df <- data_frame(site = g[, 1], sex = g[, 2], z = ll_grid)
# ggplot(df) + geom_contour(aes(x = site, y = sex, z = z), breaks = c(qchisq(.95, 1))) + 
#   geom_hline(yintercept = glm_mod$coefficients[3], linetype = 2) +
#   geom_vline(xintercept = glm_mod$coefficients[2], linetype = 2)
# ##CEL
# CEL_fixed_fun <- function(X, y, beta_fix_val, fix_index, ...) {
#   res <- CEL_glm_beta_fixed(beta_fix_val, fix_index, X, y, ...)
#   return(res$logelr_stat)
# }
# CEL_grid <- apply(g, 1, function(inp){CEL_fixed_fun(predictors, resp, c(inp[1], inp[2]), c(2, 3), family = 'poisson', dual_step = 20, 
#                                                     wts_beta_rep = 2, RB_tau = 1.5, outer_rel_eps = 1e-3, 
#                                                     outer_tol_type = 'primal_dual', mirror_arg = list(mirror_eps = 1e-5, maxit = 2000))})
# CEL_df <- data_frame(site = g[, 1], sex = g[, 2], z = CEL_grid)
# ggplot(CEL_df) + geom_contour(aes(x = site, y = sex, z = z), breaks = c(qchisq(.95, 1))) + 
#   geom_hline(yintercept = glm_mod$coefficients[3], linetype = 2) +
#   geom_vline(xintercept = glm_mod$coefficients[2], linetype = 2)
# ##
# #######

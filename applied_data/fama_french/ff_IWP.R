rm(list = ls())

set.seed(123)

pacman::p_load(tidyverse, stringr, readr, janitor, doMC)

source('lm_admm.R')
source('lm_nested_concord.R')
source('lm_profiler.R')


ff_df <- read_csv('applied_data/fama_french/US_5_Factors_Daily.csv', skip = 3) %>% 
  clean_names() %>% 
  rename(date = x1)

ff_df2 <- ff_df %>%
  mutate(date = as.character(date)) %>% 
  filter(str_detect(date, '2016') | str_detect(date, '2017')) %>% 
  mutate(date = as.Date(date, '%Y%m%d'))


fund_df <- read_csv('applied_data/fama_french/IWP.csv') %>% 
  clean_names()

fund_df

fund_pr <- fund_df %>% 
  pull(adj_close)

fund_len <- length(fund_pr)

fund_rr <- (fund_pr[2:fund_len] - fund_pr[1:(fund_len - 1)]) / fund_pr[1:(fund_len - 1)] * 100
# fund_rr <- (fund_pr[2:fund_len] - fund_pr[1:(fund_len - 1)]) ##return, not rate of return

fund_df2 <- fund_df %>%
  mutate(rr = c(NA, fund_rr)) %>% 
  mutate(date = as.Date(date, '%m/%d/%y')) %>% 
  dplyr::select(date, fund_rr = rr) %>% 
  filter(!is.na(fund_rr))

comb_df <- inner_join(ff_df2, fund_df2, by = 'date') %>% 
  mutate(fund_rr = fund_rr - rf)

mod_lm <- lm(fund_rr ~ mkt_rf + smb + hml + rmw + cma, data = comb_df)
summary(mod_lm)

mat_form <- comb_df %>%
  mutate(intercept = 1) %>% 
  dplyr::select(intercept, everything(), -date, -rf)

mat_form

design_mat <- mat_form %>% 
  dplyr::select(-fund_rr) %>% 
  as.matrix()

y <- mat_form %>% 
  pull(fund_rr)

##heavy tails?
qqnorm(mod_lm$residuals)
qqline(mod_lm$residuals)

qqnorm(y)
qqline(y)
source('qqplot_ggplot.R')
qqplot.data(mod_lm$residuals) +
  theme_bw()
##

####EL testing

##for each coefficient - bootstrap .95 statistic, CI using chi2, CI using bootstrap

##boot strap .95 statistic
B <- 200
registerDoMC(parallel::detectCores())
boot_trace_ADMM <- mclapply(1:ncol(design_mat), CEL_lm_beta_boot, design_mat, y, B = B, parallel = F,
                            algorithm = 'ADMM', dual_step = 50, RB_tau = 1.5, wts_beta_rep = 3,
                            outer_tol_type = 'primal_dual', outer_rel_eps = 1e-3, 
                            mirror_arg = list(maxit = 2000, mirror_eps = 1e-7))
boot_stat_ADMM <- map_dbl(boot_trace_ADMM, 'boot_stat')
boot_trace_nested <- mclapply(1:ncol(design_mat), CEL_lm_beta_boot, design_mat, y, B = B, parallel = F,
                              algorithm = 'nested')
boot_stat_nested <- map_dbl(boot_trace_nested, 'boot_stat')
##

##CI using chi2
CI_trace_chi2_ADMM <- mclapply(1:ncol(design_mat), CEL_lm_beta_profiler, design_mat, y,
                               algorithm = 'ADMM', dual_step = 50, RB_tau = 1.5, wts_beta_rep = 3,
                               outer_tol_type = 'primal_dual', outer_rel_eps = 1e-3,
                               mirror_arg = list(maxit = 2000, mirror_eps = 1e-7))
CI_trace_chi2_nested <- mclapply(1:ncol(design_mat), CEL_lm_beta_profiler, design_mat, y,
                                 algorithm = 'nested')

CI_chi2_ADMM <- map_df(CI_trace_chi2_ADMM, function (x) {upper <- x$upper_int
lower <- x$lower_int
mean_estimate <- x$mean_val
data_frame(upper = upper, lower = lower, mean_estimate = mean_estimate)})
CI_chi2_nested <- map_df(CI_trace_chi2_nested, function (x) {upper <- x$upper_int
lower <- x$lower_int
mean_estimate <- x$mean_val
data_frame(upper = upper, lower = lower, mean_estimate = mean_estimate)})
##

##CI using boot stat
CI_trace_boot_ADMM <- mcmapply(function(idx, stat) {CEL_lm_beta_profiler(idx, design_mat, y, test_thresh = stat,
                                                                         algorithm = 'ADMM', dual_step = 50, RB_tau = 1.5, wts_beta_rep = 3,
                                                                         outer_tol_type = 'primal_dual', outer_rel_eps = 1e-3,
                                                                         mirror_arg = list(maxit = 2000, mirror_eps = 1e-7))},
                               1:ncol(design_mat), boot_stat_ADMM, SIMPLIFY = F)
CI_trace_boot_nested <- mcmapply(function(idx, stat) {CEL_lm_beta_profiler(idx, design_mat, y, test_thresh = stat,
                                                                           algorithm = 'nested')},
                                 1:ncol(design_mat), boot_stat_nested, SIMPLIFY = F)
CI_boot_ADMM <- map_df(CI_trace_boot_ADMM, function (x) {upper <- x$upper_int
lower <- x$lower_int
mean_estimate <- x$mean_val
data_frame(upper = upper, lower = lower, mean_estimate = mean_estimate)})
CI_boot_nested <- map_df(CI_trace_boot_nested, function (x) {upper <- x$upper_int
lower <- x$lower_int
mean_estimate <- x$mean_val
data_frame(upper = upper, lower = lower, mean_estimate = mean_estimate)})
##

##normal theory
CI_normal <- confint(mod_lm) %>% 
  as.data.frame() %>% 
  rownames_to_column('coef') %>% 
  as_tibble() %>% 
  clean_names() %>% 
  rename(upper = x97_5_percent, lower = x2_5_percent) %>% 
  mutate(mean_estimate = CI_boot_ADMM$mean_estimate)
##

##combine DFs
CI_ELs <- map(list(CI_chisq_ADMM = CI_chi2_ADMM,
                   CI_chisq_nested = CI_chi2_nested,
                   CI_boot_ADMM = CI_boot_ADMM,
                   CI_boot_nested = CI_boot_nested), ~ mutate(., coef = CI_normal$coef))
CI_comb <- bind_rows(list(CI_normal = CI_normal), CI_ELs, .id = 'type') %>% 
  filter(! (type %in% c('CI_chisq_nested', 'CI_boot_nested'))) %>% 
  mutate(type = factor(type, labels = c('BS-CEL', 'CS-CEL', 'Normal'))) %>% 
  mutate(coef = factor(coef, levels = c('(Intercept)', 'mkt_rf', 'smb', 'hml', 'rmw', 'cma'), labels = c('Intercept', 'MKT', 'SMB', 'HML', 'RMW', 'CMA')))
##

##plot
ggplot(CI_comb, aes(x = type, y = mean_estimate, group = 1)) +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = '', y = '') +
  facet_wrap(~ coef, scales = 'free_y') +
  theme_bw()
##

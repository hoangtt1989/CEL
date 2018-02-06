rm(list = ls())

pacman::p_load(tidyverse, phyloseq)

source('applied_data/HMP/phylo_utility.R')

phy <- readRDS('applied_data/HMP/HMPv35.rds')

table(sample_data(phy)$HMPbodysubsite)

X_var <- 'Tongue_dorsum'
Y_var <- 'Throat'

phy_site <- prune_samples(sample_data(phy)$HMPbodysubsite %in% c(X_var, Y_var) & sample_data(phy)$visitno %in% c(1, 2), phy)

site_tbl <- sample_data(phy_site) %>% 
  as_tibble()

##subset to subjects with both visits in both sites
site_both <- site_tbl %>%
  nest(-HMPbodysubsite) %>%
  mutate(site_both = map(data, ~ group_by(., RSID) %>%
                           summarise(num_visit1 = length(which(visitno == 1)),
                                     num_visit2 = length(which(visitno == 2))) %>% 
                           filter(num_visit1 == 1 & num_visit2 == 1)))
site_both_tbl <- site_tbl %>% 
  filter(RSID %in% intersect(site_both$site_both[[1]]$RSID, site_both$site_both[[2]]$RSID))

table(site_both_tbl$HMPbodysubsite, site_both_tbl$visitno)

## sample id for Y
sam_Y <- site_both_tbl %>%
  filter(HMPbodysubsite == Y_var) %>% 
  nest(-visitno) %>% 
  mutate(data = map(data, ~ arrange(., RSID)))
identical(sam_Y$data[[1]]$RSID, sam_Y$data[[2]]$RSID)
  
## sample id for X
sam_X <- site_both_tbl %>%
  filter(HMPbodysubsite == X_var) %>% 
  nest(-visitno) %>% 
  mutate(data = map(data, ~ arrange(., RSID)))
identical(sam_X$data[[1]]$RSID, sam_X$data[[2]]$RSID)

identical(sam_Y$data[[1]]$RSID, sam_X$data[[1]]$RSID)
identical(sam_X$data[[1]]$RSID, sam_X$data[[1]]$RSID)

sam_Y2 <- sam_Y %>% 
  unnest()
sam_X2 <- sam_X %>%
  unnest()

phy_site_both <- prune_samples(sample_data(phy_site)$RSID %in% site_both_tbl$RSID, phy_site)

##filter for prevalence
phy_site_filt <- phyloseq_filter(phy_site_both, .25)
phy_site_filt

##for differential abundance - extract the 10 most abundant taxa from time 1
phy_time1 <- prune_samples(sample_data(phy_site_filt)$visitno == 1, phy_site_filt)
OTU_time1 <- phyloseq_OTU(phy_time1)
OTU_top10_ord <- order(colSums(OTU_time1))
OTU_top10 <- OTU_time1[, OTU_top10_ord[1:10]]
##

##extract response and sex
sex <- sample_data(phy_time1)$sex
site <- sample_data(phy_time1)$HMPbodysubsite
##

saveRDS(OTU_top10, 'OTU_top10.rds')
saveRDS(sex, 'sample_data_sex.rds')
saveRDS(site, 'sample_data_site.rds')

# 
# 
# ##extract OTU table and process
# OTU_X_all <- phyloseq_OTU(prune_samples(sample_data(phy_site_filt)$HMPbodysubsite == X_var, phy_site_filt))
# OTU_Y_all <- phyloseq_OTU(prune_samples(sample_data(phy_site_filt)$HMPbodysubsite == Y_var, phy_site_filt))
# 
# OTU_X_all_asinh <- asinh(OTU_X_all)
# OTU_X_all_std <- scale(OTU_X_all_asinh, scale = T, center = T)
# 
# OTU_Y_all_asinh <- asinh(OTU_Y_all)
# OTU_Y_all_std <- scale(OTU_Y_all_asinh, scale = T, center = T)
# 
# ##do PCA on Y
# svd_Y <- svd(OTU_Y_all_std)
# svd_Y_top50 <- svd_Y$u[, 1:50] %*% diag(svd_Y$d[1:50])
# rownames(svd_Y_top50) <- rownames(OTU_Y_all)
# 
# 
# ###train set
# Y_train <- svd_Y_top50[sam_Y2 %>% 
#                          filter(visitno == 1) %>% 
#                          pull(X.SampleID) %>% 
#                          as.character(), ]
# X_train <- OTU_X_all_std[sam_X2 %>% 
#                            filter(visitno == 1) %>% 
#                            pull(X.SampleID) %>% 
#                            as.character(), ]
# 
# 
# ###test set
# Y_test <- svd_Y_top50[sam_Y2 %>% 
#                          filter(visitno == 2) %>% 
#                          pull(X.SampleID) %>% 
#                          as.character(), ]
# X_test <- OTU_X_all_std[sam_X2 %>% 
#                            filter(visitno == 2) %>% 
#                            pull(X.SampleID) %>% 
#                            as.character(), ]
# 
# write_csv(as.data.frame(Y_train), 'Y_top50_train.csv', col_names = F)
# write_csv(as.data.frame(X_train), 'X_train.csv', col_names = F)
# write_csv(as.data.frame(Y_test), 'Y_top50_test.csv', col_names = F)
# write_csv(as.data.frame(X_test), 'X_test.csv', col_names = F)

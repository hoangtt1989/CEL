# pacman::p_load(dplyr, tibble, purrr, phyloseq)
#extract OTU table from phyloseq
phyloseq_OTU <- function(phy, df = F){
  #set df = T to convert to data.frame
  #automatic conversion to samples X taxa
  OTU <- as(otu_table(phy), 'matrix')
  if(taxa_are_rows(phy)){
    OTU <- t(OTU)
  }
  if(df){
    OTU <- as.data.frame(OTU)
  }
  return(OTU)
}

#extract sample data from phyloseq
phyloseq_sam <- function(phy){
  sam_df <- as(sample_data(phy), 'data.frame')
  return(sam_df)
}

#clr function with special handling of zeros
clr_func <- function(x.f, base=exp(1), tol=.Machine$double.eps) {
  nzero <- (x.f >= tol)
  LOG <- log(ifelse(nzero, x.f, 1), base)
  ifelse(nzero, LOG - mean(LOG)/mean(nzero), 0.0)
}

#clr transform
phyloseq_clr <- function(phy_object, special_zero = T){
  if(!special_zero){
    otu_table(phy_object) <- otu_table(phy_object) + .001
  }
  otu_table(phy_object) <- transform_sample_counts(otu_table(phy_object), clr_func)
  return(phy_object)
}

#filter for prevalence
phyloseq_filter <- function(phy_object, filter_pct = 0.25){
  #this function filters OTUs for prevalence i.e. OTUs must be prevalent in at least "filter_pct %" of samples to be kept
  #phy_object: a phyloseq object
  #filter_pct: the prevalence percentage in decimal form
  if(filter_pct > 1){
    stop('You entered a prevalence threshold greater than 1. This number must be in decimal form.')
  }
  keep <- genefilter_sample(phy_object, filterfun_sample(function(x) x > 0), A = filter_pct * nsamples(phy_object))
  phy_new <- prune_taxa(keep, phy_object)
  return(phy_new)
}

#condense an OTU table to genus/family/etc.
#this is actually redundant because
#I just discovered phyloseq has a function called tax_glom that seems to do the same thing
taxa_comb <- function(phy, lvl_in = 'Genus'){
  OTU <- phyloseq_OTU(phy) %>% t %>%
    as.data.frame %>%
    rownames_to_column
  
  #get tax_table
  tax_tab <- tax_table(phy, errorIfNULL = T)
  tax_levels <- colnames(tax_tab)
  #post error if input is not included in tax_table
  if(!(lvl_in %in% tax_levels)){
    err_msg <- paste('Must supply a valid taxa classification. Possible values are', 
                     paste(tax_levels, collapse = ' '), sep = ' ')
    stop(err_msg)
  }
  lvl_ind <- which(tax_levels == lvl_in)
  #getting the "unknown" value for this level, will be something like g__
  lvl_unk <- paste(tolower(strsplit(lvl_in, '')[[1]][1]), '__', sep = '')
  #find "unknowns"
  tax_unk <- which(tax_tab[,lvl_ind] == lvl_unk)
  #replace the unknowns with the next "lowest" taxa level
  if(lvl_ind > 1){
    tax_tab[tax_unk, lvl_ind] <- tax_tab[tax_unk, lvl_ind - 1]
  }else{
    warning('You have reached the lowest taxa classification')
  }
  
  #now create an OTU table at the genus/family/etc. level
  if(!identical(OTU$rowname, rownames(tax_tab))){
    #not necessary if OTU comes from the same phyloseq object
    stop('OTUs are not identical in OTU table and tax_table')
  }
  OTU$rowname <- as.factor(tax_tab[,lvl_ind])
  lvl_freq <- table(OTU$rowname)
  OTU_calc <- OTU %>% slice_rows('rowname') %>% by_slice(colSums, .collate = c('cols'))
  OTU_sub <- OTU_calc[,-1] %>% t
  rownames(OTU_sub) <- sample_names(phy)
  colnames(OTU_sub) <- OTU_calc$rowname
  OTU_lvl <- as.data.frame(OTU_sub)
  #OTU is a new OTU table with OTUs replaced by genus/family, etc.
  #unk_ind is the indices of OTUs that were unclassified at the input level
  #lvl_freq is a frequency table of the OTU classifications at the input level
  return(list(OTU = OTU_lvl, unk_ind = tax_unk, lvl_freq = lvl_freq))
}
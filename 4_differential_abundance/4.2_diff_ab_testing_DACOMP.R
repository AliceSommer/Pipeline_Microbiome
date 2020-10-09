library(dacomp)
library(phyloseq); packageVersion("phyloseq")

# For permutation load W matrix for randomization test
# load("/Users/alicesommer/Desktop/Bureau/DOCTORATE/data_pipeline_microbiome/W_paired_smoke_bis.Rdata")
load("/Users/alicesommer/Desktop/Bureau/DOCTORATE/data_pipeline_microbiome/W_paired_PM25.Rdata")
# W_paired <- W_paired_smoke
nrep <- 10^4

# reference sets 
# load('/Users/alicesommer/Desktop/DACOMP_cluster/dacomp_results/selected_ref.RData')
load('/Users/alicesommer/Desktop/DACOMP_cluster/dacomp_results/selected_ref_PM.RData')

# rank_names
rank_names <- levels(selected_ref$rank)

# save results
dacomp_results <- NULL

verbose = T

for(k in 1:length(rank_names)){
  
  # load data ps
  ps_work <- readRDS(paste0("/Users/alicesommer/Desktop/DACOMP_cluster/ps_filt_",k,"_PM.rds"))
  
  Y = sample_data(ps_work)$W # research variable
  # Permutated "W" matrix
  Y_perm <- cbind(as.numeric(Y == 1),W_paired[,1:nrep]) # the first column of the Y_perm matrix is W_obs
  
  # Microbiome Data
  X = as(otu_table(ps_work), "matrix")
  
  # reference set
  ind_reference_taxa = selected_ref$condition[selected_ref$rank == rank_names[k]]
  
  ###############################################################################################
  #### TEST WITH DACOMP-RATIO METHOD ####
  
  p_test = ncol(X)
  n_test = nrow(X)
  stats_matrix_test = matrix(NA, ncol = p_test, nrow = nrep+1)
  ratio_matrix_test = matrix(NA, nrow = n_test, ncol = 1)
  ratio_matrix_test_raw = matrix(NA, nrow = p_test, ncol = 1)
  
  # Compute reference values
  reference_values = apply(X[,ind_reference_taxa], 1, sum)
  
  print(paste(rank_names[k],'ref. sanity check:', sum(reference_values == 0)))
  
  #iterate over taxa and test
  for(t in 1:p_test){
    
    if(verbose)
      if(t%% ceiling(p_test/10) == 1)
        cat(paste0('Testing taxon : ',t,'/',p_test,' \n\r'))
    
    nom_test = X[,t]
    dnom_test = reference_values
    
    #no need to test reference taxa
    if(t %in% ind_reference_taxa){
      print(t)
      # nom_test = apply(X_test[,ind_reference_taxa[-which(ind_reference_taxa == t)]], 1, sum)
      next
    }
    
    ratio_matrix_test[,1] = nom_test/(dnom_test+nom_test)
    
    ratio_matrix_test_raw[t,] = mean(log10(ratio_matrix_test[Y_perm[,1]==0,1]+1)) - mean(log10(ratio_matrix_test[Y_perm[,1]==1,1]+1))
    
    stats_matrix_test[,t] = dacomp:::Compute.resample.test(ratio_matrix_test,Y_perm, 
                                                           statistic = DACOMP.TEST.NAME.LOG_FOLD_DIFFERENCE_IN_MEANS)
    
  }
  
  # computes p-values:
  p.values.ratio.normalization = apply(stats_matrix_test, 2, function(x) mean(x >= x[1]))
  
  # head(sort(p.values.ratio.normalization),20)
  # which(p.adjust(p.values.ratio.normalization,method = 'BH') <= 0.15)
  # unname(tax_table(ps_work)[which(p.adjust(p.values.ratio.normalization,method = 'BH')<=0.2),])
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Multiple comparison adjustment 
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  ### Lee, Forastiere, Miratix, Pillai (2017) method  for multiple comparison adjustment ### 
  
  # STEP 1 to 3: recorded in "stats_matrix"
  dim(stats_matrix_test) # the first row is the observed
  
  # AT STEP 4: consider at each rep. T_rep as T_obs to calculate nrep p-values 
  # for the hypothetical test statistics
  
  # hyp_matrix <- stats_matrix_test # remove first row (obs.)
  hyp_p_value <- matrix(NA, ncol = p_test, nrow = dim(stats_matrix_test)[1])
  
  # based on value (hyp_obs) of each row
  for (r in 1:nrep){
    if(verbose)
      if(r%% ceiling(nrep/10) == 1)
        cat(paste0('Testing rep : ',r,'/',nrep,' \n\r'))
    # calc. hypothetical p_value on each column of the matrix 
    hyp_p_value[r,] <- apply(stats_matrix_test, 2, function(x) mean(x >= x[r]))
  }
  
  # for each rep. take the min. p_value
  min_p_nrep <- apply(hyp_p_value, 1, function(x) min(x, na.rm = TRUE))
  
  # calculate the proportion of min_p_nrep that is sm/eq. p_value (for obs.)
  p_value_adj <- sapply(p.values.ratio.normalization, function(x) mean(min_p_nrep <= x))
  
  # p_adj_rejections <- which(p_value_adj <= 0.2)  
  # unname(tax_table(ps_work)[p_adj_rejections,])
  
  rank <- as.factor(rep(rank_names[k],p_test))
  
  one_rank_results <- data.frame(rank,
             raw_test_stat = ratio_matrix_test_raw,
             p_value_nom = p.values.ratio.normalization,
             p_value_adj)
  
  
  dacomp_results <- rbind(dacomp_results, one_rank_results)

}


# save(dacomp_results, file = '/Users/alicesommer/Desktop/DACOMP_cluster/dacomp_results/dacomp_results_PM.RData')




library(NetCoMi)
# library(filematrix)
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(reshape2)

load('/Users/alicesommer/Desktop/Bureau/DOCTORATE/data_pipeline_microbiome/ps_to_net_Gen_smoke.RData')
# load('/Users/alicesommer/Desktop/Bureau/DOCTORATE/data_pipeline_microbiome/ps_to_net_Gen.RData')

load('/Users/alicesommer/Desktop/Bureau/DOCTORATE/data_pipeline_microbiome/net_W_output_Gen_smoke_Oct.RData')
# load('/Users/alicesommer/Desktop/Bureau/DOCTORATE/data_pipeline_microbiome/net_W_output_Gen.RData')

# load('/Users/alicesommer/Desktop/Bureau/DOCTORATE/data_pipeline_microbiome/props_W_Gen_smoke_Oct.RData')


### netCompare ###
load('/Users/alicesommer/Desktop/Bureau/DOCTORATE/data_pipeline_microbiome/comp_W_Gen_10K.RData')
# load('/Users/alicesommer/Desktop/Bureau/DOCTORATE/data_pipeline_microbiome/comp_W_Gen_10K_PM.RData')

summary(comp_W)

metrics <- c('diffClustCoef', 'diffModul', 'diffpnRatio')
Tarray_metrics <- comp_W$permDiffGlobal[,metrics]
colnames(Tarray_metrics) <- c('Clustering coefficient', 'Modularity', '+ edges/- edges')
Tarray_metrics_melt <- melt(Tarray_metrics)
colnames(Tarray_metrics_melt) <- c('','variable','value')

dat_text_lab <- data.frame(variable = colnames(Tarray_metrics))
dat_text_lab$obs_stat <- c(comp_W$diffGlobal$diffClustCoef, comp_W$diffGlobal$diffModul, comp_W$diffGlobal$diffpnRatio)

plot <- ggplot(data = Tarray_metrics_melt, aes(x = value)) + 
  # geom_histogram(aes(y = ..density..) , binwidth=.2) +  
  geom_histogram(binwidth=.01) +  
  facet_wrap(~variable, ncol = 3, scales = 'free') +
  geom_vline(data = dat_text_lab, mapping = aes(xintercept = obs_stat), 
             linetype = "dashed", colour = "red", size = .3) 


### diffnet ###
# load('/Users/alicesommer/Desktop/Bureau/DOCTORATE/data_pipeline_microbiome/diff_net_Gen_10K.RData')
load('/Users/alicesommer/Desktop/Bureau/DOCTORATE/data_pipeline_microbiome/diff_net_Gen_10K_PM.RData')

plot(diff_net, adjusted = FALSE, cexNodes = 0.8, cexLegend = 0.7, mar = c(7,7,7,10))

head(table(diff_net$pvalsVec))
head(table(diff_net$pAdjustVec))

# # load big assoPerm ##
# setwd('/Users/alicesommer/Desktop/')
# 
# diff_net <- diffnet(net_W, nPerm = 1000,
#                     fileLoadAssoPerm = "NetCoMi_cluster/assoPerm/assoPerm",
#                     adjust = "BH")
# 
# plot(diff_net, adjusted = FALSE, cexLegend = .6)
# 
# table(diff_net$pvalsVec)
# table(diff_net$pAdjustVec)

seq_to_seq <- names(diff_net$pvalsVec) 
diff_ass <- seq_to_seq[which(diff_net$pvalsVec < 0.005)]
diff_ass
p_val <- diff_net$pvalsVec[which(diff_net$pvalsVec < 0.005)]

diff_ass_split <- str_split(diff_ass, "_")

for(i in 1:length(diff_ass_split)){
  print(paste('p-value:',p_val[i]))
  # print(diff_ass_split[[i]])
  print(tax_table(ps_Genus_prune)[diff_ass_split[[i]],c("Order","Genus")])
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Multiple comparison adjustment 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
verbose = T
### Lee, Forastiere, Miratix, Pillai (2017) method  for multiple comparison adjustment ### 

# subset data zero edges in both matrices
zero_mat1 <- diff_net$assoMat1 == 0
zero_mat2 <- diff_net$assoMat2 == 0

# keep all the entries where not TRUE + TRUE
mat_FALSE <- zero_mat1 == zero_mat2
# transform matrices to vectors
lowtri <- lower.tri(mat_FALSE, diag = FALSE)
vector_FALSE <- mat_FALSE[lowtri]
# add names in same format than NetCoMi
vec_names <- NetCoMi:::get_vec_names(mat_FALSE)
names(vector_FALSE) <- vec_names

# compare to observed stat = 0 
# because could be that cor same in group 1 and 2 > stat = 0
head(table(diff_net$testStatData))
tail(names(vector_FALSE[vector_FALSE == FALSE]))
tail(names(diff_net$testStatData[diff_net$testStatData != 0]))

identical(names(vector_FALSE[vector_FALSE == FALSE]), names(diff_net$testStatData[diff_net$testStatData != 0]))

subset_edges <- names(vector_FALSE[vector_FALSE == FALSE])
length(subset_edges)

# STEP 1 to 3: recorded in "stats_matrix"
dim(diff_net$testStatPerm)
head(colnames(diff_net$testStatPerm))
# the first row is the observed

# AT STEP 4: consider at each rep. T_rep as T_obs to calculate nrep p-values 
# for the hypothetical test statistics

hyp_matrix <- rbind(diff_net$testStatData[subset_edges], diff_net$testStatPerm[,subset_edges]) 
n_col = dim(hyp_matrix)[2]
n_row = dim(hyp_matrix)[1]

hyp_p_value <- matrix(NA, ncol = n_col, nrow = n_row)

# based on value (hyp_obs) of each row
for (r in 1:n_row){
  if(verbose)
    if(r%% ceiling(n_row/10) == 1)
      cat(paste0('Testing rep : ',r,'/',n_row,' \n\r'))
  # calc. hypothetical p_value on each column of the matrix 
  hyp_p_value[r,] <- apply(hyp_matrix, 2, function(x) mean(x >= x[r]))
}

# check if last row = diff_net$pvalsVec
table(hyp_p_value[1,])
table(diff_net$pvalsVec[subset_edges])

identical(hyp_p_value[1,], as.numeric(diff_net$pvalsVec[subset_edges]))

# for each rep. take the min. p_value
min_p_nrep <- apply(hyp_p_value, 1, function(x) min(x, na.rm = TRUE))
head(min_p_nrep)

# calculate the proportion of min_p_nrep that is sm/eq. p_value (for obs.)
p_value_adj <- sapply(hyp_p_value[1,], function(x) mean(min_p_nrep <= x))
table(p_value_adj)

which(p_value_adj < .3)

diff_net$pvalsVec[subset_edges][which(p_value_adj < .3)]

p_value_adj[which(p_value_adj < .3)]

loc_diff_ass <- which(subset_edges %in% diff_ass)

for(i in 1:length(diff_ass_split)){
  print(paste('p-value:',p_val[i]))
  print(paste('p-value (adj.):',p_value_adj[loc_diff_ass][i]))
  print(tax_table(ps_Genus_prune)[diff_ass_split[[i]],c("Family","Order","Genus")])
}




library(phyloseq); packageVersion("phyloseq")
library(ggplot2)
library(gridExtra)
library(compositions)


dat_prev_var <- NULL

rank_names <- c( "ASV" , "Species" , "Genus"  , "Family" , "Order"  , "Class"  , "Phylum")
 
## create data.frames with prevelance and variance of taxa values
for(k in 1:length(rank_names)){
  
  # load data ps filtered data at 5% prevalence threshold
  ps_work <- readRDS(paste0("/Users/alicesommer/Desktop/DACOMP_cluster/ps_filt_",k,".rds"))
  
  p = ncol(otu_table(ps_work))
  n = nrow(otu_table(ps_work))
  
  ## prevalence
  taxa_prev <- apply(unname(otu_table(ps_work)), 2, function(x) sum(x > 0, na.rm = TRUE)/n)
  
  ## variance
  ps_ait <- transform_sample_counts(ps_work, function(x) {log((x + 1)/sum(x))})
  taxa_var_ait <- apply(unname(otu_table(ps_ait)), 2, function(x) var(x))
  
  ## aggregation level indicator
  rank <- as.factor(rep(rank_names[k],p))
  
  ## append to data.frames
  # prevalence data.frame
  prev_var_rank <- data.frame(rank, taxa_prev, taxa_var_ait)
  # variance data.frame
  dat_prev_var <- rbind(dat_prev_var, prev_var_rank)
  
}

#########################
### SELECT REFERENCES ###
#########################

selected_ref <- NULL

for(r in 1:length(rank_names)){
  ps_work <- readRDS(paste0("/Users/alicesommer/Desktop/DACOMP_cluster/ps_filt_",r,".rds"))
  
  if (r == 1){
    condition <- dat_prev_var$taxa_var_ait[dat_prev_var$rank == rank_names[r]] < 3 & dat_prev_var$taxa_prev[dat_prev_var$rank == rank_names[r]] > .4
    condition <- which(condition)
  } else {
    condition <- dat_prev_var$taxa_var_ait[dat_prev_var$rank == rank_names[r]] < 2 & dat_prev_var$taxa_prev[dat_prev_var$rank == rank_names[r]] > .9
    condition <- which(condition)
  }
  
  # condition <- dat_prev_var$taxa_var_ait[dat_prev_var$rank == rank_names[r]] < 2 & dat_prev_var$taxa_prev[dat_prev_var$rank == rank_names[r]] > .9
  # condition <- which(condition)
  rank <- as.factor(rep(rank_names[r],length(condition)))
  
  ref_select <- data.frame(rank, condition)
  selected_ref <- rbind(selected_ref, ref_select)
  
  print(rank_names[r])
  print(condition)
  
  X = as(otu_table(ps_work), "matrix")
  
  if (length(condition) == 1) {
    reference_values = X[,condition]
  } else {
    reference_values = apply(X[,condition], 1, sum)
  }
  
  print(paste('min. ref. value:', min(reference_values)))
  print(paste('nb. zeros:', sum(reference_values == 0)))
  # print(unname(tax_table(ps_work)[condition,]))
}

dat_prev_var$col <- dat_prev_var$taxa_var_ait < 2 & dat_prev_var$taxa_prev > .9
# special case for ASV 
dat_prev_var$col[dat_prev_var$rank == "ASV"] <- dat_prev_var$taxa_var_ait[dat_prev_var$rank == "ASV"] < 3 & dat_prev_var$taxa_prev[dat_prev_var$rank == "ASV"] > .4

## plot the reference selection process at all rank levels
g_1 <- ggplot(data = dat_prev_var) +
  geom_point(aes(x = taxa_prev, y = taxa_var_ait, color = col), size = .3) +
  facet_wrap(rank~., ncol = 2) + 
  theme(legend.title=element_blank()) +
  ylab("Variance") + xlab("Prevalence") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "darkgrey"), 
                     labels = c(c("TRUE" = "reference", "FALSE" = "test")))

ggsave(file = "/Users/alicesommer/Desktop/Bureau/DOCTORATE/plots_pipeline_microbiome/selected_reference_set.png",
     g_1,
     dpi=300,
     width = 150,
     height = 160,
     units = "mm")

# save(selected_ref, file = "/Users/alicesommer/Desktop/DACOMP_cluster/dacomp_results/selected_ref_PM.RData")

# #%%%%%%%%%%%%
# # Subsample #
# #%%%%%%%%%%%%
# 
# splits = 4
# 
# p = ncol(otu_table(ps_work))
# n = nrow(otu_table(ps_work))
# 
# prev_mat <- matrix(NA, ncol = splits, nrow = p)
# var_mat <- matrix(NA, ncol = splits, nrow = p)
# 
# for(i in 1:splits){
#   sample_pairs <- sample(unique(sample_data(ps_work)$pair_nb), n/splits)
#   ps_split <- subset_samples(ps_work, sample_data(ps_work)$pair_nb %in% sample_pairs)
#   
#   # prevalence
#   n_split <- dim(otu_table(ps_split))[1]
#   prev_mat[,i] <- apply(unname(otu_table(ps_split)), 2, function(x) sum(x > 0, na.rm = TRUE)/n_split)
#   
#   # variance aitchison
#   ps_ait_split <- transform_sample_counts(ps_split, function(x) {log((x + 1)/sum(x))})
#   var_mat[,i] <- apply(unname(otu_table(ps_ait_split)), 2, function(x) var(x))
#   
#   condition_split <- which(var_mat[,i] < 1.5 & prev_mat[,i] > .99)
#   print(condition_split)
#   print(unname(tax_table(ps_split)[condition_split,]))
# }
# 
# ################
# # plot results #
# ################
# 
# condition <- which(taxa_var_ait < 1.5 & taxa_prev > .99)
# condition
# 
# unname(tax_table(ps_work)[condition,])
# 
# par(mfrow = c(1,1))
# plot(taxa_prev, taxa_var_clr)
# par(mfrow = c(1,1))
# plot(taxa_prev, taxa_var_ait)
# abline(h = 1.5, col = 'red')
# abline(v = .99, col = 'red')
# 
# par(mfrow = c(2,2))
# plot(prev_mat[,1], var_mat[,1])
# abline(h = 1.5, col = 'red')
# abline(v = .99, col = 'red')
# plot(prev_mat[,2], var_mat[,2])
# abline(h = 1.5, col = 'red')
# abline(v = .99, col = 'red')
# plot(prev_mat[,3], var_mat[,3])
# abline(h = 1.5, col = 'red')
# abline(v = .99, col = 'red')
# plot(prev_mat[,4], var_mat[,4])
# abline(h = 1.5, col = 'red')
# abline(v = .99, col = 'red')


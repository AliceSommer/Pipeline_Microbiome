library(phyloseq); packageVersion("phyloseq")
library(reshape2)
library(ggplot2)
library(dacomp)

# set working directory
setwd('/Users/alicesommer/Desktop/Bureau/DOCTORATE')

# load microbiome data
ASV_table <- readRDS('data_pipeline_microbiome/dada2output/seqtab2020.rds')
taxon_assign <- readRDS('data_pipeline_microbiome/dada2output/taxa2020.rds')

# load sample/matched_data
load('data_pipeline_microbiome/dat_matched_PM25_bis.RData')

sample_df <- matched_df[order(matched_df$ff4_prid),]
sample_df$W <- as.factor(sample_df$W)
samples.out <- as.character(sample_df$ff4_prid)
rownames(sample_df) <- samples.out

# create a phyloseq object
ps <- phyloseq(otu_table(ASV_table, taxa_are_rows=FALSE),
               sample_data(sample_df),
               tax_table(taxon_assign))
ps

rank_names(ps)

# For permutation load W matrix for randomization test
load("data_pipeline_microbiome/W_paired_PM25.Rdata")

#### RESOLUTION ####

perc_seq <- seq(0,0.15,0.05)
ASV_cut <- matrix(NA, nrow = 1, ncol = length(perc_seq))
p_val_mat <- matrix(NA, nrow = dim(otu_table(ps))[2], ncol = length(perc_seq))
p_val_dacomp <- matrix(NA, nrow = dim(otu_table(ps))[2], ncol = length(perc_seq))
ASV_nr <- matrix(NA, nrow = 1, ncol = length(perc_seq))
ref_size <- matrix(NA, nrow = 1, ncol = length(perc_seq))
ref_abund_min <- matrix(NA, nrow = 1, ncol = length(perc_seq))
stats_matrix_test <- list()
ratio_matrix_test <- list()

vec_taxa <- apply(otu_table(ps), 2, function(x) sum(x > 0, na.rm = TRUE))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## research variable
Y = sample_data(ps)$W 
## Permutated "W" matrix
nrep <- 10^4
Y_perm <- cbind(as.numeric(Y == 1),W_paired[,1:nrep]) # the first column of the Y_perm matrix is W_obs

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Parameters for method:
samples_for_reference = round(length(Y)*0.15) # how many samples should be taken for reference

############
## DACOMP ##
############

for (i in 1:length(perc_seq)){
  print(i)
     
  # how many tax are obs. in at least x% of samples
  perc <- perc_seq[i] # x%
  ASV_cut[i] <- trunc(dim(sample_df)[1]*perc)
  
  ASV_nr[i] <- length(which(vec_taxa > ASV_cut[i]))
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Data
  X = as(otu_table(ps), "matrix")
  # add NA for all taxa that do not satisfy the filtering 
  filt = unname(which(vec_taxa <= ASV_cut[i]))
  X[,filt] <- matrix(NA, nrow = nrow(X), ncol = length(filt))
  
  taxa_to_test <- which(vec_taxa > ASV_cut[i])
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  set.seed(16)
  condition = sample(1:samples_for_reference)

  X_reference_select = X[condition,]
  X_test = X[-c(condition),]
  Y_reference_select = Y[condition]
  Y_test = Y[-c(condition)]
  Y_perm_reference_select = Y_perm[condition,]
  Y_perm_test = Y_perm[-c(condition),]

  #definitions and allocations:
  p = ncol(X_reference_select)
  n = nrow(X_reference_select)

  ratio_matrix = matrix(NA, nrow = n, ncol = 1)
  stats_matrix = matrix(NA, ncol = p, nrow = nrep+1)

  verbose = T

  # TEST = DACOMP.TEST.NAME.WILCOXON
  # TEST = DACOMP.TEST.NAME.DIFFERENCE_IN_MEANS
  TEST = DACOMP.TEST.NAME.LOG_FOLD_DIFFERENCE_IN_MEANS
  # TEST = DACOMP.TEST.NAME.WELCH
  # TEST = DACOMP.TEST.NAME.WELCH_LOGSCALE

  #%%%%%%%%%%%%%%%%%%%%
  # SELECT REFERENCES
  #%%%%%%%%%%%%%%%%%%%

  #iterate over taxa and test
  for(a in 1:length(taxa_to_test) ){

    if(verbose)
     if(a %% ceiling(p/10) == 1)
       cat(paste0('Testing taxon : ',a,'/',length(taxa_to_test),' \n\r'))

    t = taxa_to_test[a]

    nom = X_reference_select[,t]
    dnom = apply(X_reference_select, 1, function(x) sum(x, na.rm = TRUE)) # because we don't have a reference set yet

    #perform subsample and test
    ratio_matrix[,1] = nom/(dnom)
    stats_matrix[,t] = dacomp:::Compute.resample.test(ratio_matrix, Y_perm_reference_select, statistic = TEST)

    }

   #computes pvalues:
   p_val_mat[,i] = apply(stats_matrix, 2, function(x) mean(x >= x[1], na.rm = TRUE))

   # Select the reference set of taxa using the above P-values
   ind_reference_taxa = which((p_val_mat[,i] >= .5))

   # Sanity checks:
   ref_size[i] <- length(ind_reference_taxa)
   # abundance check
   counts_ref <- apply(X[,ind_reference_taxa], 1,function(x) sum(x, na.rm = TRUE))
   ref_abund_min[i] <- summary(counts_ref)[1]

   #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   # TEST WITH DACOMP-RATIO METHOD
   #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   p_test = ncol(X_test)
   n_test = nrow(X_test)
   
   ratio_matrix_test[[i]] = matrix(NA, ncol = p_test, nrow = n_test)
   ratio_matrix_test_function = matrix(NA, nrow = n_test, ncol = 1)
   
   stats_matrix_test[[i]] = matrix(NA, ncol = p_test, nrow = nrep+1)

   # Compute reference values
   reference_values = apply(X_test[,ind_reference_taxa], 1, sum)

   #iterate over taxa and test
   for(a in 1:length(taxa_to_test)){

     if(verbose)
       if(a %% ceiling(p/10) == 1)
         cat(paste0('Testing taxon : ',a,'/',length(taxa_to_test),' \n\r'))
     
     t = taxa_to_test[a]

     nom_test = X_test[,t]
     dnom_test = reference_values

     ratio_matrix_test[[i]][,t] = nom_test/(dnom_test+nom_test)
     ratio_matrix_test_function[,1] = nom_test/(dnom_test+nom_test)

     stats_matrix_test[[i]][,t] = dacomp:::Compute.resample.test(ratio_matrix_test_function, 
                                                                 Y_perm_test, statistic = TEST)

   }

   #computes pvalues:
   p_val_dacomp[,i] = apply(stats_matrix_test[[i]], 2, function(x) mean(x >= x[1]))
   
}

head(sort(p_val_mat),20)
which(p_val_mat[,1] < 0.06)

## multiple comparison adjustment
head(p.adjust(p_val_dacomp[,1],method = 'BH'))
head(sort(p.adjust(p_val_dacomp[,1],method = 'BH')))
which(p.adjust(p_val_dacomp[,1],method = 'BH') < .1)
len_1 <- length(which(p.adjust(p_val_dacomp[,1],method = 'BH') < .1))

head(sort(p.adjust(p_val_dacomp[,2],method = 'BH')))
which(p.adjust(p_val_dacomp[,2],method = 'BH') < .1)

head(sort(p.adjust(p_val_dacomp[,3],method = 'BH')))
which(p.adjust(p_val_dacomp[,3],method = 'BH') < .1)

head(sort(p.adjust(p_val_dacomp[,4],method = 'BH')))
which(p.adjust(p_val_dacomp[,4],method = 'BH') < .1)

unname(tax_table(ps)[286,])

dat_hist_1 <- melt(ratio_matrix_test[[1]][,which(p.adjust(p_val_dacomp[,1],method = 'BH') < .1)])
dat_hist_1$W <- rep(Y_test, len_1)

ggplot(data=dat_hist_1, aes(x = value, fill = W)) +
  geom_histogram() +
  facet_wrap(. ~ Var2, scales = "free")

dat_hist <- data.frame(asv = ratio_matrix_test[[4]][,286],
                       W = Y_test)

ggplot(data=dat_hist, aes(x = asv, fill = W)) +
  geom_histogram(binwidth = 0.001)

###########################
#### PLOTS FOR PVALUES ####
###########################

low_high <- c(order(p_val_dacomp[,2])[1:10], 
              which(p_val_dacomp[,2] > .9)[1:10])

dat_low_high <- data.frame(p_val_dacomp[low_high,2:4])
tail(dat_low_high)

dat_low_high$facet <- as.numeric(dat_low_high[,1] < 0.03)
dat_low_high$id <- low_high
dim(dat_low_high)

dat_plot_low_high <- melt(dat_low_high, measure.vars = 1:3, 
                          id.vars = c("facet","id"))

prev_labels <- list('1'="10 high p-values",'0'="10 lowest p-values")

prev_labeller <- function(variable,value){
  return(prev_labels[value])
}


g_sens <- ggplot(dat_plot_low_high, aes(y = value, x = as.factor(variable), 
                         color = as.factor(id), group = as.factor(id))) +
  geom_line(size = .5) + geom_point(size = .5) +
  scale_color_discrete(name = "ASV") +
  scale_x_discrete(name ="prevalence filtering threshold",
                   labels = c("X1" = "5%", 
                              "X2" = "10%",
                              "X3" = "15%")) +
  facet_wrap(~ as.factor(facet), ncol = 2, 
             scales = "free_y", labeller=prev_labeller) +
  theme(legend.position="bottom")

# ggsave(file = 'Pipeline_Microbiome/4_differential_abundance/sensitivity_diff_ab_pval.jpeg',
#        g_sens,
#        dpi=300,
#        width = 100,
#        height = 90,
#        units = "mm")

### check the low p-values ASVs

col_num <- which(p.adjust(p_val_dacomp[,2],method = 'BH') < .1)

unname(tax_table(ps)[col_num,])

rownames(tax_table(ps)[col_num,])

## reference size plot
g_ref <- ggplot() + 
  geom_point(aes(y = ref_size[1,]/ASV_nr[1,], x = perc_seq)) +
  geom_line(aes(y = ref_size[1,]/ASV_nr[1,], x = perc_seq)) +
  xlab('prevalence filtering threshold') +
  ylab('# references/# ASVs') + 
  scale_x_continuous(name ="prevalence filtering threshold", 
                   breaks = c(0.00,0.05,0.10,0.15),
                   labels = c('0%', '5%', '10%', '15%'))

# ggsave(file = 'Pipeline_Microbiome/4_differential_abundance/sensitivity_diff_ab_ref.jpeg',
#        g_ref,
#        dpi=300,
#        width = 80,
#        height = 60,
#        units = "mm")

hist(stats_matrix_test[[2]][-1,col_num], breaks = 30, main = "", xlab = "log fold diff. in means")
abline(v = stats_matrix_test[[2]][1,col_num], col = 'red', lwd = 2, lty = 2)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Multiple comparison adjustment 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### Lee, Forastiere, Miratix, Pillai (2017) method  for multiple comparison adjustment ### 

# STEP 1 to 3: recorded in "stats_matrix"
dim(stats_matrix_test[[2]])
# the first row is the observed

# AT STEP 4: consider at each rep. T_rep as T_obs to calculate nrep p-values 
# for the hypothetical test statistics

hyp_matrix <- stats_matrix_test[[2]][-1,] # remove first row (obs.)
hyp_p_value <- matrix(NA, ncol = p_test, nrow = nrep)

# based on value (hyp_obs) of each row
for (r in 1:nrep){
  if(verbose)
    if(r%% ceiling(nrep/100) == 1)
      cat(paste0('Testing rep : ',r,'/',nrep,' \n\r'))

  # calc. hypothetical p_value on each column of the matrix 
  hyp_p_value[r,] <- apply(hyp_matrix, 2, function(x) mean(x >= x[r]))
}

# for each rep. take the min. p_value
min_p_nrep <- apply(hyp_p_value, 1, function(x) min(x, na.rm = TRUE))
head(min_p_nrep,20)

hist(min_p_nrep, breaks = 60, main = "", xlab = "minimun p-value for each rep.")
abline(v = p.values.ratio.normalization[117], col = "red")
text(0.02, y = 1190, 0.00019998, col = "red")
text(0.06, y = 800, "adj. p-value: 0.0064", col = "darkgreen")

# calculate the proportion of min_p_nrep that is sm/eq. p_value (for obs.)
p_value_adj <- sapply(p_val_dacomp[,2], function(x) mean(min_p_nrep <= x))
head(p_value_adj)

head(sort(p_value_adj),50) 
p_adj_rejections <- which(p_value_adj <= 0.2)  
p_adj_rejections

unname(tax_table(ps_work)[p_adj_rejections,])
  

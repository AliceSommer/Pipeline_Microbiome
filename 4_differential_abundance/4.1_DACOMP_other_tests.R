# library(coin)
library(dacomp)
# library(Hmisc)
# library(doRNG)
library(phyloseq); packageVersion("phyloseq")
# library(openxlsx)

# set working directory
setwd('/Users/alicesommer/Desktop/Bureau/DOCTORATE')

# This script gives an example for split testing using dacomp.
# source('Pipeline_Microbiome/misc/Functions_for_univariate_tests.R')
# source('Pipeline_Microbiome/misc/dacomp_testing_and_reference_selection_by_split.R')

###############################################################################

# load microbiome data
ASV_table <- readRDS('data_pipeline_microbiome/dada2output/seqtab2020.rds')
taxon_assign <- readRDS('data_pipeline_microbiome/dada2output/taxa2020.rds')

# load sample/matched_data
load('data_pipeline_microbiome/dat_matched_PM25_bis.RData')

# # load original dataset (try other Ws)
# df <- read.csv('/Volumes/GoogleDrive/My\ Drive/DOCTORATE/Thesis/KORA\ DATA/Microbiome_data/KORA_microbiome_variables.csv')

# matched AP data
# matched_df <- matched_data
sample_df <- matched_df[order(matched_df$ff4_prid),]
sample_df$W <- as.factor(sample_df$W)

# # "other" data
# sample_df <- df[order(df$ff4_prid),]
# sample_df$W <- as.factor(as.numeric(sample_df$u3tcigsmk == 1))
samples.out <- as.character(sample_df$ff4_prid)
rownames(sample_df) <- samples.out

# For permutation load W matrix for randomization test
load("data_pipeline_microbiome/W_paired_PM25.Rdata")

################################################################################

# create a phyloseq object
ps <- phyloseq(otu_table(ASV_table, taxa_are_rows=FALSE),
               sample_data(sample_df),
               tax_table(taxon_assign))
ps

rank_names(ps)

# locate the species that are totally absent in the matched data 
# 1. have a vector of nr. of observed samples per taxa
vec_taxa <- apply(otu_table(ps), 2, function(x) sum(x > 0, na.rm = TRUE))
length(which(vec_taxa == 0))
# 2. how many tax are obs. in at least x% of samples
perc <- 0.05 # x%
samp_perc <- trunc(dim(sample_df)[1]*perc)
# samp_perc <- 0
length(which(vec_taxa > samp_perc))

ps_prune <- prune_taxa(vec_taxa >= 1, ps)

## agglomerate to Species ##
ps_Species <- tax_glom(ps, taxrank = "Species", NArm = FALSE)
vec_taxa_Sp <- apply(otu_table(ps_Species), 2, function(x) sum(x > 0, na.rm = TRUE))
length(which(vec_taxa_Sp > samp_perc))
ps_Species <- prune_taxa(vec_taxa_Sp >= samp_perc, ps_Species)
ps_Species

## agglomerate to Genus ##
ps_Genus <- tax_glom(ps_prune, taxrank = "Genus", NArm = FALSE)
vec_taxa_Gen <- apply(otu_table(ps_Genus), 2, function(x) sum(x > 0, na.rm = TRUE))
table(vec_taxa_Gen)
ps_Genus_prune <- prune_taxa(vec_taxa_Gen >= samp_perc, ps_Genus)
ps_Genus_prune

## agglomerate to Family ##
ps_Family <- tax_glom(ps_prune, taxrank = "Phylum", NArm = FALSE)

ps_work <- ps_Species

################################################################################

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Data
X = as(otu_table(ps_work), "matrix")
Y = sample_data(ps_work)$W # research variable
# Permutated"W" matrix
nrep <- 10^4
# dim(W_paired)
# length(Y)
Y_perm <- cbind(as.numeric(Y == 1),W_paired[,1:nrep]) # the first column of the Y_perm matrix is W_obs

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set.seed(16)
# Parameters for method:
samples_for_reference = round(nrow(X)*0.15) # how many samples should be taken for reference

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

stats_matrix = matrix(NA, ncol = p, nrow = nrep+1)
ratio_matrix = matrix(NA, nrow = n, ncol = 1)

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
for(i in 1:p){
  
  if(verbose)
    if(i%% ceiling(p/10) == 1)
      cat(paste0('Testing taxon : ',i,'/',p,' \n\r'))
  
  nom = X_reference_select[,i]
  dnom = apply(X_reference_select,1,sum) # because we don't have a reference set yet
  
  #perform subsample and test
  ratio_matrix[,1] = nom/(dnom)
  stats_matrix[,i] = dacomp:::Compute.resample.test(ratio_matrix, Y_perm_reference_select, statistic = TEST)
  
}

#computes pvalues:
p.values = apply(stats_matrix, 2, function(x) mean(x >= x[1], na.rm = TRUE))

head(sort(p.values),20)

# Select the reference set of taxa using the above P-values
ind_reference_taxa = which((p.values >= .5))

# Sanity checks:
length(ind_reference_taxa)
# abundance check
counts_ref <- apply(otu_table(ps_work)[,ind_reference_taxa], 1,function(x) sum(x, na.rm = TRUE))
summary(counts_ref)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# TEST WITH DACOMP-RATIO METHOD 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_test = ncol(X_test)
n_test = nrow(X_test)
stats_matrix_test = matrix(NA, ncol = p_test, nrow = nrep+1)
ratio_matrix_test = matrix(NA, nrow = n_test, ncol = 1)

# Compute reference values
reference_values = apply(X_test[,ind_reference_taxa], 1, sum)

#iterate over taxa and test
for(t in 1:p_test){
  
  if(verbose)
    if(t%% ceiling(p/10) == 1)
      cat(paste0('Testing taxon : ',t,'/',p,' \n\r'))
  
  #no need to test reference taxa
  if(t %in% ind_reference_taxa){
    print(t)
    next
  }

  nom_test = X_test[,t]
  dnom_test = reference_values
  
  ratio_matrix_test[,1] = nom_test/(dnom_test+nom_test)
  
  stats_matrix_test[,t] = dacomp:::Compute.resample.test(ratio_matrix_test,Y_perm_test, statistic = TEST)
  
}

#computes pvalues:
p.values.ratio.normalization = apply(stats_matrix_test, 2, function(x) mean(x >= x[1]))

head(sort(p.values.ratio.normalization),20)

which(p.adjust(p.values.ratio.normalization,method = 'BH')<=0.2)

head(sort(p.adjust(p.values.ratio.normalization,method = 'BH')),20)

unname(tax_table(ps_work)[which(p.adjust(p.values.ratio.normalization,method = 'BH')<=0.21),])
# tax_table(ps_work)[which(p.adjust(p.values.ratio.normalization,method = 'BH')<=0.1),]

hist(stats_matrix_test[,117], xlab = "T_obs", main = "")
abline(v = stats_matrix_test[1,117], col = "red")

counts_discov <- apply(otu_table(ps_work)[,which(p.adjust(p.values.ratio.normalization,method = 'BH')<=0.21)], 2,function(x) sum(x > 0, na.rm = TRUE))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Multiple comparison adjustment 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### Lee, Forastiere, Miratix, Pillai (2017) method  for multiple comparison adjustment ### 

# STEP 1 to 3: recorded in "stats_matrix"
dim(stats_matrix_test)
# the first row is the observed

# AT STEP 4: consider at each rep. T_rep as T_obs to calculate nrep p-values 
# for the hypothetical test statistics

hyp_matrix <- stats_matrix_test[-1,] # remove first row (obs.)
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
p_value_adj <- sapply(p.values.ratio.normalization, function(x) mean(min_p_nrep <= x))
head(p_value_adj[117])

head(sort(p_value_adj),50) 
p_adj_rejections <- which(p_value_adj <= 0.1)  
p_adj_rejections

unname(tax_table(ps_work)[p_adj_rejections,])

# # compute DS-FDR: >>>>>>>> NOT WORKING YET 
# DSFDR_Filter = rep(T,ncol(X))
# disable_DSFDR = F
# 
# ## ADD RELEVANT CODE FOR DSFDR to work

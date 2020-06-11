library(dacomp)
library(corncob)

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
empty_species <- colSums(otu_table(ps))
length(which(empty_species == 0))

ps_prune <- prune_taxa(empty_species != 0, ps)

## agglomerate to Genus ##
ps_Genus <- tax_glom(ps_prune, taxrank = "Genus", NArm = FALSE)
ps_Genus 

## agglomerate to Family ##
ps_Family <- tax_glom(ps_prune, taxrank = "Family", NArm = FALSE)
ps_Family 

ps_work = ps_Genus

################################################################################

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Data
X = as(otu_table(ps_work), "matrix")
Y = sample_data(ps_work)$W # research variable
# Permutated"W" matrix
nrep <- 10^3
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
p.values = apply(stats_matrix, 2, function(x) mean(x >= x[1]))

head(sort(p.values),20)

# Select the reference set of taxa using the above P-values
ind_reference_taxa = which(p.values >= .5)

# Sanity checks:
length(ind_reference_taxa)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# TEST WITH DACOMP-RATIO METHOD 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_test = ncol(X_test)
n_test = nrow(X_test)

stats_matrix_test = matrix(NA, ncol = p_test, nrow = nrep+1)

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
  
  for (l in 1:(nrep+1)){ # the first iteration is to calc. the observed p_value
    # corncob 
    df_cb <- data.frame("W" = nom_test,
                        "M" = nom_test+dnom_test,
                        "X1" = as.factor(Y_perm_test[,l]))
    
    model = bbdml(formula = cbind(W, M - W) ~ X1,
                  phi.formula = ~ X1,
                  data = df_cb);
    model_summary = summary(model)
    
    
    stats_matrix_test[l,t] = model_summary$coefficients["mu.X11","Estimate"]
  }

}

#computes pvalues:
p.values.ratio.normalization = apply(stats_matrix_test, 2, function(x) mean(x >= x[1]))

head(sort(p.values.ratio.normalization),20)

which(p.adjust(p.values.ratio.normalization,method = 'BH')<=0.1)

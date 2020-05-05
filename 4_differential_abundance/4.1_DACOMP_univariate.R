library(coin)
library(dacomp)
library(Hmisc)
library(doRNG)
library(phyloseq); packageVersion("phyloseq")
# library(openxlsx)

# set working directory
setwd('/Users/alicesommer/Desktop/Bureau/DOCTORATE')

# This script gives an example for split testing using dacomp.
source('Pipeline_Microbiome/misc/Functions_for_univariate_tests.R')
source('Pipeline_Microbiome/misc/dacomp_testing_and_reference_selection_by_split.R')

###############################################################################

# load microbiome data
ASV_table <- readRDS('data_pipeline_microbiome/dada2output/seqtab2020.rds')
taxon_assign <- readRDS('data_pipeline_microbiome/dada2output/taxa2020.rds')

# load sample/matched_data
load('data_pipeline_microbiome/dat_matched_PM25.RData')

sample_df <- matched_df[order(matched_df$ff4_prid),]
sample_df$W <- as.factor(sample_df$W)
samples.out <- as.character(sample_df$ff4_prid)
rownames(sample_df) <- samples.out

################################################################################

# create a phyloseq object
ps <- phyloseq(otu_table(ASV_table, taxa_are_rows=FALSE),
               sample_data(sample_df),
               tax_table(taxon_assign))
ps

rank_names(ps)

## agglomerate to Genus ##
ps_Genus <- tax_glom(ps, taxrank = "Genus")
ps_Genus # 269 taxa

## agglomerate to Genus ##
ps_Family <- tax_glom(ps, taxrank = "Family")
ps_Family

ps_work = ps_Family

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set.seed(16)
# Parameters for method:
samples_for_reference = 50 # how many samples should be taken for reference

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Data
X = as(otu_table(ps_Genus), "matrix")
Y = sample_data(ps_Genus)$W # research variable
# Z = as.factor(sample_data(ps)$u3csex) # strata
Z = as.factor(sample_data(ps_Genus)$pair_nb) # use Z as pairs for the paired test

# we sample some pairs to keep a balanced set
pair_id = unique(sample_data(ps_Genus)$pair_nb)
X = as(otu_table(ps_work), "matrix")
Y = sample_data(ps_work)$W # research variable
# Z = as.factor(sample_data(ps)$u3csex) # strata
Z = as.factor(sample_data(ps_work)$pair_nb) # use Z as pairs for the paired test

# we sample some pairs to keep a balanced set
pair_id = unique(sample_data(ps_work)$pair_nb)
pair_id_sample = sample(pair_id, replace = F, size = samples_for_reference) 

## Note: code only if you want to test specific taxa
# ####### DATA/TAXA to test #######
# med_IFG <- read.xlsx('Microbiome\ 2020/4.\ April/DACOMP_Split_AJS_Edition-master/meditation_analysis_Liu.xlsx', sheet = 1, rows = 3:82)
# head(med_IFG)
# 
# med_DT2 <- read.xlsx('Microbiome\ 2020/4.\ April/DACOMP_Split_AJS_Edition-master/meditation_analysis_Liu.xlsx', sheet = 2, rows = 3:87)
# head(med_DT2)
# str(med_DT2)
# 
# med_species <- med_IFG[,"Species"]
# med_species_D <- med_DT2[,"Species"]
# regex_sub <- '[a-z]__[aA-zZ]'
# 
# locate_sp <- tax_table(ps)[,"Species"] %in% sub('[a-z]__', "", med_species[grep(regex_sub, med_species)])
# taxa_id <- which(locate_sp == TRUE)
# taxa_id
# 
# locate_sp_D <- tax_table(ps)[,"Species"] %in% sub('[a-z]__', "", med_species_D[grep(regex_sub, med_species_D)])
# taxa_id_unique <- unique(c(taxa_id, which(locate_sp_D == TRUE)))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Data analysis starts from here:

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step 1: Split data to reference selection set and test set

condition = sample_data(ps_Genus)$pair_nb %in% pair_id_sample # warning: can only do that because X, Y, Z have same sample_id order

X_reference_select = X[condition,]
X_test = X[!condition,]
Y_reference_select = Y[condition]
Y_test = Y[!condition]
Z_reference_select = Z[condition]
Z_test = Z[!condition]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step 2: Run marginal tests (TSS normalization) on the reference selection data
start.time = Sys.time()
pvals_marginal_result = dacomp.test.with.strata(X = X_reference_select,
                                                Y = Y_reference_select,
                                                Z = Z_reference_select,
                                                Method = 'Wilcoxon-Paired',  # 'Wilcoxon-Strata-Asymp' for asymptotic test over strata, 'C_Wilcoxon' for (block) permutation based
                                                do.block.mean.normalization = T,
                                                nr.perm = 10000, Minimum_Block_Size = 2,
                                                normalize_by_DACOMP_ratio = T,
                                                run.in.parallel = T)
end.time = Sys.time()
end.time - start.time
plot(-log(pvals_marginal_result$P.values),main = 'log Pval for reference selection cohort')
abline(h = -log(0.5),col = 'red')

## Note: code only if you want to test specific taxa
# avoid selecting the taxa to test in my reference set
# set to zero
# pvals_marginal_result$P.values[taxa_id_unique] = 0

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step 3: Select the reference set of taxa using the above P-values
reference_obj = dacomp.select_references.by.split(P_vals_marginal = pvals_marginal_result,
                                                  Counts_Mat_testing = X_test,
                                                  REFERENCE_SELECTION_PVAL_THRESHOLD = 0.5,
                                                  MINIMAL_NR_SUBJECTS_FOR_REFERENCE_TEST = 5,
                                                  MINIMAL_TA = 50,
                                                  MAXIMAL_TA = 1000)

# Sanity checks:
reference_obj$selected_MinAbundance
length(reference_obj$selected_references)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step 4: Test using the selected reference set
start.time = Sys.time()
pvals_DACOMP_RATIO = dacomp.test.with.strata(X = X_test,
                                             Y = Y_test,
                                             Z = Z_test,
                                             taxa_to_normalize_by = reference_obj,
                                             Method = 'Wilcoxon-Paired', # 'Wilcoxon-Strata-Asymp' for asymptotic test over strata, 'C_Wilcoxon' for (block) permutation based
                                             do.block.mean.normalization = F,
                                             nr.perm = 10000, Minimum_Block_Size = 2,
                                             normalize_by_DACOMP_ratio = T,
                                             run.in.parallel = T)

end.time = Sys.time()
end.time - start.time

DACOMP_discoveries = which(p.adjust(pvals_DACOMP_RATIO$P.values,method = 'BH')<=0.1)
DACOMP_discoveries

DSFDR_rejections = which(pvals_DACOMP_RATIO$DSFDR.AdjustedPvalues<=0.1)
DSFDR_rejections

# sanity check:
# FDR = sum(DACOMP_discoveries>m1)/max(length(DACOMP_discoveries),1) ; FDR
# TP = sum(DACOMP_discoveries<=m1) ; TP

## Note: code only if you want to test specific taxa
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Step 5: Test, but only over 20 of the taxa.
# #The first 10 taxa will have lower DSFDR adjusted P-values, if we can focus our reseach question
# #to a smaller familiy of hypotheses with less signals
# 
# DSFDR_filter = rep(F,ncol(X))
# DSFDR_filter[taxa_id_unique] = T
# pvals_DACOMP_RATIO_with_Filter = dacomp.test.with.strata(X = X_test,
#                                                          Y = Y_test,
#                                                          Z = Z_test,
#                                                          taxa_to_normalize_by = reference_obj,
#                                                          Method = 'Wilcoxon-Paired',
#                                                          do.block.mean.normalization = T,
#                                                          nr.perm = 10000, Minimum_Block_Size = 2,
#                                                          normalize_by_DACOMP_ratio = T,
#                                                          run.in.parallel = F,
#                                                          select.taxa.for.DSFDR = DSFDR_filter)
# 
# #Compare DS-FDR p-values when testing:
# pvals_DACOMP_RATIO$DSFDR.AdjustedPvalues[taxa_id_unique] # over all taxa
# pvals_DACOMP_RATIO_with_Filter$P.values[taxa_id_unique]
# pvals_DACOMP_RATIO_with_Filter$DSFDR.AdjustedPvalues[taxa_id_unique] # over a specific set of taxa
# 
# #############################
# # Subset some taxa to test #
# ############################
# 
# head(tax_table(ps))
# dim(tax_table(ps))
# colnames(tax_table(ps))
# 
# # species only
# head(tax_table(ps)[,"Order"],50)
# 
# grep('copri', tax_table(ps)[,"Species"])
# 
# pvals_DACOMP_RATIO$P.values[grep('Prevotella', tax_table(ps)[,"Genus"])]
# 
# sub_small_p <- subset(tax_table(ps),pvals_DACOMP_RATIO$P.values < .059)
# dim(sub_small_p)
# 
# ##############################
# 
# pvals_DACOMP_RATIO$P.values[unique(taxa_id)]
# 
# unique(taxa_id)[which(p.adjust(pvals_DACOMP_RATIO$P.values[unique(taxa_id)],method = 'BH')<=0.1)]
# 
# tax_table(ps)[230,]
# 
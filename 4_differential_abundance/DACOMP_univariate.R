library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(dacomp)

###############################################################################

# set working directory
setwd('/Users/alicesommer/Desktop/Bureau/DOCTORATE/data_pipeline_microbiome')

# load microbiome data
ASV_table <- readRDS('dada2output/seqtab2020.rds')
taxon_assign <- readRDS('dada2output/taxa2020.rds')

# load sample/matched_data
load('dat_matched_PM25.RData')

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

otu_matrix = as(otu_table(ps), "matrix")

result.selected.references = dacomp.select_references(X = otu_matrix,
                                                      median_SD_threshold = 0.5, #APPLICATION SPECIFIC
                                                      verbose = F)
print(result.selected.references)


dacomp.plot_reference_scores(result.selected.references)
?dacomp.select_references

### DACOMP

#multiplicity correction levels for the BH and DS-FDR methods
q_BH = q_DSFDR = 0.1

#Perform testing:
result.test = dacomp.test(X = otu_matrix, #counts data
                          y = sample_data(ps)$W, #phenotype in y argument
                          # obtained from dacomp.select_references(...):
                          ind_reference_taxa = result.selected.references,
                          test = DACOMP.TEST.NAME.WILCOXON, #constant, name of test
                          verbose = F, q = q_DSFDR, nr_perm = 10000, # multiplicity adjustment level
                          compute_ratio_normalization = T)

#These are the indices of taxa discoverted as differentially abundant:
# by applying a BH multiplicity adjustment on the P-values:
rejected_BH = which(p.adjust(result.test$p.values.test.ratio.normalization, method = 'BH') <= q_BH)
rejected_BH = which(result.test$p.values.test.ratio.normalization <= .05) # not adjusted BH, nor FDR

#by applying a DS-FDR multiplicity adjustment on the P-values:
rejected_DSFDR = result.test$dsfdr_rejected

print(result.test)

### Corncob
### don't know if I want to use that for my project.................

# taxa_to_test = c(1:ncol(otu_matrix))
# 
# #We can have two types of models:
# #"DACOMP type models" - models for the number of counts belonging to the tested taxon, when observing lambda counts belonging to the tested taxon and the reference
# #(i.e., rarefy subvector (X_j,X_reference), j being the tested taxon, to equal number of reads acroos samples)
# P.values.rarefaction = rep(NA,length(taxa_to_test))
# # "DACOMP-ratiotype models" - models for the ratio of reads belong to the tested taxon, when observing the total number of reads available to the tested taxon and the reference set.
# #(i.e., model subvector (X_j,X_reference) directly, j being the tested taxon, with a different number of reads per subject)
# P.values.ratio = rep(NA,length(taxa_to_test))
# 
# # Iterate over taxa:
# for(taxon_id in 1:length(taxa_to_test)){
#   print(paste0(' Model for taxon ',taxon_id))
#   current_taxon = taxa_to_test[taxon_id]
#   
#   if(current_taxon %in% result.selected.references$selected_references){
#     print('taxon in reference taxa,skipping')
#     next()
#   }
#   
#   #perform hypergeomtric sampling:
#   nom = otu_matrix[,current_taxon]
#   dnom = apply(otu_matrix[,result.selected.references$selected_references,drop = F],1,sum)
#   lambda_j = min(nom + dnom)
#   X_tilde = rhyper(nrow(otu_matrix),m = nom,n = dnom,k = lambda_j)
#   # build model - rarefaction, equivlant of DACOMP
#   data_Frame <- data.frame("W" = X_tilde,
#                            "M" = lambda_j,
#                            "X1" = sample_data(ps)$W)
#   try({
#     model = bbdml(formula = cbind(W, M - W) ~ X1,
#                   phi.formula = ~ X1,
#                   data = data_Frame);
#     model_summary = summary(model);
#     
#     P.values.rarefaction[taxon_id] = model_summary$coefficients["mu.X11","Pr(>|t|)"]  
#   })
#   
#   
#   # build model - rarefaction, equivlant of DACOMP-ratio
#   data_Frame <- data.frame("W" = nom,
#                            "M" = nom+dnom,
#                            "X1" = sample_data(ps)$W)
#   try({
#     model = bbdml(formula = cbind(W, M - W) ~ X1,
#                   phi.formula = ~ X1,
#                   data = data_Frame);
#     model_summary = summary(model);
#     
#     P.values.ratio[taxon_id] = model_summary$coefficients["mu.X11","Pr(>|t|)"]
#   })
#   
# }


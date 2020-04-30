library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(SpiecEasi)

# set working directory
setwd('/Volumes/GoogleDrive/My\ Drive/DOCTORATE/Thesis')

# load microbiome data
ASV_table_filtered <- readRDS('KORA\ DATA/Microbiome_data/dada2output/dada2output2020/seqtab2020_filtered.rds')
taxon_assign <- readRDS('KORA\ DATA/Microbiome_data/dada2output/dada2output2020/taxa2020.rds')
# load phylogenetic information
load("KORA\ DATA/Microbiome_data/dada2output/dada2output2020/phylotree2020.phy")
# load sample/matched_data
load('Microbiome\ 2020/1.\ January/environment_matching/dat_matched_PM25.RData')

sample_df <- matched_df[order(matched_df$ff4_prid),]
sample_df$W <- as.factor(sample_df$W)
samples.out <- as.character(sample_df$ff4_prid)
rownames(sample_df) <- samples.out

# create a phyloseq object
ps <- phyloseq(ASV_table_filtered,
               sample_data(sample_df),
               tax_table(taxon_assign),
               phy_tree(tGTR$tree))
ps


##############
# SPIEC-EASI #
##############

# se.gl.ps_overall <- spiec.easi(ps, method='glasso', lambda.min.ratio=1e-2,
#                                nlambda=20, pulsar.params=list(rep.num=50))

# # treated
ps_W1 <- subset_samples(ps, W == 1)
# se.gl.ps_W1 <- spiec.easi(ps_W1, method='glasso', lambda.min.ratio=1e-2,
#                           nlambda=20, pulsar.params=list(rep.num=50))
# # control
ps_W0 <- subset_samples(ps, W == 0)
# se.gl.ps_W0 <- spiec.easi(ps_W0, method='glasso', lambda.min.ratio=1e-2,
#                           nlambda=20, pulsar.params=list(rep.num=50))
# 
# ### SAVE SpiecEasi fit ###
setwd('/Volumes/GoogleDrive/My\ Drive/DOCTORATE/Thesis/Microbiome\ 2020/2.\ February')
# save(se.gl.ps_W1, file = 'SE_output_W1_PM25_ASV20.RData')
# save(se.gl.ps_W0, file = 'SE_output_W0_PM25_ASV20.RData')
# save(se.gl.ps_overall, file = 'SE_overall_PM25_ASV20.RData')

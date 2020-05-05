library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

# set working directory
setwd('/Volumes/GoogleDrive/My\ Drive/DOCTORATE/Thesis')

# load microbiome data
ASV_table <- readRDS('KORA\ DATA/Microbiome_data/dada2output/dada2output2020/seqtab2020.rds')
taxon_assign <- readRDS('KORA\ DATA/Microbiome_data/dada2output/dada2output2020/taxa2020.rds')
# load phylogenetic information
load("KORA\ DATA/Microbiome_data/dada2output/dada2output2020/phylotree2020.phy")

# create a phyloseq object
ps <- phyloseq(otu_table(ASV_table, taxa_are_rows=FALSE),
               #sample_data(sample_df),
               tax_table(taxon_assign),
               phy_tree(tGTR$tree))
ps

########################
# Prevalence Filtering #
########################

# Explore the relationship of prevalence and total read count for each feature. 
# Sometimes this reveals outliers that should probably be removed, 
# and also provides insight into the ranges of either feature that might be useful.

prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps, "Phylum"))

ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + 
  geom_point(size = 2, alpha = 0.7) + scale_x_log10() + 
  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") + 
  facet_wrap(~Phylum, nrow = 2) + theme(legend.position="none")

#  Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps)
prevalenceThreshold

# Execute prevalence  filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps)

table(tax_table(ps2)[, "Phylum"], exclude = NULL)

# saveRDS(otu_table(ps2), "/Volumes/GoogleDrive/My\ Drive/DOCTORATE/Thesis/KORA\ DATA/Microbiome_data/dada2output/dada2output2020/seqtab2020_filtered.rds")


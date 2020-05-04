library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

###############################################################################

# set working directory
setwd('/Users/alicesommer/Desktop/Bureau/DOCTORATE/data_pipeline_microbiome')

# load microbiome data
ASV_table <- readRDS('dada2output/seqtab2020.rds')
taxon_assign <- readRDS('dada2output/taxa2020.rds')
# load phylogenetic information
load("dada2output/phylotree2020.phy")

###############################################################################

# create a phyloseq object
ps <- phyloseq(otu_table(ASV_table, taxa_are_rows=FALSE),
               #sample_data(sample_df),
               tax_table(taxon_assign),
               phy_tree(tGTR$tree))
ps

par(mfrow=c(2,2))

# nr. of ASV present
locate_ASV_in_sample <- apply(otu_table(ps), 1, function(x) sum(x != 0))
hist(locate_ASV_in_sample, breaks = 30, main = "", xlab = "# ASV (sample)", xlim = c(0,400))

## sequencing depth (count statistics accross each n)
hist(sample_sums(ps), breaks = 50, main = "", xlab = "Sequencing Depth (sample)")
min(sample_sums(ps))

## count statistics accross each p
hist(taxa_sums(ps), breaks = 7000, main = "", xlab = "Total Counts (taxa)", xlim=c(1,7000))

## zero dist accross each p
zero_p <- apply(otu_table(ps), 2, function(x) sum(x == 0))
hist(zero_p, breaks = 100, main = "", xlab = "# Zeros (taxa)")

# apply function on every row of taxa table
names(taxon_assign[1,])[is.na(taxon_assign[1,])][1]

locate_NA_taxa <- apply(taxon_assign, 1, function(x) which(is.na(x))[1])
table(locate_NA_taxa)
sum(table(locate_NA_taxa))

locate_NA_taxa <- apply(taxon_assign, 1, function(x) sum(is.na(x)))
table(locate_NA_taxa)
sum(table(locate_NA_taxa))


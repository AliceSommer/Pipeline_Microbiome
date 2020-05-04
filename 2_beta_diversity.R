
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(plyr); packageVersion("plyr")
library(pldist)

###############################################################################

# set working directory
setwd('/Users/alicesommer/Desktop/Bureau/DOCTORATE/data_pipeline_microbiome')

# load microbiome data
ASV_table <- readRDS('dada2output/seqtab2020.rds')
taxon_assign <- readRDS('dada2output/taxa2020.rds')
# load phylogenetic information
load("dada2output/phylotree2020.phy")

# load sample/matched_data
load('dat_matched_PM25.RData')

###############################################################################

sample_df <- matched_df[order(matched_df$ff4_prid),]
sample_df$W <- as.factor(sample_df$W)
samples.out <- as.character(sample_df$ff4_prid)
rownames(sample_df) <- samples.out

# create a phyloseq object
ps <- phyloseq(otu_table(ASV_table, taxa_are_rows=FALSE),
               sample_data(sample_df),
               tax_table(taxon_assign),
               phy_tree(tGTR$tree))
ps

# Input: Notice that row names are sample IDs 
paired.otus <- as(otu_table(ps), "matrix") 
paired.otus[1:4,1:4]

paired.meta <- sample_data(ps)[,c("pair_nb","ff4_prid","W")]
colnames(paired.meta) <- c("subjID", "sampID", "time")
paired.meta[1:4,]

# Transformation function 
otu.data <- data_prep(paired.otus, paired.meta, paired = TRUE, pseudoct = NULL)
otu.data$otu.props[1:3,1:3]  # OTU proportions 
otu.data$otu.clr[1:3,1:3]    # CLR-transformed proportions
res <- pltransform(otu.data, paired = TRUE, norm = TRUE)


# Binary transformation 
# 0.5 indicates OTU was present at Time 2, absent at Time 1
# -0.5 indicates OTU was present at Time 1, absent at Time 2 
# Row names are now subject IDs 
res$dat.binary[1:3,1:3]

# Quantitative transformation (see details in later sections)
round(res$dat.quant.prop[1:3,1:3], 2)
round(res$dat.quant.clr[1:3,1:3], 2)
# Average proportion per OTU per subject 
round(res$avg.prop[1:3,1:3], 2)
# This was a paired transformation 
res$type

# LUniFrac
D.unifrac <- LUniFrac(otu.tab = paired.otus, metadata = paired.meta, tree = tGTR$tree,
                           gam = c(0, 0.5, 1), paired = TRUE, check.input = TRUE)
head(D.unifrac[, , "d_1"]) # gamma = 1 (quantitative paired transformation)
head(D.unifrac[, , "d_UW"])  # unweighted LUniFrac (qualitative/binary paired transf.)]

fit_uni <- cmdscale(D.unifrac[, , "d_UW"],eig=TRUE, k=2) # k is the number of dim
fit_uni

# plot solution
x_uni <- fit_uni$points[,1]
y_uni <- fit_uni$points[,2]
plot(x_uni, y_uni, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS - Unifrac (unweighted) - pldist", col = sample_data(ps)$W)

############################
##### Other distances ######
############################

gower_dist <- pldist(paired.otus, paired.meta, paired = TRUE, binary = FALSE, 
                     method = "gower", clr = TRUE)$D
head(gower_dist)

fit <- cmdscale(gower_dist,eig=TRUE, k=2) # k is the number of dim
fit

# plot solution
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS - Gower", col = sample_data(ps)$W)

bray_dist <- pldist(paired.otus, paired.meta, paired = TRUE, binary = FALSE, 
                     method = "bray", clr = TRUE)$D

fit_b <- cmdscale(bray_dist,eig=TRUE, k=2) # k is the number of dim
fit_b

# plot solution
x_b <- fit_b$points[,1]
y_b <- fit_b$points[,2]
plot(x_b, y_b, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS - Bray", col = sample_data(ps)$W)


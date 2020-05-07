
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(plyr); packageVersion("plyr")
library(pldist)
library(MiRKAT)
library(compositions)

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

# locate the species that are totally absent in the matched data
empty_species <- colSums(otu_table(ps))
length(which(empty_species == 0))

ps_prune <- prune_taxa(empty_species != 0, ps)

###########################################
##### Calculate distances with pldist #####
###########################################

# Input: Notice that row names are sample IDs 
paired.otus <- as(otu_table(ps_prune), "matrix") 
paired.otus[1:4,1:4]

paired.meta <- sample_data(ps_prune)[,c("pair_nb","ff4_prid","W")]
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
# head(D.unifrac[, , "d_1"]) # gamma = 1 (quantitative paired transformation)
# head(D.unifrac[, , "d_UW"])  # unweighted LUniFrac (qualitative/binary paired transf.)]

#### attention !! #### cannot test MiRKAT with W as trait/exposure because W used for the pairs 
## dimension of D.unifrac is n_pair x n_pair

# try to check if influence on diabetes
K.unweighted <- D2K(D.unifrac[, , "d_UW"])

id_pair_var <- data.frame(pair_nb = rownames(D.unifrac[, , "d_UW"]))
head(id_pair_var)

colnames(sample_data(ps_prune))[grep("glu",colnames(sample_data(ps)))]
table(sample_data(ps_prune)$u3tdiabet)

dat_dia_pair <- data.frame(sample_data(ps_prune)[sample_data(ps_prune)$W == 1,c("pair_nb","u3tdiabet","u3csex")])

id_pair_var <- merge(id_pair_var, dat_dia_pair, sort = FALSE,
                     by.x = "pair_nb", by.y = "pair_nb",all.x = TRUE)

head(id_pair_var)
dim(K.unweighted)
# testing using a single Kernel
MiRKAT(y = id_pair_var$u3tdiabet, X = NULL, Ks = K.unweighted, out_type = "D", 
       method = "davies", returnKRV = TRUE, returnR2 = TRUE)

#### maybe interesting for microbiome -> AP -> diabetes ####

#####################################
##### Global hypothesis testing #####
#####################################

## MiRKAT
unifracs <- UniFrac(ps_prune, weighted=FALSE, parallel=FALSE, fast=TRUE)
unifracs_mat <- as.matrix(unifracs)
# transform distance matrix to kernel
K.unweighted_Uuni <- D2K(unifracs_mat)

head(sample_data(ps_prune)$W)
outcome <- as.numeric(sample_data(ps_prune)$W == 1)
head(outcome)

## testing using a single Kernel
MiRKAT(y = outcome, X = NULL, Ks = K.unweighted_Uuni, out_type = "D", 
       method = "davies", returnKRV = TRUE, returnR2 = TRUE)

## omnibus test if multiple distance matrices ("Optimal MiRKAT")
ps_clr <- transform_sample_counts(ps_prune, function(x){clr(x)})
ps_comp <- transform_sample_counts(ps, function(x){x/sum(x)})

# the available distance methods coded in distance
dist_methods <- unlist(distanceMethodList)
print(dist_methods)

# Euclidean
euclidean_dist <- distance(ps_clr, method="euclidean")
K.euclidean_dist <- D2K(as.matrix(euclidean_dist))

# Bray
bray_dist <- distance(ps_comp, method="bray") # RA
K.bray_dist<- D2K(as.matrix(bray_dist))

# Jaccard 
jaccard_dist <- distance(ps_comp, method="jaccard") # RA
K.jaccard_dist <- D2K(as.matrix(jaccard_dist))

# Kulczynski
kulczynski_dist <- distance(ps_comp, method="kulczynski") # RA
K.kulczynski_dist <- D2K(as.matrix(kulczynski_dist))

# Gower
gower_dist <- distance(ps_comp, method="gower") # RA
K.gower_dist <- D2K(as.matrix(gower_dist))

## testing using a several Kernels
Ks = list(K.unweighted_Uuni = K.unweighted_Uuni, K.euclidean_dist = K.euclidean_dist, 
          K.bray_dist = K.bray_dist, K.jaccard_dist = K.jaccard_dist, 
          K.kulczynski_dist = K.kulczynski_dist, K.gower_dist = K.gower_dist)
MiRKAT(y = outcome, Ks = Ks, X = NULL, out_type = "D", 
       method = "davies", returnKRV = TRUE, returnR2 = TRUE)


###########################################
##### Low-dimensional representation ######
###########################################
fit_eu <- cmdscale(euclidean_dist,eig=TRUE, k=2) # k is the number of dim
# fit_eu

# plot solution
x_eu <- fit_eu$points[,1]
y_eu <- fit_eu$points[,2]
plot(x_uni, y_uni, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS - Unifrac (unweighted) - pldist", col = sample_data(ps)$W)

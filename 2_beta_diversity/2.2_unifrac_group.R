library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(glmnet)
library(compositions)
library(pheatmap)

###############################################################################

# set working directory
setwd('/Users/alicesommer/Desktop/Bureau/DOCTORATE')

# load microbiome data
ASV_table <- readRDS('data_pipeline_microbiome/dada2output/seqtab2020.rds')
taxon_assign <- readRDS('data_pipeline_microbiome/dada2output/taxa2020.rds')
# load phylogenetic information
load("data_pipeline_microbiome/dada2output/phylotree2020.RData")

# load sample/matched_data
# load('data_pipeline_microbiome/dat_matched_PM25_bis.RData')
load('data_pipeline_microbiome/dat_matched_smoke_bis.RData')

# load W matrix for randomization test
# load("data_pipeline_microbiome/W_paired_PM25.Rdata")
load("data_pipeline_microbiome/W_paired_smoke_bis.Rdata")

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
# 1. have a vector of nr. of observed samples per taxa
vec_taxa <- apply(otu_table(ps), 2, function(x) sum(x > 0, na.rm = TRUE))
length(which(vec_taxa == 0))
# 2. how many tax are obs. in at least x% of samples
perc <- 0.05 # x%
# samp_perc <- trunc(dim(sample_df)[1]*perc)
samp_perc <- 0
length(which(vec_taxa > samp_perc))

ps_prune <- prune_taxa(vec_taxa > samp_perc, ps)

###############################################################################
iDist <- distance(ps_prune, method='unifrac')
iMDS  <- ordinate(ps_prune, "MDS", distance=iDist)

sample_data(ps_prune)$seq_depth <- apply(otu_table(ps_prune), 1, function(x) sum(x, na.rm = TRUE))

plot_ordination(ps_prune, iMDS, color="W") 

sample_data(ps_prune)$uni_group <- as.factor(as.numeric(iMDS$vectors[,2] > -0.095))
table(sample_data(ps_prune)$uni_group)

###############################################################################
#### GLM ####

grep('bmi', colnames(sample_data(ps_prune)))

Y <- sample_data(ps_prune)$uni_group
X <- data.matrix(sample_data(ps_prune)[,c(10:138,148,150:152,2283)])
X <- X[,-which(apply(X, 2, function(x) sum(is.na(x))) > 0)]

fit <- cv.glmnet(X, Y, family = "binomial", type.measure = "class")

plot(fit)

fit$lambda.min

coef(fit, s = "lambda.min")

# u3tbmi
# u3tmlg04bd : Urologische Spasmolytika --- (ATC = G04BD)
# u3tmlm03b : Muskelrelaxantien, zentral wirkend Mittel --- (ATC = M03B)
# u3tmlr03 : Sympatomimetika: respiratorisch --- (ATC = R03)
# u3tmbbl : Einnahme von Beta-Blockern
# u3tmglyk : Einnahme von Herzglykosiden
# u3tphact1 : Physical Activity in Kategorien (1: regelmäßig mehr als 2 Std. in der Woche)

apply((sample_data(ps_prune)[, c('u3tmlg04bd','u3tmlm03b','u3tmlr03','u3tmbbl','u3tmglyk')]), 2, sum)

sample_data(ps_prune)[,c('u3tmbbl', 'u3tphact1','W',"year")]

table(sample_data(ps_prune)$u3tmbbl, sample_data(ps_prune)$W)

table(sample_data(ps_prune)$u3tphact1, sample_data(ps_prune)$W)

table(sample_data(ps_prune)$year, sample_data(ps_prune)$W)

predict(fit, newx = X[sample(1:530,100),], s = "lambda.min", type = "class")

###############################################################################
#### Heatmap ####

# clr transform data
level <- "Genus"
ps_agg <- tax_glom(ps_prune, taxrank = level, NArm = FALSE)
ps_work <- transform_sample_counts(ps_agg, function(x) {clr(x + .5)})

X_heat <- as(otu_table(ps_work), "matrix") # unname to remove sequences as "title"
X_heat[otu_table(ps_agg) == 0] <- 0 # replace the pseudo-counted zeros to zeros again

X_split <- X_heat[order(sample_data(ps_work)$uni_group),]

data.prop <- t(X_split)
rownames(data.prop) <- tax_table(ps_work)[,level]

# Data frame with column annotations.
new_W <- sample_data(ps_work)$uni_group[order(sample_data(ps_work)$uni_group)]
mat_col <- data.frame(uni_group = as.factor(new_W))
rownames(mat_col) <- colnames(data.prop)

# List with colors for each annotation.
mat_colors <- list(group = brewer.pal(3, "Set1"))
names(mat_colors$group) <- unique(sample_df$uni_group)

# pheatmap(
#   mat               = data.prop,
#   # color             = inferno(length(mat_breaks) - 1),
#   # border_color      = NA,
#   gaps_col = sum(new_W == 0),
#   cluster_rows = TRUE,
#   cluster_cols = FALSE,
#   show_colnames     = FALSE,
#   show_rownames     = TRUE,
#   annotation_col    = mat_col,
#   annotation_colors = list(uni_group = c('0' = "#E41A1C", '1' = "#4DAF4A")),
#   drop_levels       = FALSE,
#   fontsize_row      = 10,
#   cellwidth = 1.6
#   # main              = ""
#   # cutree_cols = 3
# )

mat_col_clust <- data.frame(uni_group = as.factor(sample_data(ps_work)$uni_group),
                            bmi = sample_data(ps_work)$u3tbmi,
                            phys_act = as.factor(sample_data(ps_work)$u3tphact1), 
                            beta_block = as.factor(sample_data(ps_work)$u3tmbbl),
                            W = as.factor(sample_data(ps_work)$u3tmbbl))

rownames(mat_col_clust) <- rownames(sample_data(ps_work))

pheatmap(
  mat               = data.prop,
  border_color      = NA,
  show_colnames     = FALSE,
  show_rownames     = TRUE,
  annotation_col    = mat_col_clust,
  annotation_colors = list(W = c('0' = "#E41A1C", '1' = "#4DAF4A")),
  drop_levels       = TRUE,
  fontsize          = 14,
  # cellwidth = 1.6,
  main              = ""
)

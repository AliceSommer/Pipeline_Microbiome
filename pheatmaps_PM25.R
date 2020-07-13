library(phyloseq); packageVersion("phyloseq")
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(compositions)

# set working directory
setwd('/Users/alicesommer/Desktop/Bureau/DOCTORATE')

# load microbiome data
ASV_table <- readRDS('data_pipeline_microbiome/dada2output/seqtab2020.rds')
taxon_assign <- readRDS('data_pipeline_microbiome/dada2output/taxa2020.rds')

# load sample/matched_data
load('data_pipeline_microbiome/dat_matched_PM25_bis.RData')

sample_df <- matched_df[order(matched_df$ff4_prid),]
sample_df$W <- as.factor(sample_df$W)
samples.out <- as.character(sample_df$ff4_prid)
rownames(sample_df) <- samples.out

# create a phyloseq object
ps <- phyloseq(otu_table(ASV_table, taxa_are_rows=FALSE),
               sample_data(sample_df),
               tax_table(taxon_assign))
ps

rank_names(ps)

vec_taxa <- apply(otu_table(ps), 2, function(x) sum(x > 0, na.rm = TRUE))
ps_prune <- prune_taxa(vec_taxa >= 9 , ps)

level = "Species"

## agglomerate ##
ps_agg <- tax_glom(ps_prune, taxrank = level, NArm = FALSE)
# ps_work <- transform_sample_counts(ps_agg, function(x) {x/sum(x)})
ps_work <- transform_sample_counts(ps_agg, function(x) {clr(x + .5)})

############
# PHEATMAP #
############

X <- as(otu_table(ps_work), "matrix") # unname to remove sequences as "title"

data.prop <- t(X)
# data.prop[1:3, 1:3]
rownames(data.prop) <- tax_table(ps_work)[,level]

# Data frame with column annotations.
mat_col <- data.frame(W = as.factor(sample_df$W))
rownames(mat_col) <- colnames(data.prop)

# List with colors for each annotation.
mat_colors <- list(group = brewer.pal(3, "Set1"))
names(mat_colors$group) <- unique(sample_df$W)

pheatmap(
  mat               = data.prop,
  # color             = inferno(length(mat_breaks) - 1),
  # border_color      = NA,
  show_colnames     = FALSE,
  show_rownames     = TRUE,
  annotation_col    = mat_col,
  annotation_colors = list(W = c('0' = "#E41A1C", '1' = "#4DAF4A")),
  drop_levels       = FALSE,
  fontsize_row      = 6,
  cellwidth = 1.3
  # main              = ""
  # cutree_cols = 3
)








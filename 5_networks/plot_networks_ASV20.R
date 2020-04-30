library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(SpiecEasi)
library(RColorBrewer)
library(ForceAtlas2)
library(igraph)


library(EnvStats)
library(gridExtra)
library(Matrix)
library(ForceAtlas2)


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

# # treated
ps_W1 <- subset_samples(ps, W == 1)

# # control
ps_W0 <- subset_samples(ps, W == 0)

load('Microbiome\ 2020/2.\ February/SE_output_W1_PM25_ASV20.RData')
load('Microbiome\ 2020/2.\ February/SE_output_W0_PM25_ASV20.RData')
load('Microbiome\ 2020/2.\ February/SE_overall_PM25_ASV20.RData')

rm(otu_table, samples.out, taxa_final, tGTR, sample_df, matched_df)

########## OVERALL ##########
ig.gl_overall <- adj2igraph(getRefit(se.gl.ps_overall),
                            vertex.attr=list(Phylum = tax_table(ps)[, 'Phylum'], 
                                             Class = tax_table(ps)[, 'Class'],
                                             Specie = tax_table(ps)[, 'Species']))

# color
pal <- brewer.pal(length(unique(vertex_attr(ig.gl_overall)$Phylum)), "Paired")
# add color as vertex attribute
ig.gl_overall <- set_vertex_attr(ig.gl_overall, "color_palette", index = V(ig.gl_overall), pal[as.factor(vertex_attr(ig.gl_overall, "Phylum"))])
# calculate weights
secor <- cov2cor(getOptCov(se.gl.ps_overall))
elist.gl <- summary(triu(secor*getRefit(se.gl.ps_overall), k=1))
elist.gl <- elist.gl[order(elist.gl$i),] 
# add weights as edges attribute
ig.gl_overall <- set_edge_attr(ig.gl_overall, "weight", index = E(ig.gl_overall), elist.gl$x)

# set.seed(87)
# gl.coord <- layout.forceatlas2(ig.gl_overall, iterations=2000, plotstep=100)
# V(ig.gl_overall)$layout_x <- gl.coord[,1]
# V(ig.gl_overall)$layout_y <- gl.coord[,2]

# color of weights
E(ig.gl_overall)$color <- "blue"
E(ig.gl_overall)[weight<0]$color <- "pink"

par(mfrow=c(1,1))
plot(ig.gl_overall, layout=gl.coord, 
     vertex.size=4, edge.width = 5*abs(E(ig.gl_overall)$weight),
     vertex.label=NA,
     vertex.color = vertex_attr(ig.gl_overall, 'color_palette'))


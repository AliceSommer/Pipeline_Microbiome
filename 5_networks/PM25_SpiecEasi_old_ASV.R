library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(EnvStats)
library(gridExtra)
library(dplyr)
library(plyr)
library(SpiecEasi)
library(igraph)
library(Matrix)
library(RColorBrewer)
library(ForceAtlas2)
library(orca)
library(corrplot)
library(MASS)

# set working directory
setwd('/Volumes/GoogleDrive/My\ Drive/DOCTORATE/Thesis')

# load Microbiome data
otu_table <- readRDS('/Volumes/GoogleDrive/My\ Drive/DOCTORATE/Thesis/KORA\ DATA/Microbiome_data/dada2output/otu_table_filtered_D.rds')
taxa_final <- readRDS('/Volumes/GoogleDrive/My\ Drive/DOCTORATE/Thesis/KORA\ DATA/Microbiome_data/dada2output/taxa_final.rds')
load("/Volumes/GoogleDrive/My Drive/DOCTORATE/Thesis/KORA DATA/Microbiome_data/dada2output/tGTR.phy")

# load matched data Low vs. High AP
load('Microbiome\ 2020/1.\ January/environment_matching/dat_matched_PM25.RData')

sample_df <- matched_df[order(matched_df$ff4_prid),]
sample_df$W <- as.factor(sample_df$W)
samples.out <- as.character(sample_df$ff4_prid)
rownames(sample_df) <- samples.out

# create a phyloseq object
ps <- phyloseq(otu_table, ## BIG OTU TABLE
               sample_data(sample_df),
               tax_table(taxa_final),
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
# setwd('/Volumes/GoogleDrive/My\ Drive/DOCTORATE/Thesis/Microbiome\ 2020/1.\ January')
# save(se.gl.ps_W1, file = 'SE_output_W1_PM25.RData')
# save(se.gl.ps_W0, file = 'SE_output_W0_PM25.RData')
#save(se.gl.ps_overall, file = 'SE_output_overall_PM25.RData')

locate_NA_taxa_filt <- apply(tax_table(ps), 1, function(x) sum(is.na(x)))
sum(table(locate_NA_taxa_filt))
table(locate_NA_taxa_filt)


#############
# Diversity #
#############

ps_rich <- data.frame(estimate_richness(ps, measures=c("Observed", "InvSimpson", "Shannon", "Chao1", "Simpson")))
head(ps_rich)
sample_data(ps)[,c("Observed_Rich_phyl", "Chao1_phyl", "se.chao1_phyl", "Shannon_phyl", "Simpson_phyl", "InvSimpson_phyl")] <- ps_rich

g_PM <- ggplot(sample_data(ps), aes(color = factor(W), y = Chao1_phyl)) +
  geom_boxplot(alpha = .5) + ylab('Chao 1 diversity measure') +
  scale_x_discrete(name = "") +
  scale_colour_manual(values = c("gray","blue4"), limits=c("1","0"), name ="Long-term PM2.5", labels = c("Low (<= 11)","High (>= 12)")) +
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) +
  guides(color=guide_legend(nrow=2,byrow=TRUE))

################
# Plot network #
################

load('/Volumes/GoogleDrive/My\ Drive/DOCTORATE/Thesis/Microbiome\ 2020/1.\ January/SE_output_W1_PM25.RData')
load('/Volumes/GoogleDrive/My\ Drive/DOCTORATE/Thesis/Microbiome\ 2020/1.\ January/SE_output_W0_PM25.RData')
load('/Volumes/GoogleDrive/My\ Drive/DOCTORATE/Thesis/Microbiome\ 2020/1.\ January/SE_output_overall_PM25.RData')

rm(otu_table, samples.out, taxa_final, tGTR, sample_df, matched_df)

########## OVERALL ##########
ig.gl_overall <- adj2igraph(getRefit(se.gl.ps_overall),
                            vertex.attr=list(Phylum = tax_table(ps)[, 'Phylum'], 
                                             Class = tax_table(ps)[, 'Class'],
                                             Specie = tax_table(ps)[, 'Species']))

# color
pal <- brewer.pal(length(unique(vertex_attr(ig.gl_overall)$Phylum)), "Spectral")
# add color as vertex attribute
ig.gl_overall <- set_vertex_attr(ig.gl_overall, "color_palette", index = V(ig.gl_overall), pal[as.factor(vertex_attr(ig.gl_overall, "Phylum"))])
# calculate weights
secor <- cov2cor(getOptCov(se.gl.ps_overall))
elist.gl <- summary(triu(secor*getRefit(se.gl.ps_overall), k=1))
elist.gl <- elist.gl[order(elist.gl$i),] 
# add weights as edges attribute
ig.gl_overall <- set_edge_attr(ig.gl_overall, "weight", index = E(ig.gl_overall), elist.gl$x)

set.seed(87)
gl.coord <- layout.forceatlas2(ig.gl_overall, iterations=2000, plotstep=100)
V(ig.gl_overall)$layout_x <- gl.coord[,1]
V(ig.gl_overall)$layout_y <- gl.coord[,2]

# color of weights
E(ig.gl_overall)$color <- "blue"
E(ig.gl_overall)[weight<0]$color <- "pink"

par(mfrow=c(1,1))
plot(ig.gl_overall, layout=gl.coord, 
     vertex.size=3, edge.width = 5*abs(E(ig.gl_overall)$weight),
     vertex.label=NA,
     vertex.color = vertex_attr(ig.gl_overall, 'color_palette'))

########## TREATED ##########
## Create igraph objects
ig.gl <- adj2igraph(getRefit(se.gl.ps_W1), 
                    vertex.attr=list(ASV = rownames(tax_table(ps_W0)),
                                     Phylum = tax_table(ps_W1)[, 'Phylum'], 
                                     Class = tax_table(ps_W1)[, 'Class'],
                                     Specie = tax_table(ps_W1)[, 'Species']))

# color
pal <- brewer.pal(length(unique(vertex_attr(ig.gl)$Phylum)), "Spectral")
# add color as vertex attribute
ig.gl <- set_vertex_attr(ig.gl, "color_palette", index = V(ig.gl), pal[as.factor(vertex_attr(ig.gl, "Phylum"))])
# calculate weights
secor <- cov2cor(getOptCov(se.gl.ps_W1))
elist.gl <- summary(triu(secor*getRefit(se.gl.ps_W1), k=1))
elist.gl <- elist.gl[order(elist.gl$i),] 
# add weights as edges attribute
ig.gl <- set_edge_attr(ig.gl, "weight", index = E(ig.gl), elist.gl$x)

# layout of the vertices
# set.seed(87)
# gl.coord <- layout.forceatlas2(ig.gl, iterations=2000, plotstep=100)
# V(ig.gl)$layout_x <- gl.coord[,1]
# V(ig.gl)$layout_y <- gl.coord[,2]

# layout of the vertices
V(ig.gl)$layout_x <- gl.coord[,1]
V(ig.gl)$layout_y <- gl.coord[,2]

# color of weights
E(ig.gl)$color <- "blue"
E(ig.gl)[weight<0]$color <- "pink"

plot(ig.gl, layout=gl.coord, 
     vertex.size=3, edge.width = 5*abs(E(ig.gl)$weight),
     vertex.label=NA,
     vertex.color = vertex_attr(ig.gl, 'color_palette'))

# # edge weights plot
# hist(elist.gl[,3], main='', xlab='edge weights')
# 
# # degree statistic
# dd.gl <- degree.distribution(ig.gl)
# # plot
# plot(0:(length(dd.gl)-1), dd.gl, ylim=c(0,.35), type='b',
#      ylab="Frequency", xlab="Degree", main="Degree Distributions")

########## CONTROL ##########

## Create igraph objects
ig.gl_C <- adj2igraph(getRefit(se.gl.ps_W0), 
                      vertex.attr=list(ASV = rownames(tax_table(ps_W1)),
                                       Phylum = tax_table(ps_W1)[, 'Phylum'], 
                                       Class = tax_table(ps_W1)[, 'Class'],
                                       Specie = tax_table(ps_W1)[, 'Species']))

# add color as vertex attribute
ig.gl_C <- set_vertex_attr(ig.gl_C, "color_palette", index = V(ig.gl_C), pal[as.factor(vertex_attr(ig.gl_C, "Phylum"))])
# calculate weights
secor_C <- cov2cor(getOptCov(se.gl.ps_W0))
elist.gl_C <- summary(triu(secor_C*getRefit(se.gl.ps_W0), k=1))
elist.gl_C <- elist.gl[order(elist.gl_C$i),] 
# add weights as edges attribute
ig.gl_C <- set_edge_attr(ig.gl_C, "weight", index = E(ig.gl_C), elist.gl_C$x)

# layout of the vertices
V(ig.gl_C)$layout_x <- gl.coord[,1]
V(ig.gl_C)$layout_y <- gl.coord[,2]

# color of weights
E(ig.gl_C)$color <- "blue"
E(ig.gl_C)[weight<0]$color <- "pink"

par(mfrow=c(1,2))
plot(ig.gl, layout=cbind(V(ig.gl)$layout_x, V(ig.gl)$layout_y), 
     vertex.size=3, edge.width = 5*abs(E(ig.gl)$weight), 
     vertex.label=NA, 
     vertex.color = vertex_attr(ig.gl, 'color_palette'),
     main="treated, low PM")
legend(.75,-.6, bty = "n", cex = .85,
       legend=levels(as.factor(vertex_attr(ig.gl_C)$Phylum)),
       fill=pal, border=NA)
plot(ig.gl_C, layout=cbind(V(ig.gl_C)$layout_x, V(ig.gl_C)$layout_y), 
     vertex.size=3, edge.width = 5*abs(E(ig.gl_C)$weight),
     vertex.label=NA, 
     vertex.color = vertex_attr(ig.gl, 'color_palette'),
     main="control, high PM")

##############
# Difference #    
##############

### intersection ###
int <- graph.intersection(ig.gl,ig.gl_C, keep.all.vertices = FALSE) ### the core graph 
length(E(int)) 
### weights??????
E(int)$weight_avg = (E(int)$weight_1 + E(int)$weight_2)/2

# color of weights
E(int)$color <- "blue"
E(int)[weight_avg<0]$color <- "pink"

# plot the core graph
int_no_empty_vertice = delete.vertices(int,degree(int)<1)

par(mfrow=c(1,1))
plot(int_no_empty_vertice, 
     layout = cbind(V(int_no_empty_vertice)$layout_x_1, V(int_no_empty_vertice)$layout_y_1),
     vertex.size=3, 
     edge.width = 7*abs(E(int_no_empty_vertice)$weight_avg),
     vertex.label=NA, 
     vertex.color = vertex_attr(int_no_empty_vertice, 'color_palette_1'),
     main="core graph")
legend(.75,-.6, bty = "n", cex = .85,
       legend=levels(as.factor(vertex_attr(ig.gl_C)$Phylum)),
       fill=pal, border=NA)


#######################
# Only positive edges #
#######################

# keep original vertex ID
V(ig.gl)$vertex_ID <- V(ig.gl)
V(ig.gl_C)$vertex_ID <- V(ig.gl_C)

positive_g <- subgraph.edges(ig.gl, eids = which(E(ig.gl)$weight > 0), delete.vertices = TRUE)
positive_g_C <- subgraph.edges(ig.gl_C, eids = which(E(ig.gl_C)$weight > 0), delete.vertices = TRUE)

##########################
# Fully connected graph #
#########################

dg_treated <- decompose.graph(positive_g)
ig.gl_connect <- dg_treated[[1]]

dg_control <- decompose.graph(positive_g_C)
ig.gl_C_connect <- dg_control[[1]]

#######################
# Community structure #
#######################

##### FAST GREEDY 2004 (Newman)
## treated
fc_treated <- cluster_fast_greedy(positive_g, weights = NULL)
membership(fc_treated)
sizes(fc_treated)
modularity(fc_treated)
length(fc_treated)

V(positive_g)$community <- membership(fc_treated)

## control
fc_control <- cluster_fast_greedy(positive_g_C, weights = NULL)
membership(fc_control)
sizes(fc_control)
modularity(fc_control)
length(fc_control)

V(positive_g_C)$community <- membership(fc_control)

##########################
## Compare memberships ##
#########################

length(V(positive_g)$vertex_ID)
length(V(positive_g_C)$vertex_ID)

treated_com <- data.frame(id = V(positive_g)$vertex_ID, com_t = V(positive_g)$community)
head(treated_com)
dim(treated_com)

control_com <- data.frame(id = V(positive_g_C)$vertex_ID, com_c = V(positive_g_C)$community)
head(control_com)
dim(control_com)

merged_com <- merge(treated_com, control_com, by = "id", all = T)
dim(merged_com)
head(merged_com)

# every edge is a microbe
inc_com <- merge(treated_com, control_com, by = "id")
inc_com[,2] <- paste0(inc_com[,2], "_T")
inc_com[,3] <- paste0(inc_com[,3], "_C")

g <- graph.data.frame(inc_com[,c(2,3)], directed = F)
V(g)$type <- V(g)$name %in% inc_com[,2] #the first column of edges is TRUE type

V(g)$color <- V(g)$type
V(g)$color=gsub("FALSE","orange",V(g)$color)
V(g)$color=gsub("TRUE","lightblue",V(g)$color)

layout_bi <- layout_with_sugiyama(g, hgap = 30)

plot(g, vertex.size = 15, layout=layout_bi$layout)

# every edge is a bunch of microbes
inc_com_sum <- ddply(data.frame(inc_com), .(com_t, com_c), summarise, sum = length(id))

g_sum <- graph.data.frame(inc_com_sum[,c(1,2)], directed = F)
V(g_sum)$type <- V(g_sum)$name %in% inc_com_sum[,1] #the first column of edges is TRUE type

# node size represents the number of OTUs in each module
sizes(fc_treated)
sizes(fc_control)
# add community name (without T or C)
V(g_sum)$com_nr <- as.numeric(unlist(strsplit(V(g_sum)$name, "_"))[c(TRUE,FALSE)])
V(g_sum)$com_size[V(g_sum)$type == TRUE] <- sizes(fc_treated)[V(g_sum)$com_nr[V(g_sum)$type == TRUE]] # treated sizes
V(g_sum)$com_size[V(g_sum)$type == FALSE] <- sizes(fc_control)[V(g_sum)$com_nr[V(g_sum)$type == FALSE]] # control sizes

V(g_sum)$color <- V(g_sum)$type
V(g_sum)$color=gsub("FALSE","orange",V(g_sum)$color)
V(g_sum)$color=gsub("TRUE","lightblue",V(g_sum)$color)
E(g_sum)$weight <- as.numeric(inc_com_sum[,3])

layout_bi <- layout_with_sugiyama(g_sum, hgap = 30)
layout_bi_rows <- layout_as_bipartite(g_sum, hgap = 50, vgap = 60, maxiter = 1000)

# order by vertex size
par(mfrow=c(1,1))
plot(g_sum, vertex.size = 15, edge.width=E(g_sum)$weight, layout=layout_bi$layout)
plot(g_sum, vertex.size = V(g_sum)$com_size/3, edge.width=E(g_sum)$weight, 
     layout=layout_bi_rows[,2:1], vertex.label.dist=5, vertex.label.degree = pi*V(g_sum)$type, 
     margin = -.3, asp = 2, edge.arrow.size=.7)

##############################################
## Dig into the treated module 2 separation ##
##############################################

# color
sizes(fc_treated) # keep only above 4
pal_mod <- brewer.pal(7, "Dark2")

V(positive_g)$color_community[!V(positive_g)$community %in% c(1:7)] <- "grey"
com_big <- V(positive_g)$community[V(positive_g)$community %in% c(1:7)]
V(positive_g)$color_community[V(positive_g)$community %in% c(1:7)] <- pal_mod[as.factor(com_big)]

# color control
sizes(fc_control) # keep only above 3
pal_mod_C <- brewer.pal(10, "Set3")

V(positive_g_C)$color_community[!V(positive_g_C)$community %in% c(1:10)] <- "grey"
com_big_C <- V(positive_g_C)$community[V(positive_g_C)$community %in% c(1:10)]
V(positive_g_C)$color_community[V(positive_g_C)$community %in% c(1:10)] <- pal_mod_C[as.factor(com_big_C)]

# grey weigth color for visibility
E(positive_g)$color <- "grey"
E(positive_g_C)$color <- "grey"

par(mfrow=c(1,2))
plot(positive_g, vertex.size=3, vertex.label=NA, layout = cbind(V(positive_g)$layout_x, V(positive_g)$layout_y), 
     edge.width = 2*abs(E(positive_g)$weight), vertex.color = vertex_attr(positive_g, 'color_community'),
     main="treated, low PM")
legend(1,-.4, bty = "n", cex = .85,
       legend=levels(as.factor(com_big)),
       fill=pal_mod, border=NA)

plot(positive_g_C, vertex.size=3, vertex.label=NA, layout = cbind(V(positive_g_C)$layout_x, V(positive_g_C)$layout_y), 
     edge.width = 2*abs(E(positive_g_C)$weight), vertex.color = vertex_attr(positive_g_C, 'color_community'),
     main="control, high PM")
legend(1,-.4, bty = "n", cex = .85,
       legend=levels(as.factor(com_big_C)),
       fill=pal_mod_C, border=NA)

g2 = induced_subgraph(graph=positive_g, vids = which(V(positive_g)$community == 2))

plot(g2, layout = cbind(V(g2)$layout_x, V(g2)$layout_y), vertex.size=3, vertex.label=NA, 
     edge.width = 3*abs(E(g2)$weight), vertex.color = vertex_attr(g2, 'color_community'))

V(positive_g_C)$community_control <- NA
V(positive_g_C)$community_control[V(positive_g_C)$vertex_ID %in% V(positive_g)$vertex_ID] <- V(positive_g)$community[V(positive_g)$vertex_ID %in% V(positive_g_C)$vertex_ID]

g2_C = induced_subgraph(graph=positive_g_C, vids = which(V(positive_g_C)$community_control == 2))

plot(g2_C, layout = cbind(V(g2_C)$layout_x, V(g2_C)$layout_y), vertex.size=3, vertex.label=NA, 
     edge.width = 3*abs(E(g2_C)$weight), vertex.color = vertex_attr(g2_C, 'color_community'))

library(ggtree)
load('/Volumes/GoogleDrive/My\ Drive/DOCTORATE/Thesis/Microbiome\ -\ DEC19/small_tree.RData')

tree_small <- tree_small + theme(legend.position="bottom")

tree_small + geom_tiplab(aes(label = Class, subset=(label %in% V(g2)$ASV), angle=angle), size = 2.3)



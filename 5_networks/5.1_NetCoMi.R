# devtools::install_github("stefpeschel/NetCoMi", dependencies = TRUE, 
#                          repos = c("https://cloud.r-project.org/",
#                                    BiocManager::repositories()))

library(NetCoMi)
library(metagMisc)
library(phyloseq); packageVersion("phyloseq")
library(igraph)
library(plyr)

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

# relabel the ASV table to have numbers consistent with taxa table
taxa_names(ps) <- paste0("Seq", seq(ntaxa(ps)))

# filter taxa at 10% prevalence
vec_taxa <- apply(otu_table(ps), 2, function(x) sum(x > 0, na.rm = TRUE))
cut_sample <- trunc(dim(sample_df)[1]*.1)
ps_prune <- prune_taxa(vec_taxa > cut_sample, ps)

######################################################################################
# NETWORK ANALYSIS ###################################################################
######################################################################################

########################
### FIT TWO NETWORKS ###
########################

# gut_split <- metagMisc::phyloseq_sep_variable(ps_prune, "W")

# net_W <- netConstruct(gut_split$`0`, gut_split$`1`, verbose = 2,
#                            # filtTax = "highestVar",
#                            # filtTaxPar = list(highestVar = 50),
#                            measure = "spieceasi",
#                            measurePar = list(method = "mb",
#                                              nlambda=20,
#                                              pulsar.params=list(rep.num=50)),
#                            normMethod = "none", zeroMethod = "none",
#                            sparsMethod = "none", seed = 123456)

# save(net_W, file = 'Pipeline_Microbiome/5_networks/net_W_output.RData')

load('data_pipleine_microbiome/net_W_output.RData')
rm(ps, sample_df, taxon_assign, matched_df, ASV_table)

#################################
### ANALYZE AND PLOT NETWORKS ###
#################################

props_W <- netAnalyze(net_W, clustMethod = "cluster_fast_greedy")
# go for the fast_greedy_algo

phyl_ps_prune <- as.factor(tax_table(ps_prune)[, 'Phylum'])
names(phyl_ps_prune) <- rownames(tax_table(ps_prune)[, 'Phylum'])
phylcol <- rainbow(length(unique(phyl_ps_prune)))

plot(props_W, sameLayout = TRUE, layoutGroup = 1, labels = FALSE, 
     featVecCol = phyl_ps_prune, 
     nodeColor = "feature", colorVec = phylcol,
     nodeTransp = 40,
     borderCol = "gray40", highlightHubs = FALSE,
     nodeSize = "normCounts", cexNodes = 1, 
     edgeTranspLow = 80, edgeTranspHigh = 50,
     groupNames = c("Clean air", "Polluted air"))
legend("bottomright",legend=levels(phyl_ps_prune),
       col=phylcol, bty="n", text.col = phylcol,
       cex = .5) 

## modularity graph      
nclust <- as.numeric(max(names(table(props_W$clustering$clust2))))
col <- topo.colors(nclust)

plot(props_W, sameLayout = TRUE, layoutGroup = 1, 
     labels = list(props_W$clustering$clust1, props_W$clustering$clust2), 
     labelFont = 1, cexLabels = 2,
     featVecCol = phyl_ps_prune, nodeColor = "cluster", 
     colorVec = col,
     nodeTransp = 40, 
     borderCol = "gray40", highlightHubs = FALSE,
     nodeSize = "normCounts", cexNodes = 1, 
     edgeTranspLow = 80, edgeTranspHigh = 50,
     groupNames = c("Clean air", "Polluted air"))

## modularity "split" graph for one cluster    
names_clust <- as.character(names(which(props_W$clustering$clust1 == 7)))
labels_phyl <- substr(as.character(phyl_ps_prune[names_clust]), 1, 3)

plot(props_W, sameLayout = TRUE, layoutGroup = 1, 
     labels = labels_phyl, 
     labelFont = 1, cexLabels = 1.5,
     # featVecCol = phyl_ps_prune, 
     nodeColor = "cluster", 
     colorVec = col,
     nodeTransp = 40, 
     borderCol = "gray40", highlightHubs = FALSE,
     # nodeSize = "normCounts", cexNodes = 1, 
     nodeFilter = "names", 
     nodeFilterPar = names_clust,
     edgeTranspLow = 80, edgeTranspHigh = 50,
     groupNames = c("Clean air", "Polluted air"))

## compare membership
dat_graph <- data.frame(id = names(props_W$clustering$clust1), 
                        com_c = props_W$clustering$clust1,
                        com_p = props_W$clustering$clust2)

dat_graph[,2] <- paste0(dat_graph[,2], "_C")
dat_graph[,3] <- paste0(dat_graph[,3], "_P")

g <- graph.data.frame(dat_graph[,c(2,3)], directed = F)
V(g)$type <- V(g)$name %in% dat_graph[,2] #the first column of edges is TRUE type

V(g)$color <- V(g)$type
V(g)$color=gsub("FALSE","orange",V(g)$color)
V(g)$color=gsub("TRUE","lightblue",V(g)$color)

layout_bi <- layout_with_sugiyama(g, hgap = 30)

plot(g, vertex.size = 15, layout=layout_bi$layout)

# every edge is a bunch of microbes
inc_com_sum <- ddply(data.frame(dat_graph), .(com_c, com_p), 
                     summarise, sum = length(id))

g_sum <- graph.data.frame(inc_com_sum[,c(1,2)], directed = F)
V(g_sum)$type <- V(g_sum)$name %in% inc_com_sum[,1] #the first column of edges is TRUE type

# node size represents the number of OTUs in each module
table(props_W$clustering$clust1)
table(props_W$clustering$clust2)
# add community name (without T or C)
V(g_sum)$com_nr <- as.numeric(unlist(strsplit(V(g_sum)$name, "_"))[c(TRUE,FALSE)])
V(g_sum)$com_size[V(g_sum)$type == TRUE] <- table(props_W$clustering$clust1)[V(g_sum)$com_nr[V(g_sum)$type == TRUE]] # treated sizes
V(g_sum)$com_size[V(g_sum)$type == FALSE] <- table(props_W$clustering$clust2)[V(g_sum)$com_nr[V(g_sum)$type == FALSE]] # control sizes

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

#########################
### COMPARE NETWORKS ###
########################

comp_W <- netCompare(props_W, permTest = FALSE, verbose = FALSE)

summary(comp_W, showCentr = c("degree", "eigen"), numbTaxa = 5)

library(phyloseq)
library(ggplot2)
library(ggtree)
library(RColorBrewer)

setwd('/Users/alicesommer/Desktop/DACOMP_cluster/')
rank_names <- c( "ASV" , "Species" , "Genus"  , "Family" , "Order"  , "Class"  , "Phylum")

###########
# SMOKING #
###########
load('dacomp_results/dacomp_results_PM.RData')
head(sort(dacomp_results$p_value_adj[dacomp_results$rank == "Genus"]))

load('dacomp_results/dacomp_results_smoke.RData')

for (d in 1:7){
  condition <- dacomp_results$p_value_adj[dacomp_results$rank == rank_names[d]] <= 0.4
  print(paste(d, ":", sum(condition, na.rm = TRUE)))
  print(dacomp_results$p_value_adj[dacomp_results$rank == rank_names[d]][which(condition)])
  print(dacomp_results$raw_test_stat[dacomp_results$rank == rank_names[d]][which(condition)])
  ps <- readRDS(paste0("ps_filt_", d, ".rds"))
  print(unname(tax_table(ps)[which(condition),]))
}


#### Plot tree ####
load('/Users/alicesommer/Desktop/Bureau/DOCTORATE/data_pipeline_microbiome/dada2output/phylotree2020.RData')
# load microbiome data
# ps_smoke <- readRDS(paste0("/Users/alicesommer/Desktop/DACOMP_cluster/ps_filt_3_PM.rds"))
ps_smoke <- readRDS(paste0("/Users/alicesommer/Desktop/DACOMP_cluster/ps_filt_3.rds"))

ps_smoke_tree <- merge_phyloseq(ps_smoke, phy_tree(tGTR$tree))

tip_labels <- unname(tax_table(ps_smoke)[,"Genus"])
tip_labels[is.na(tip_labels)] <- paste0("Order_", unname(tax_table(ps_smoke)[is.na(tip_labels),"Order"]))

# for effect size
tip_effect <- data.frame(label = rownames(tax_table(ps_smoke)), 
                         effect = abs(dacomp_results$raw_test_stat[dacomp_results$rank == "Genus"]),
                         p_value_adj = dacomp_results$p_value_adj[dacomp_results$rank == "Genus"],
                         tip_labels = tip_labels)

tree_smoke <- ggtree(ps_smoke_tree, color = 'grey', 
                     # ladderize = FALSE,
                     branch.length='none', 
                     layout = "circular") 

# for PM change scale for phylum
phylcol <- c(brewer.pal(9, "Set1")[1:7], brewer.pal(9, "Greens")[4], brewer.pal(9, "Set1")[8:9])

tip_effect$tip_size <- NA
tip_effect$tip_size[tip_effect$p_value_adj < .1] <- 6
tip_effect$tip_size[tip_effect$p_value_adj >= .1] <- 5
tip_effect$tip_size[tip_effect$p_value_adj >= .2] <- 4
tip_effect$tip_size[tip_effect$p_value_adj >= .3] <- 3
tip_effect$tip_size[tip_effect$p_value_adj >= .4] <- 2
tip_effect$tip_size[tip_effect$p_value_adj >= .5] <- 1

# change p-value adj for Marvinbryantia
# tip_effect$p_value_adj[tip_effect$tip_labels == "Marvinbryantia"] <- 0.01
tip_effect$tip_size2 <- NA
tip_effect$tip_size2[tip_effect$tip_labels == "Marvinbryantia"] <- 6

tree_smoke_effect <- tree_smoke %<+% tip_effect +
                    scale_color_brewer(palette = "Set1") +
                    # scale_color_manual(values = phylcol) +
                    geom_tippoint(aes(color=Phylum, size = tip_size), alpha=0.7) +
                    geom_tippoint(aes(color=Phylum, size = tip_size2), alpha=0.7, 
                                  pch = 21, stroke = 1,
                                  bg="#FF7F00", col='black',
                                  show.legend = FALSE) +
                    scale_size(name = "p-value adj.", 
                                labels = c('>= .5', '< .5', '< .4', 
                                           '< .3', '< .2', '< .1') ) +
                    geom_tiplab(aes(label = tip_labels, color = Phylum), size = 3, offset = 2) +
                    xlim(0, 40) 
                    # theme(legend.position="bottom", legend.box="vertical", legend.margin=margin())

# p <- gheatmap(tree_smoke_effect, tip_p_value, offset=.5, width=.03) 

ggsave(tree_smoke_effect, file = "/Users/alicesommer/Desktop/Bureau/DOCTORATE/plots_pipeline_microbiome/phylo_tree_DA_smoke.png",
              dpi=300,
              width = 230,
              height = 210,
              units = "mm")
 

#############
# POLLUTION #
#############

load('dacomp_results/dacomp_results_PM.RData')

for (d in 1:7){
  condition <- dacomp_results$p_value_adj[dacomp_results$rank == rank_names[d]] <= 0.2
  print(paste(d, ":", sum(condition, na.rm = TRUE)))
  print(dacomp_results$p_value_adj[dacomp_results$rank == rank_names[d]][which(condition)])
  print(dacomp_results$raw_test_stat[dacomp_results$rank == rank_names[d]][which(condition)])
  ps <- readRDS(paste0("ps_filt_", d, "_PM.rds"))
  if(sum(condition > 0, na.rm = TRUE)){print(unname(tax_table(ps)[which(condition),]))}
}

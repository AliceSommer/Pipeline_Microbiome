library(phyloseq)
library(ggplot2)
library(ggtree)
library(pals)
library(DECIPHER)
library(phangorn)

setwd('/Users/alicesommer/Desktop/DACOMP_cluster/')
rank_names <- c( "ASV" , "Species" , "Genus"  , "Family" , "Order"  , "Class"  , "Phylum")

###########
# SMOKING #
###########

load('dacomp_results/dacomp_results_smoke.RData')

for (d in 1:7){
  condition <- dacomp_results$p_value_adj[dacomp_results$rank == rank_names[d]] <= 0.2
  print(paste(d, ":", sum(condition, na.rm = TRUE)))
  print(dacomp_results$p_value_adj[dacomp_results$rank == rank_names[d]][which(condition)])
  print(dacomp_results$raw_test_stat[dacomp_results$rank == rank_names[d]][which(condition)])
  ps <- readRDS(paste0("ps_filt_", d, ".rds"))
  print(unname(tax_table(ps)[which(condition),]))
}


#### Plot tree ####
load('/Users/alicesommer/Desktop/Bureau/DOCTORATE/data_pipeline_microbiome/dada2output/phylotree2020.RData')
# load microbiome data
ps_smoke <- readRDS(paste0("/Users/alicesommer/Desktop/DACOMP_cluster/ps_filt_3.rds"))

ps_smoke_tree <- merge_phyloseq(ps_smoke, phy_tree(tGTR$tree))

tip_labels <- unname(tax_table(ps_smoke)[,"Genus"])
tip_labels[is.na(tip_labels)] <- paste0("Order_", unname(tax_table(ps_smoke)[is.na(tip_labels),"Order"]))

# for effect size
tip_effect <- data.frame(label = rownames(tax_table(ps_smoke)), 
                         effect = abs(dacomp_results$raw_test_stat[dacomp_results$rank == "Genus"]),
                         p_value_adj = dacomp_results$p_value_adj[dacomp_results$rank == "Genus"],
                         tip_labels = tip_labels)

# # for p-value
# tip_p_value_heat <- data.frame(p_value = dacomp_results$p_value_adj[dacomp_results$rank == "Genus"])
# rownames(tip_p_value) <- rownames(tax_table(ps_smoke))

tree_smoke <- ggtree(ps_smoke_tree, color = 'grey', 
                     # ladderize = FALSE,
                     branch.length='none', 
                     layout = "circular") 
                     
tree_smoke_effect <- tree_smoke %<+% tip_effect +
                    geom_tippoint(aes(color=Phylum, size = p_value_adj), alpha=0.7) +
                    scale_size(trans = 'reverse') +
                      # aes(color=Phylum)
                      # scale_color_viridis_d()
                    scale_color_brewer(palette = "Set1") +
                    # scale_color_manual(values=phylcol, labels=levels(tree_smoke$data$Phylum))
                    geom_tiplab(aes(label = tip_labels, color = Phylum), size = 3, offset = 2) +
                    xlim(0, 50) 
                    # theme(legend.position="bottom")

# p <- gheatmap(tree_smoke_effect, tip_p_value, offset=.5, width=.03) 

ggsave(tree_smoke_effect, file = "/Users/alicesommer/Desktop/phylo_tree_DA_smoke.png",
              dpi=300,
              width = 300,
              height = 350,
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

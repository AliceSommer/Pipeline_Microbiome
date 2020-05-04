library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(ggtree)

###############################################################################

# set working directory
setwd('/Users/alicesommer/Desktop/Bureau/DOCTORATE/data_pipeline_microbiome')

# load microbiome data
ASV_table <- readRDS('dada2output/seqtab2020.rds')
ASV_table_filtered <- readRDS('dada2output/seqtab2020_filtered.rds')
taxon_assign <- readRDS('dada2output/taxa2020.rds')
# load phylogenetic information
load("dada2output/phylotree2020.phy")

###############################################################################

ps_big <- phyloseq(otu_table(ASV_table, taxa_are_rows = FALSE), ## small OTU TABLE
                   tax_table(taxon_assign),
                   phy_tree(tGTR$tree))
ps_big

tree_big <- ggtree(ps_big, aes(color=Phylum), branch.length='none', layout = "circular") 
tree_big <- tree_big + theme(legend.position="bottom")

# save(tree_big, file='big_phylo_tree20.rds') 

# ggsave(file = 'big_phylo_plot20.jpeg',
#               tree_big,
#               dpi=300,
#               width = 300,
#               height = 200,
#               units = "mm")

ps_small <- phyloseq(ASV_table_filtered, 
                   tax_table(taxon_assign),
                   phy_tree(tGTR$tree))
ps_small

tree_small <- ggtree(ps_small, aes(color=Phylum), branch.length='none', layout = "circular") 
tree_small <- tree_small + theme(legend.position="bottom")

# save(tree_small, file='small_phylo_tree20.rds') 

# ggsave(file = 'figures/small_phylo_plot20.jpeg',
#               tree_small,
#               dpi=300,
#               width = 300,
#               height = 200,
#               units = "mm")

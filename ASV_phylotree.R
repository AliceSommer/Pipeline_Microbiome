library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(ggtree)

# set working directory
setwd('/Volumes/GoogleDrive/My\ Drive/DOCTORATE/Thesis')

# load microbiome data
ASV_table <- readRDS('KORA\ DATA/Microbiome_data/dada2output/dada2output2020/seqtab2020.rds')
ASV_table_filtered <- readRDS('KORA\ DATA/Microbiome_data/dada2output/dada2output2020/seqtab2020_filtered.rds')
taxon_assign <- readRDS('KORA\ DATA/Microbiome_data/dada2output/dada2output2020/taxa2020.rds')
# load phylogenetic information
load("KORA\ DATA/Microbiome_data/dada2output/dada2output2020/phylotree2020.phy")

ps_big <- phyloseq(otu_table(ASV_table, taxa_are_rows = FALSE), ## small OTU TABLE
                   tax_table(taxon_assign),
                   phy_tree(tGTR$tree))
ps_big

tree_big <- ggtree(ps_big, aes(color=Phylum), branch.length='none', layout = "circular") 
tree_big <- tree_big + theme(legend.position="bottom")

# save(tree_big, file='/Volumes/GoogleDrive/My\ Drive/DOCTORATE/Thesis/Microbiome\ 2020/2.\ February/big_phylo_tree20.rds') 

# ggsave(file = '/Volumes/GoogleDrive/My\ Drive/DOCTORATE/Thesis/Microbiome\ 2020/2.\ February/big_phylo_plot20.jpeg',
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

# save(tree_small, file='/Volumes/GoogleDrive/My\ Drive/DOCTORATE/Thesis/Microbiome\ 2020/2.\ February/small_phylo_tree20.rds') 

# ggsave(file = '/Volumes/GoogleDrive/My\ Drive/DOCTORATE/Thesis/Microbiome\ 2020/2.\ February/small_phylo_plot20.jpeg',
#               tree_small,
#               dpi=300,
#               width = 300,
#               height = 200,
#               units = "mm")

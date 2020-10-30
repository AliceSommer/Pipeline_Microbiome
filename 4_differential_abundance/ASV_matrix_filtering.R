library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(reshape2)

###############################################################################

# set working directory
setwd('/Users/alicesommer/Desktop/Bureau/DOCTORATE/data_pipeline_microbiome')

# load microbiome data
ASV_table <- readRDS('dada2output/seqtab2020.rds')
taxon_assign <- readRDS('dada2output/taxa2020.rds')
# load phylogenetic information
load("dada2output/phylotree2020.phy")

###############################################################################

# create a phyloseq object
ps <- phyloseq(otu_table(ASV_table, taxa_are_rows=FALSE),
               #sample_data(sample_df),
               tax_table(taxon_assign),
               phy_tree(tGTR$tree))
ps

par(mfrow=c(2,2))

# nr. of ASV present
locate_ASV_in_sample <- apply(otu_table(ps), 1, function(x) sum(x != 0))
hist(locate_ASV_in_sample, breaks = 30, main = "", xlab = "# ASV (sample)", xlim = c(0,400))

## sequencing depth (count statistics accross each n)
hist(sample_sums(ps), breaks = 50, main = "", xlab = "Sequencing Depth (sample)")
min(sample_sums(ps))

## count statistics accross each p
hist(taxa_sums(ps), breaks = 7000, main = "", xlab = "Total Counts (taxa)", xlim=c(1,7000))

## zero dist accross each p
zero_p <- apply(otu_table(ps), 2, function(x) sum(x == 0))
hist(zero_p, breaks = 100, main = "", xlab = "# Zeros (taxa)")

# apply function on every row of taxa table
names(taxon_assign[1,])[is.na(taxon_assign[1,])][1]

locate_NA_taxa <- apply(taxon_assign, 1, function(x) which(is.na(x))[1])
table(locate_NA_taxa)
sum(table(locate_NA_taxa))

locate_NA_taxa <- apply(taxon_assign, 1, function(x) sum(is.na(x)))
table(locate_NA_taxa)
sum(table(locate_NA_taxa))

########################
# PREVALENCE FILTERING #
########################

# load sample/matched_data
# load('dat_matched_PM25_bis.RData')
load('dat_matched_smoke_bis.RData')

sample_df <- matched_df[order(matched_df$ff4_prid),]
sample_df$W <- as.factor(sample_df$W)
samples.out <- as.character(sample_df$ff4_prid)
rownames(sample_df) <- samples.out

ps_PM25 <- phyloseq(otu_table(ASV_table, taxa_are_rows=FALSE),
                    sample_data(sample_df),
                    tax_table(taxon_assign))
ps_PM25

# Data
X = as(otu_table(ps_PM25), "matrix")

#### RESOLUTION ####

perc_seq <- seq(0,0.15,0.05)
ASV_cut <- matrix(NA, nrow = 1, ncol = length(perc_seq))
ASV_nr <- matrix(NA, nrow = 1, ncol = length(perc_seq))
count_per_sample <- matrix(NA, nrow = dim(X)[1], ncol = length(perc_seq))
ASV_per_sample <- matrix(NA, nrow = dim(X)[1], ncol = length(perc_seq))

vec_taxa <- apply(otu_table(ps_PM25), 2, function(x) sum(x > 0, na.rm = TRUE))

for (i in 1:length(perc_seq)){
  print(i)
  
  # how many tax are obs. in at least x% of samples
  perc <- perc_seq[i] # x%
  ASV_cut[i] <- trunc(dim(sample_df)[1]*perc)
  ASV_nr[i] <- length(which(vec_taxa > ASV_cut[i]))
  
  # add NA for all taxa that do not satisfy the filtering 
  ref <- which(vec_taxa <= ASV_cut[i])
  
  for (col in 1:length(ref)){
    c <- ref[col]
    X[,c] <- rep(NA, length(sample_df$W))
  }
  
  count_per_sample[,i] <- apply(X, 1, function(x) sum(x, na.rm = TRUE))
  ASV_per_sample[,i] <- apply(X, 1, function(x) sum(x > 0, na.rm = TRUE))
  
}

dim(count_per_sample)

dat_plot <- melt(count_per_sample)
colnames(dat_plot) <- c('sample', 'filter', 'values')
dat_plot$filter <- as.factor(dat_plot$filter)

dat_plot_2 <- melt(ASV_per_sample)
colnames(dat_plot_2) <- c('sample', 'filter', 'values')
dat_plot_2$filter <- as.factor(dat_plot_2$filter)

# prev_labels <- list(
#   '1'="0 % - p = 4,370",'2'="5 % - p = 515",
#   '3'="10 % - p = 311",'4'="15 % - p = 222"
# )

prev_labels <- list(
  '1'="0 % - p = 7,409",'2'="5 % - p = 483",
  '3'="10 % - p = 296",'4'="15 % - p = 209"
)

prev_labeller <- function(variable,value){
  return(prev_labels[value])
}

dat_text_lab <- data.frame(filter = as.factor(unique(dat_plot$filter)))
dat_text_lab$min_count <- NA
dat_text_lab$min_ASV <- NA

for (l in 1:4){
  print(l)
  print(summary(count_per_sample[,l]))
  print(summary(count_per_sample[,l])[1])
  
  dat_text_lab$min_count[l] <- summary(count_per_sample[,l])[1]
  dat_text_lab$min_ASV[l] <- summary(ASV_per_sample[,l])[1]
}

g_1 <- ggplot(dat_plot, aes(x = values)) +
  geom_histogram(fill="white",colour="black") +
  scale_x_continuous(limits = c(0, 60000), breaks = seq(0, 80000, 5000)) +
  facet_wrap(~filter, labeller=prev_labeller) +
  theme(axis.text.x = element_text(angle=90)) +
  geom_text(data = dat_text_lab, mapping = aes(x = 5000, y = 120, label = min_count), 
            colour = "red") + 
  geom_vline(data = dat_text_lab, mapping = aes(xintercept = min_count), 
             linetype = "dashed", colour = "red", size = .3) +
  xlab("ASV counts/sample")

g_2 <- ggplot(dat_plot_2, aes(x = values)) +
  geom_histogram(fill="white",colour="black") +
  scale_x_continuous(limits = c(0, 300), breaks = seq(0, 300, 20)) +
  facet_wrap(~filter, labeller=prev_labeller) +
  theme(axis.text.x = element_text(angle=90)) +
  geom_text(data = dat_text_lab, mapping = aes(x = 10, y = 90, label = min_ASV),
            colour = "red") +
  geom_vline(data = dat_text_lab, mapping = aes(xintercept = min_ASV),
             linetype = "dashed", colour = "red", size = .3) +
  xlab("observed ASVs/sample")

# ggsave(file = '/Users/alicesommer/Desktop/Bureau/DOCTORATE/plots_pipeline_microbiome/filter_count.jpeg', g_1,
#        dpi=300,
#        width = 170,
#        height = 170,
#        units = "mm")
# 
# ggsave(file = '/Users/alicesommer/Desktop/Bureau/DOCTORATE/plots_pipeline_microbiome/filter_ASV_nr.jpeg', g_2,
#        dpi=300,
#        width = 170,
#        height = 170,
#        units = "mm")

# ggsave(file = '/Users/alicesommer/Desktop/Bureau/DOCTORATE/plots_pipeline_microbiome/filter_count_smoke.jpeg', g_1,
#        dpi=300,
#        width = 170,
#        height = 170,
#        units = "mm")
# 
# ggsave(file = '/Users/alicesommer/Desktop/Bureau/DOCTORATE/plots_pipeline_microbiome/filter_ASV_nr_smoke.jpeg', g_2,
#        dpi=300,
#        width = 170,
#        height = 170,
#        units = "mm")




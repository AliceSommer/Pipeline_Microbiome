library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(DivNet)
library(doParallel)
library(doSNOW)
library(gridExtra)

###############################################################################

# set working directory
setwd('/Users/alicesommer/Desktop/Bureau/DOCTORATE/data_pipeline_microbiome')

# load microbiome data
ASV_table <- readRDS('dada2output/seqtab2020.rds')
taxon_assign <- readRDS('dada2output/taxa2020.rds')
# load phylogenetic information
load("dada2output/phylotree2020.RData")

# load sample/matched_data
# load('dat_matched_PM25_bis.RData')
load('dat_matched_smoke_bis.RData')

# load W matrix for randomization test
# load("W_paired_PM25.Rdata")
load("W_paired_smoke_bis.Rdata")

sample_df <- matched_df[order(matched_df$ff4_prid),]
sample_df$W <- as.factor(sample_df$W)
samples.out <- as.character(sample_df$ff4_prid)
rownames(sample_df) <- samples.out

################################################################################

# create a phyloseq object
ps <- phyloseq(otu_table(ASV_table, taxa_are_rows=FALSE),
               sample_data(sample_df),
               tax_table(taxon_assign),
               phy_tree(tGTR$tree))
ps

# locate the species that are totally absent in the matched data
empty_species <- colSums(otu_table(ps))
length(which(empty_species == 0))

ps_prune <- prune_taxa(empty_species != 0, ps)

# ps_prune <- subset_samples(ps_prune, pair_nb != 46)

### 1. ESTIMATE THE TOTAL ALPHA DIVERSITY FOR EACH SAMPLE ###

# agglomerate data to family level
ps_fam <- tax_glom(ps_prune, taxrank="Genus", NArm = FALSE)

shannon_plug <- sample_shannon(ps_fam)

base_abundant_taxa <- rownames(tax_table(ps_fam)[which(tax_table(ps_fam)[,"Genus"] == "Blautia")])

divnet_phylum <- divnet(ps_fam, base = base_abundant_taxa,
                         ncores = 4)
divnet_phylum

sample_data(ps_prune)[,"DivNet_W"] <- summary(divnet_phylum$shannon)$estimate
sample_data(ps_prune)[,"DN_error_W"] <- summary(divnet_phylum$shannon)$error
sample_data(ps_prune)[,"lower"] <- summary(divnet_phylum$shannon)$lower
sample_data(ps_prune)[,"upper"] <- summary(divnet_phylum$shannon)$upper
# add sequencing depth as variable
sample_data(ps_prune)[,"sample_counts"] <- sample_sums(ps_prune)
sample_data(ps_prune)[,"shannon_plug"] <- summary(shannon_plug)$estimate

sample_data(ps_prune)[,"dim"] <- as.factor(1:dim(sample_df)[1])

ggplot(sample_data(ps_prune), aes(x = dim)) + 
  geom_point(aes(y = DivNet_W), colour = "red") + 
  geom_errorbar(aes(ymin=lower, ymax=upper)) +
  geom_point(aes(y = shannon_plug), colour = "blue", alpha = .4) +
  xlab('Samples') + ylab("Shannon estimate")


#############
## SHANNON ##
#############

g_PM <- ggplot(sample_data(ps_prune), aes(color = factor(W), y = DivNet_W)) +
  geom_boxplot(alpha = .5) + ylab('DivNet shannon index') +
  scale_x_discrete(name = "") +
  # scale_colour_manual(values = c("darkgreen","blue4"), limits=c("1","0"), name ="Long-term PM2.5", labels = c("Low","High")) +
  scale_colour_manual(values = c("darkgreen","blue4"), limits=c("1","0"), name ="Smoking", labels = c("No","Yes")) +
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) +
  guides(color=guide_legend(nrow=2,byrow=TRUE))  + 
  annotate(geom="text",x=.7, y=1.7, label="stat. = 0.1037") +
  annotate(geom="text",x=.7, y=1.8, label="p-value = 0.0388") 

# PM: p-value = 0.0388; test-statistic = 0.1036517 

# ggsave(g_PM, file = "/Users/alicesommer/Desktop/Bureau/DOCTORATE/plots_pipeline_microbiome/box_shan_PM.png",
#        dpi=300,
#        width = 85,
#        height = 120,
#        units = "mm")

g_PM_sex <- ggplot(sample_data(ps_prune), aes(x = factor(u3csex), y = DivNet_W, color = factor(W))) +
  geom_boxplot(alpha = .5) + ylab('DivNet shannon index') +
  scale_x_discrete(name = "Sex", limits=c("0","1"), labels = c("Female", "Male")) +
  scale_colour_manual(values = c("gray","blue4"), limits=c("1","0"), name ="Long-term PM2.5", labels = c("Low","High")) +
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) +
  guides(color=guide_legend(nrow=2,byrow=TRUE))

g_arrange <- grid.arrange(g_PM,g_PM_sex, nrow = 1)

### 2. USE BETTA FUNCTION ###
x <- cbind(1, 
           # sample_data(ps_prune)$u3tcigsmk1, 
           # sample_data(ps_prune)$u3csex
           sample_data(ps_prune)$W )

reg <- betta(summary(divnet_phylum$shannon)$estimate,
             summary(divnet_phylum$shannon)$error, X = x)
reg$table
estim_obs <- reg$table[2,1]

### 3. PERFORM A RANDOMIZATION TEST ###
W_paired <- W_paired_smoke
dim(W_paired)

# set the number of randomizations
nrep <- ncol(W_paired)/100

# create a matrix where the t_rand will be saved
t_array <- NULL

for(i in 1:nrep){
  print(i)
  x = cbind(1, 
            # sample_data(ps_prune)$u3tcigsmk1,
            # sample_data(ps_prune)$u3csex
            W_paired[,i] )
  
  reg = betta(summary(divnet_phylum$shannon)$estimate,
              summary(divnet_phylum$shannon)$error, X = x)
  
  # fill t_array
  t_array[i] = reg$table[2,1] 
}

## calculate p_value
p_value <- mean(t_array >= estim_obs)
p_value

# 1. Open png file
png("/Users/alicesommer/Desktop/Bureau/DOCTORATE/plots_pipeline_microbiome/beta_shan_PM.png", width = 350, height = 350)
# 2. Create the plot
hist(t_array, breaks = 30, main = "", xlab = "shannon beta (PM2.5)")
abline(v = estim_obs, col = 'red', lwd = 2, lty = 2)
# 3. Close the file
dev.off()




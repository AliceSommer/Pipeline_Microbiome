library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(breakaway)
library(dplyr)
library(gridExtra)

#############
# load data #
#############

## set working directory
setwd('/Users/alicesommer/Desktop/Bureau/DOCTORATE/data_pipeline_microbiome')

## load microbiome data
ASV_table <- readRDS('dada2output/seqtab2020.rds')
taxon_assign <- readRDS('dada2output/taxa2020.rds')
## load phylogenetic information
load("dada2output/phylotree2020.RData")

## load sample/matched_data
load('dat_matched_PM25_bis.RData')
# load('dat_matched_smoke_bis.RData')

## load W matrix for randomization test
load("W_paired_PM25.Rdata")
# load("W_paired_smoke_bis.Rdata")

############################
# create a phyloseq object # 
############################

## order covariates data frame
sample_df <- matched_df[order(matched_df$ff4_prid),]
sample_df$W <- as.factor(sample_df$W)
samples.out <- as.character(sample_df$ff4_prid)
rownames(sample_df) <- samples.out

# create a phyloseq object
ps <- phyloseq(otu_table(ASV_table, taxa_are_rows=FALSE),
               sample_data(sample_df),
               tax_table(taxon_assign))
ps

## locate the species that are totally absent in the matched data
empty_species <- colSums(otu_table(ps))
length(which(empty_species == 0))
# remove them
ps_prune <- prune_taxa(empty_species != 0, ps)

################################################
### 1. ESTIMATE THE RICHNESS FOR EACH SAMPLE ###
################################################

# calculate plug-in richness
rich <- sample_richness(ps_prune)

## estimate richness with breakaway function  
ba <- breakaway(ps_prune)
ba[[1]]
# plot(ba, ps_prune, color = "W")

## retrieve statistics for inference
sample_data(ps_prune)[,"breakaway_W"] <- summary(ba)$estimate
sample_data(ps_prune)[,"ba_error_W"] <- summary(ba)$error
sample_data(ps_prune)[,"lower"] <- summary(ba)$lower
sample_data(ps_prune)[,"upper"] <- summary(ba)$upper
sample_data(ps_prune)[,"sample_counts"] <- sample_sums(ps_prune)
sample_data(ps_prune)[,"richness"] <- summary(rich)$estimate

ggplot(sample_data(ps_prune), aes(x = dim)) + 
  geom_point(aes(y = breakaway_W), colour = "red") + 
  # geom_errorbar(aes(ymin=lower, ymax=upper)) +
  geom_point(aes(y = richness), colour = "blue", alpha = .5) +
  xlab('Samples') + ylab("BA estimate")

# hist(sample_data(ps_prune)$ba_error_W, breaks = 50)
# head(sample_data(ps_prune)[sample_data(ps_prune)$breakaway_W > 400, c("ba_error_W","sample_counts","breakaway_W","pair_nb")])
# head(sample_data(ps_prune)[sample_data(ps_prune)$pair_nb == 62, c("ba_error_W","sample_counts","breakaway_W","pair_nb")])
# 
# summary(sample_data(ps_prune)$sample_counts)
# 
# # remove the outlier for plot ?
# ps_prune_out <- subset_samples(ps_prune, pair_nb != 62)

g_PM <- ggplot(sample_data(ps_prune), aes(color = factor(W), y = breakaway_W)) +
  geom_boxplot(alpha = .5) + ylab('breakaway richness measure') +
  scale_x_discrete(name = "") +
  scale_colour_manual(values = c("darkgreen","blue4"), limits=c("1","0"), name ="Long-term PM2.5", labels = c("Low","High")) +
  # scale_colour_manual(values = c("darkgreen","blue4"), limits=c("1","0"), name ="Smoking", labels = c("No","Yes")) +
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) +
  guides(color=guide_legend(nrow=2,byrow=TRUE)) + 
  annotate(geom="text",x=.7, y=45, label="stat. = 19.1") +
  annotate(geom="text",x=.7, y=70, label=expression(p-value %~~% "0.0008")) 

# PM: p-value = 0.0008; test-statistic = 19.0775
# Smoke: p-value = 0.1518; test-statistic = 3.9576

# ggsave(g_PM, file = "/Users/alicesommer/Desktop/Bureau/DOCTORATE/plots_pipeline_microbiome/box_break_PM.png",
#        dpi=300,
#        width = 85,
#        height = 120,
#        units = "mm")

# plot data stratified by sex
g_PM_sex <- ggplot(sample_data(ps_prune), aes(x = factor(u3csex), y = breakaway_W, color = factor(W))) +
  geom_boxplot(alpha = .5) + ylab('breakaway richness measure') +
  scale_x_discrete(name = "Sex", limits=c("0","1"), labels = c("Female", "Male")) +
  scale_colour_manual(values = c("gray","blue4"), limits=c("1","0"), name ="Long-term PM2.5", labels = c("Low","High")) +
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) +
  guides(color=guide_legend(nrow=2,byrow=TRUE))

g_arrange <- grid.arrange(g_PM,g_PM_sex, nrow = 1)

#######################################################
### 2. RANDOMIZATION-TEST WITH BETTA TEST STATISTIC ###
#######################################################
x <- cbind(1, 
           # sample_data(ps_prune)$u3tcigsmk,
           # sample_data(ps_prune)$u3csex
          sample_data(ps_prune)$W)
head(sample_data(ps_prune)$breakaway_W)
head(sample_data(ps_prune)$ba_error_W)

reg <- betta(sample_data(ps_prune)$breakaway_W,
             sample_data(ps_prune)$ba_error_W, X = x)
reg$table
estim_obs <- reg$table[2,1]

## W matrix 
dim(W_paired_smoke)
# W_paired_smoke <- W_paired_smoke[-which(sample_data(ps_prune)$pair_nb == 62), ]
# dim(W_paired)

# set the number of randomizations
nrep <- ncol(W_paired_smoke)/100
# nrep <- ncol(W_paired)/100

# create a matrix where the t_rand will be saved
t_array <- NULL

for(i in 1:nrep){
  print(i)
  x = cbind(1, 
            # sample_data(ps_prune)$u3tcigsmk, 
            # sample_data(ps_prune)$u3csex
            W_paired_smoke[,i])
            # W_paired[,i]) 
            
  
  reg = betta(sample_data(ps_prune)$breakaway_W,
              sample_data(ps_prune)$ba_error_W, X = x)
  
  # fill t_array
  t_array[i] = reg$table[2,1] 
}

## calculate p_value
p_value <- mean(t_array >= estim_obs, na.rm = TRUE)
p_value

## plot distribution of test-statistic
# 1. Open png file
png("/Users/alicesommer/Desktop/Bureau/DOCTORATE/plots_pipeline_microbiome/beta_break_smoke.png", width = 350, height = 350)
# 2. Create the plot
hist(t_array, breaks = 30, main = "", xlab = "breakaway beta (Smoking)")
abline(v = estim_obs, col = 'red', lwd = 2, lty = 2)
# 3. Close the file
dev.off()

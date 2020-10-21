library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

###############################################################################

# set working directory
setwd('/Users/alicesommer/Desktop/Bureau/DOCTORATE/data_pipeline_microbiome')

# load microbiome data
ASV_table <- readRDS('dada2output/seqtab2020.rds')
taxon_assign <- readRDS('dada2output/taxa2020.rds')

# load sample/matched_data
# load('dat_matched_PM25_bis.RData')
load('dat_matched_smoke_bis.RData')

# # load original dataset (try other Ws)
# df <- read.csv('/Volumes/GoogleDrive/My\ Drive/DOCTORATE/Thesis/KORA\ DATA/Microbiome_data/KORA_microbiome_variables.csv')

# load W matrix for randomization test
# load("W_paired_PM25.Rdata")
load("W_paired_smoke_bis.Rdata")

# matched AP data
sample_df <- matched_df[order(matched_df$ff4_prid),]
sample_df$W <- as.factor(sample_df$W)
# # "other" data
# sample_df <- df[order(df$ff4_prid),]
# sample_df$W <- as.factor(as.numeric(sample_df$u3tcigsmk == 1))

samples.out <- as.character(sample_df$ff4_prid)
rownames(sample_df) <- samples.out
################################################################################

# create a phyloseq object
ps <- phyloseq(otu_table(ASV_table, taxa_are_rows=FALSE),
               sample_data(sample_df),
               tax_table(taxon_assign))
ps

# locate the species that are totally absent in the matched data 
# 1. have a vector of nr. of observed samples per taxa
vec_taxa <- apply(otu_table(ps), 2, function(x) sum(x > 0, na.rm = TRUE))
length(which(vec_taxa == 0))
# 2. how many tax are obs. in at least x% of samples
perc <- 0.10 # x%
# samp_perc <- trunc(dim(sample_df)[1]*perc)
samp_perc <- 0
length(which(vec_taxa > samp_perc))

ps_prune <- prune_taxa(vec_taxa > samp_perc, ps)

# ## agglomerate to Species ##
# ps_Species <- tax_glom(ps, taxrank = "Species", NArm = FALSE)
# vec_taxa_Sp <- apply(otu_table(ps_Species), 2, function(x) sum(x > 0, na.rm = TRUE))
# length(which(vec_taxa_Sp > samp_perc))
# ps_Species <- prune_taxa(vec_taxa_Sp >= samp_perc, ps_Species)
# ps_Species

## agglomerate to Genus ##
ps_Genus <- tax_glom(ps_prune, taxrank = "Genus", NArm = FALSE)
vec_taxa_Gen <- apply(otu_table(ps_Genus), 2, function(x) sum(x > 0, na.rm = TRUE))
table(vec_taxa_Gen)
ps_Genus_prune <- prune_taxa(vec_taxa_Gen >= samp_perc, ps_Genus)
ps_Genus_prune

## agglomerate to Order ##
ps_Order <- tax_glom(ps_prune, taxrank = "Order", NArm = FALSE)

ps <- ps_Genus_prune

# Filter the data
n<-dim(otu_table(ps))[1]
p<-dim(otu_table(ps))[2]

# Add 0.5 on zeros
otu_table(ps)[otu_table(ps)==0] <- 0.5

# Separate the sample into two groups
ps_control = subset_samples(ps, sample_data(ps)$W == 0)
ps_treated = subset_samples(ps, sample_data(ps)$W == 1)

X_c = as(otu_table(ps_control), "matrix"); dim(X_c)
X_t = as(otu_table(ps_treated), "matrix"); dim(X_t)

nx <- dim(X_c)[1]
ny <- dim(X_t)[1]

x <- X_c/(rowSums(X_c)%*%matrix(1,1,p))
y <- X_t/(rowSums(X_t)%*%matrix(1,1,p))

log_x<-log(x)
log_y<-log(y)

clog_TX <- log_x-1/p*rowSums(log_x)%*%matrix(1,1,p)
clog_TY <- log_y-1/p*rowSums(log_y)%*%matrix(1,1,p)

# Mean
# x_mean <- colSums(x)/nx
# y_mean <- colSums(y)/ny
# 
# log_x_mean <- colSums(log_x)/nx
# log_y_mean <- colSums(log_y)/ny

clr_x_mean <- colSums(clog_TX)/nx
clr_y_mean <- colSums(clog_TY)/ny

# Variance
# x_var <- diag(var(x))*(nx-1)/nx
# y_var <- diag(var(y))*(ny-1)/ny
# x_stat_var <- (x_var*nx + y_var*ny)/(nx*ny)
# 
# log_x_var <- diag(var(log_x))*(nx-1)/nx
# log_y_var <- diag(var(log_y))*(ny-1)/ny
# log_x_stat_var <- (log_x_var*nx + log_y_var*ny)/(nx*ny)

clr_x_var <- diag(var(clog_TX))*(nx-1)/nx
clr_y_var <- diag(var(clog_TY))*(ny-1)/ny
clr_x_stat_var <- (clr_x_var*nx + clr_y_var*ny)/(nx*ny)

# p-values
# x_stat <- max((x_mean/sqrt(x_stat_var) - y_mean/sqrt(x_stat_var))^2)
# p_x <- 1-exp(-1/sqrt(pi)*exp(-(x_stat-(2*log(p)-log(log(p))))/2))
# log_x_stat <- max((log_x_mean/sqrt(log_x_stat_var) - log_y_mean/sqrt(log_x_stat_var))^2)
# p_logx <- 1-exp(-1/sqrt(pi)*exp(-(log_x_stat-(2*log(p)-log(log(p))))/2))
clr_x_stat <- max((clr_x_mean/sqrt(clr_x_stat_var) - clr_y_mean/sqrt(clr_x_stat_var))^2)
p_clrx <- 1-exp(-1/sqrt(pi)*exp(-(clr_x_stat-(2*log(p)-log(log(p))))/2))
clr_x_stat; p_clrx

####
#### FISHER
####

matrix_otu_bind <- rbind(clog_TX, clog_TY)
dim(matrix_otu_bind)

# set the number of randomizations
W_paired <- W_paired_smoke
nrep <- ncol(W_paired)/100

Tarray = NULL

for (i in 1:nrep){
  W_suffle <- W_paired[,i]
  
  clog_TX_rand <- matrix_otu_bind[W_suffle == 0,]
  clog_TY_rand <- matrix_otu_bind[W_suffle == 1,]
  
  clr_x_mean_rand <- colSums(clog_TX_rand)/nx
  clr_y_mean_rand <- colSums(clog_TY_rand)/ny
  
  clr_x_var_rand <- diag(var(clog_TX_rand))*(nx-1)/nx
  clr_y_var_rand <- diag(var(clog_TY_rand))*(ny-1)/ny
  clr_x_stat_var_rand <- (clr_x_var_rand*nx + clr_y_var_rand*ny)/(nx*ny)
  
  clr_x_stat_rand <- max((clr_x_mean_rand/sqrt(clr_x_stat_var_rand) - clr_y_mean_rand/sqrt(clr_x_stat_var_rand))^2)
  
  Tarray = c(Tarray, clr_x_stat_rand)
  
  print(i)
}

hist(Tarray, main = NULL, xlab = NULL, breaks = 100)
abline(v=clr_x_stat, lty = 2, lwd = 2)
pval = mean(Tarray >= clr_x_stat); pval


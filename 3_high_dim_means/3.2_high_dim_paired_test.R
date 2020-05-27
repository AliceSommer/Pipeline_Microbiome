library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

###############################################################################

# set working directory
setwd('/Users/alicesommer/Desktop/Bureau/DOCTORATE/data_pipeline_microbiome')

# load microbiome data
ASV_table <- readRDS('dada2output/seqtab2020.rds')
taxon_assign <- readRDS('dada2output/taxa2020.rds')

# load sample/matched_data
load('dat_matched_PM25_bis.RData')

# # load original dataset (try other Ws)
# df <- read.csv('/Volumes/GoogleDrive/My\ Drive/DOCTORATE/Thesis/KORA\ DATA/Microbiome_data/KORA_microbiome_variables.csv')

# load W matrix for randomization test
load("W_paired_PM25.Rdata")

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
empty_species <- colSums(otu_table(ps))
length(which(empty_species == 0))

ps_prune <- prune_taxa(empty_species != 0, ps)

# ps = ps_prune
ps = tax_glom(ps_prune, "Genus", NArm = FALSE)

# Add 0.5 on zeros
otu_table(ps)[otu_table(ps)==0] <- 0.5

# Separate the sample into two groups
ps_control = subset_samples(ps, sample_data(ps)$W == 0)
ps_treated = subset_samples(ps, sample_data(ps)$W == 1)

X_c = as(otu_table(ps_control), "matrix"); dim(X_c)
X_c = X_c[order(sample_data(ps_control)$pair_nb),] # reorder to have pairs matched !!!!!
X_t = as(otu_table(ps_treated), "matrix"); dim(X_t)

p = dim(X_c)[2]
n = dim(X_c)[1]

# Covert to the composition
X_c <- X_c/(rowSums(X_c)%*%matrix(1,1,p))
X_t <- X_t/(rowSums(X_t)%*%matrix(1,1,p))

# Logratio and clr
logX_c <- log(X_c)
logX_t <- log(X_t)

clrX_c <- logX_c-1/p*rowSums(logX_c)%*%matrix(1,1,p)
clrX_t <- logX_t-1/p*rowSums(logX_t)%*%matrix(1,1,p)

# Squared standardized difference
SqStandizeDiffX_t_c <- ((colSums(X_c-X_t))/n)^2/((diag(var(X_c-X_t)))*(n-1)/n^2)

SqStandizeDiffLogX_t_c <- ((colSums(logX_c-logX_t))/n)^2/((diag(var(logX_c-logX_t)))*(n-1)/n^2)

SqStandizeDiffclrX_t_c <- ((colSums(clrX_c-clrX_t))/n)^2/((diag(var(clrX_c-clrX_t)))*(n-1)/n^2)

# p-values
StatX_t_c <- max(SqStandizeDiffX_t_c)

StatLogX_t_c <- max(SqStandizeDiffLogX_t_c)

StatclrX_t_c <- max(SqStandizeDiffclrX_t_c)

pX_t_c <- 1-exp(-1/sqrt(pi)*exp(-(StatX_t_c-(2*log(p)-log(log(p))))/2))

pLogX_t_c <- 1-exp(-1/sqrt(pi)*exp(-(StatLogX_t_c-(2*log(p)-log(log(p))))/2))

pclrX_t_c <- 1-exp(-1/sqrt(pi)*exp(-(StatclrX_t_c-(2*log(p)-log(log(p))))/2))

StatclrX_t_c; pclrX_t_c

####
#### FISHER
####

# set the number of randomizations
nrep <- ncol(W_paired)/100

Tarray = NULL

big_X = as(otu_table(ps), "matrix")
big_X = big_X[order(sample_df$pair_nb),]

W_paired_order = W_paired[order(sample_df$pair_nb),] # reorder to have pairs matched !!!!!
head(W_paired_order[,1:4])

for (i in 1:nrep){
  
  W_shuffle <- W_paired[,i]
  
  X_c = big_X[W_shuffle == 0,]
  X_t = big_X[W_shuffle == 1,]
  
  # Logratio and clr
  logX_c <- log(X_c)
  logX_t <- log(X_t)
  
  clrX_c <- logX_c-1/p*rowSums(logX_c)%*%matrix(1,1,p)
  clrX_t <- logX_t-1/p*rowSums(logX_t)%*%matrix(1,1,p)
  
  SqStandizeDiffclrX_t_c <- ((colSums(clrX_c-clrX_t))/n)^2/((diag(var(clrX_c-clrX_t)))*(n-1)/n^2)
  
  Tarray[i] <- max(SqStandizeDiffclrX_t_c)
  
  print(i)
}

hist(Tarray, main = NULL, xlab = NULL, breaks = 100)
abline(v=StatclrX_t_c, lty = 2, lwd = 2)
pval = mean(Tarray >= StatclrX_t_c); pval

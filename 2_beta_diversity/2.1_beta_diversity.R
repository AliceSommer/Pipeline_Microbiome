
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(plyr); packageVersion("plyr")
library(MiRKAT)
library(compositions)
library(reshape2)
library(RColorBrewer)

#############
# load data #
#############
wd <- '/Users/alicesommer/Desktop/Bureau/DOCTORATE/data_pipeline_microbiome'

## set working directory
setwd(wd)

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
               tax_table(taxon_assign),
               phy_tree(tGTR$tree))
ps

## locate the species that are totally absent in the matched data
empty_species <- colSums(otu_table(ps))
length(which(empty_species == 0))
# remove them
ps_prune <- prune_taxa(empty_species != 0, ps)

##################
### 1. MiRKAT ###
#################

## MiRKAT
u_unifrac <- UniFrac(ps_prune, weighted=FALSE, parallel=FALSE, fast=TRUE)
unifrac_mat <- as.matrix(u_unifrac)
# transform distance matrix to kernel
K.unweighted_uni <- D2K(unifrac_mat)

head(sample_data(ps_prune)$W)
outcome <- as.numeric(sample_data(ps_prune)$W == 1)
head(outcome)

# ## testing using a single Kernel
# MiRKAT(y = outcome, X = NULL, Ks = K.unweighted_uni, out_type = "D", 
#                      method = "davies", returnKRV = TRUE, returnR2 = TRUE)

## omnibus test if multiple distance matrices ("Optimal MiRKAT")

## tranform the data
ps_clr <- transform_sample_counts(ps_prune, function(x){x <- x + 1; clr(x)})
ps_comp <- transform_sample_counts(ps_prune, function(x){x <- x + 1; x/sum(x)})

# the available distance methods coded in distance
dist_methods <- unlist(distanceMethodList)
print(dist_methods)

# Euclidean
euclidean_dist <- distance(ps_clr, method="euclidean")
K.euclidean_dist <- D2K(as.matrix(euclidean_dist))

# # Bray
# bray_dist <- distance(ps_prune, method="bray") 
# K.bray_dist <- D2K(as.matrix(bray_dist))

# Jaccard 
jaccard_dist <- distance(ps_prune, method="jaccard", binary = TRUE) 
K.jaccard_dist <- D2K(as.matrix(jaccard_dist))

# Gower
gower_dist <- distance(ps_clr, method="gower")
K.gower_dist <- D2K(as.matrix(gower_dist))

## testing using a several Kernels
Ks = list(u_unifrac = K.unweighted_uni, euclidean_dist = K.euclidean_dist, jaccard_dist = K.jaccard_dist, gower_dist = K.gower_dist)
MiRKAT_obs <- MiRKAT(y = outcome, Ks = Ks, X = NULL, out_type = "D", 
       method = "permutation", returnKRV = TRUE, returnR2 = TRUE)

#######################################
### 2. PERFORM A RANDOMIZATION TEST ###
######################################

#### ADAPT MiRKAT to get Q-stat ####
y = outcome; Ks = Ks; X = NULL; family = "binomial"

n <- length(y)
if (is.null(X)) {
  X1 <-  matrix(rep(1, length(y)), ncol=1)
} else {
  X1 <- model.matrix(~. , as.data.frame(X))
}

qX1 <- qr(X1)
## Take care of aliased variables and pivoting in rhs
X1 <- X1[, qX1$pivot, drop=FALSE]
X1 <- X1[, 1:qX1$rank, drop=FALSE]
options(warn=2)  # make sure this model is correct
mod <- glm(y ~ X1-1, family = family)
options(warn=1)

px  = NCOL(X1)
mu  = mod$fitted.values
res = y - mu  

# Continuous or binary outcome Q statistic 
getQ = function(K, res, s2){    
  Q = c(1 / s2 * res %*% K %*% res)
}

Qs_obs = lapply(Ks, getQ, res, s2 = 1)

## W matrix
W_paired <- W_paired_smoke
dim(W_paired)

# set the number of randomizations
nrep <- ncol(W_paired)/100

# create a matrix where the t_rand will be saved
t_arrays <- matrix(NA, ncol=length(Ks), nrow=nrep)

for(j in 1:nrep){
  print(j)
  
  y_new = W_paired[,j]
  
  mod <- glm(y_new ~ X1-1, family = family)
  
  mu  = mod$fitted.values
  res = y_new - mu  
  
  Qs_new = lapply(Ks, getQ, res, s2 = 1)
  
  # fill t_arrays 
  t_arrays[j,] = unlist(Qs_new) ## not sure this is the statistic we want (KRV?)
}

## add the observed test stats
t_arr <- rbind(unlist(Qs_obs), t_arrays)
head(t_arr)

## calculate p_values
p_values <- NULL
for (p in 1:length(Ks)){
  p_values[p] <- mean(t_arr[,p] >= Qs_obs[p])
}
p_values

MiRKAT_obs$p_values

## plot randomization distributions
t_arrays_data_frame <- data.frame(t_arrays)
colnames(t_arrays_data_frame) <- c("Unifrac", "Aitchison", "Jaccard", "Gower")

t_array_melt <- melt(t_arrays_data_frame)
# t_array_melt_plot <- t_array_melt[t_array_melt$variable != "Bray",]

dat_text_lab <- data.frame(variable = colnames(t_arrays_data_frame))
dat_text_lab$obs_stat <- as.numeric(Qs_obs)
# dat_text_lab <- dat_text_lab[dat_text_lab$variable != "Bray",]

g_rand <- ggplot(t_array_melt,aes(x = value)) +
  facet_wrap(~variable, scales = "free") +
  geom_histogram(fill="white",colour="black") + 
  geom_vline(data = dat_text_lab, mapping = aes(xintercept = obs_stat), 
             linetype = "dashed", colour = "red", size = .3) +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        panel.background = element_blank())

# ggsave(file = 'plots_pipeline_microbiome/null_dist_beta_smoke.png',
#        g_rand,
#        dpi=300,
#        width = 180,
#        height = 150,
#        units = "mm")

##########################################
### 3. MULTIPLE COMPARISON ADJUSTMENT ####
##########################################
verbose = T
### Lee, Forastiere, Miratix, Pillai (2017) method  for multiple comparison adjustment ### 

# STEP 1 to 3: recorded in "stats_matrix"
dim(t_arrays)
# the first row is the observed

# AT STEP 4: consider at each rep. T_rep as T_obs to calculate nrep p-values 
# for the hypothetical test statistics

hyp_matrix <- t_arr 
hyp_p_value <- matrix(NA, ncol = dim(t_arrays)[2], nrow = nrep)

# based on value (hyp_obs) of each row
for (r in 1:nrep){
  if(verbose)
    if(r%% ceiling(nrep/100) == 1)
      cat(paste0('Testing rep : ',r,'/',nrep,' \n\r'))
  # calc. hypothetical p_value on each column of the matrix 
  hyp_p_value[r,] <- apply(t_arr, 2, function(x) mean(x >= x[r]))
}

# for each rep. take the min. p_value
min_p_nrep <- apply(hyp_p_value, 1, function(x) min(x, na.rm = TRUE))
head(min_p_nrep)

# calculate the proportion of min_p_nrep that is sm/eq. p_value (for obs.)
p_value_adj <- sapply(p_values, function(x) mean(min_p_nrep <= x))
head(p_value_adj)

col_vline <- brewer.pal(4, "Dark2")

png("/Users/alicesommer/Desktop/Bureau/DOCTORATE/plots_pipeline_microbiome/padj_beta_smoke.png", width = 450, height = 350)
hist(min_p_nrep, main = "", xlab = "min. p-value for 10,000 rep.")
abline(v = p_values, col = col_vline, lty = 4, lwd = 1.5)
legend("topright", colnames(t_arrays_data_frame), text.col = col_vline)
dev.off()

###########################################
##### Low-dimensional representation ######
###########################################
# dist_methods = c("unifrac","euclidean","bray","jaccard",gower")
dist_methods = c("unifrac","euclidean","jaccard","gower") # without bray

plist <- vector("list", length(dist_methods))
names(plist) = dist_methods
for(i in dist_methods){
  print(i)
  # Calculate distance matrix
  if(i %in% c("unifrac","bray")){
  iDist <- distance(ps_prune, method=i)
  iMDS  <- ordinate(ps_prune, "MDS", distance=iDist)
  }else if(i %in% c("jaccard")){
  iDist <- distance(ps_prune, method=i, binary = TRUE)
  iMDS  <- ordinate(ps_prune, "MDS", distance=iDist)
  }else{
  iDist <- distance(ps_clr, method=i)
  iMDS  <- ordinate(ps_clr, "MDS", distance=iDist)
  }
  
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(ps_prune, iMDS, color="W")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}

#### attention error because of the bray distance !!!!

df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"

df$distance[df$distance == "euclidean"] <- "aitchison"

p = ggplot(df, aes(Axis.1, Axis.2, color=W))
p = p + geom_point(size=1, alpha=0.5)
p = p + facet_wrap(~distance, scales="free") + 
  scale_colour_manual(values = c("darkgreen","blue4"), 
                      # limits=c("1","0"), name ="Smoking", labels = c("No","Yes")) 
                      limits=c("1","0"), name ="Long-term PM2.5", labels = c("Low","High")) 
p = p + ggtitle("MDS on various distance metrics") + 
  theme(legend.position = "bottom", legend.key.size =  unit(0.1, "in"))
p

# ggsave(file = 'plots_pipeline_microbiome/iMDS_PM.png',
#        p,
#        dpi=300,
#        width = 150,
#        height = 150,
#        units = "mm")



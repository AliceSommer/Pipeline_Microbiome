library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(breakaway)
library(DivNet)


###############################################################################

# set working directory
setwd('/Users/alicesommer/Desktop/Bureau/DOCTORATE/data_pipeline_microbiome')

# load data formated in NOV18 
load('dat_transformed_NOV18.RData')

# load microbiome data
ASV_table <- readRDS('dada2output/seqtab2020.rds')
taxon_assign <- readRDS('dada2output/taxa2020.rds')

sample_df <- dat_transformed[order(dat_transformed$ff4_prid),]
samples.out <- as.character(sample_df$ff4_prid)
rownames(sample_df) <- samples.out

################################################################################

# create a phyloseq object
ps <- phyloseq(otu_table(ASV_table, taxa_are_rows=FALSE),
               sample_data(sample_df),
               tax_table(taxon_assign))
ps

## plug-in richness estimate
rich <- sample_richness(ps)

## agglomerate data to family level
ps_spe <- tax_glom(ps, taxrank="Genus")

###############
## BREAKAWAY ##
###############

ba <- breakaway(ps)

## add the estimate and error in the dataset
sample_df[,"breakaway"] <- summary(ba)$estimate
sample_df[,"ba_error"] <- summary(ba)$error

summary(ba)[summary(ba)$estimate > 300,]
##>>>>>>>>> 3 outliers to report 

############
## DivNet ##
############

# function to find the most abundant taxa
# Goes through a phyloseq object, picks out the most abundant taxa and gives the abundance for each
# and identifies which taxa is most abundant for which sample

# find.top.taxa <- function(x,taxa){
#   require(phyloseq)
#   top.taxa <- tax_glom(x, taxa)
#   otu <- otu_table(top.taxa) # remove the transformation if using a merge_sample object
#   tax <- tax_table(top.taxa)
#   j<-apply(otu,1,which.max)
#   k <- j[!duplicated(j)]
#   l <- data.frame(tax@.Data[k,]) # note the modification here and the line below
#   m <- data.frame(otu@.Data[,k])
#   s <- as.name(taxa)
#   colnames(m) = l[,taxa]
#   n <- colnames(m)[apply(m,1,which.max)]
#   m[,taxa] <- n
#   return(m)
# }
# top_taxa <- find.top.taxa(ps,"Phylum")
# table(top_taxa$Phylum)

ps_phylum <- tax_glom(ps, taxrank="Phylum")
ps_phylum

set.seed(16)
DN <- divnet(ps_phylum, ncores = 4, variance = "parametric")

## add the estimate and error in the dataset
sample_df[,"div_shannon"] <- summary(DN$`shannon`)$estimate
sample_df[,"div_shannon_error"] <- summary(DN$`shannon`)$error

head(sample_df[,"div_shannon"])
head(sample_df[,"Shannon.effective"])

# load data formated in MAY20
save(sample_df, file = 'dat_diversity_MAY20.RData')




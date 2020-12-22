library(quantmod)

# load all the data files from KORA.passt (Anja Ludolph)
dataset <- read.csv("/Volumes/GoogleDrive/My\ Drive/DOCTORATE/Thesis/KORA\ DATA/Microbiome_data/KORA_microbiome_variables.csv")

bio_alpha <- read.csv('/Volumes/GoogleDrive/My\ Drive/DOCTORATE/Thesis/KORA\ DATA/Microbiome_data/Speziestabelle\ Bio_Nr/bio_alpha_diversity.csv')

otu_table <- read.csv('/Volumes/GoogleDrive/My\ Drive/DOCTORATE/Thesis/KORA\ DATA/Microbiome_data/Speziestabelle\ Bio_Nr/bio_otu.csv'
                      , check.names=FALSE)

# merge dataset and alpha diversity measures 
dataset2 <- merge(dataset, bio_alpha, by.x = "ff4_prid", by.y = "SampleID")

# transpose otu table
t_otu_table <- setNames(data.frame(t(otu_table[,-1])), otu_table[,1])
t_otu_table$ID <- rownames(t_otu_table)

# merge dataset2 and OTUs
dataset3 <- merge(dataset2, t_otu_table, by.x = "ff4_prid", by.y = "ID")

############################################
# Transform dataset in "analyzable" format #
############################################

first_OTU = grep('OTU', colnames(dataset3))[1]

colnames(dataset3[,1:first_OTU-1])

# KREBS
table(dataset3$lca_sf14) 
table(dataset3$u3c033a)

lca_icd = grep('lca',colnames(dataset3))

summary(dataset3[,lca_icd])

grep('edyrs',colnames(dataset3))

# variables we keep (for the moment)
# 1:3, 6, 8, 10:12, 13, 17:18, 23:24, 25:26, 28:29, 33:34, 36, 38:39, 41, 43, 60:100, 110, 
# 111:133, 134:150, 170, 172, 174, 176
# Krebs, 152:161, 184:200 (to keep as categorical original)

# CONTINUOUS: 3, 11, 13, 25:26, 111:133, 170, 172, 174, 176

# BINARIZE: 6, 18, 29, 33:34, 36, 38:39, 41, 43:44, 60:100, 134:150

cols.num <- c(6, 18, 29, 33:34, 36, 38:39, 41, 43:44, 60:100, 134:150)
# 2 becomes 0; 1 stays 1 >>>> 1 = ja; 2 = nein
dataset3[cols.num] <- sapply(dataset3[cols.num], function(x) as.numeric(x == 1))

# CATEGORIZE: 10, (12), 17, 19, 24, 28, 110

for (i in unique(dataset3$u3tphact)){
  name <- paste("u3tberufb",i,sep ="")
  dataset3[name] <- as.numeric(dataset3$u3tberufb == i)
}

for (i in unique(dataset3$u3tphact)){
  name <- paste("u3tphact",i,sep ="")
  dataset3[name] <- as.numeric(dataset3$u3tphact == i)
}

for (i in unique(dataset3$u3tcigsmk)){
  name <- paste("u3tcigsmk",i,sep ="")
  dataset3[name] <- as.numeric(dataset3$u3tcigsmk == i)
}

for (i in unique(dataset3$u3talkcat)){
  name <- paste("u3talkcat",i,sep ="")
  dataset3[name] <- as.numeric(dataset3$u3talkcat == i)
}

for (i in unique(dataset3$u3twhom)[-length(unique(dataset3$u3twhom))]){
  name <- paste("u3twhom",i,sep ="")
  dataset3[name] <- as.numeric(dataset3$u3twhom == i)
}

for (i in unique(dataset3$u3tsexual)[-length(unique(dataset3$u3tsexual))]){
  name <- paste("u3tsexual",i,sep ="")
  dataset3[name] <- as.numeric(dataset3$u3tsexual == i)
}

# add Date
dataset3$Date <- as.Date(dataset3$u3cuntdat, "%m/%d/%y") 
which(colnames(dataset3) == "Date")

# add year
dataset3$year <- factor(format(dataset3$Date, "%Y"), c("2013","2014"))

# add Seasons
yq <- as.yearqtr(as.yearmon(dataset3$Date) + 1/12)
dataset3$Season <- factor(format(yq, "%q"), levels = 1:4, 
                          labels = c("Winter", "Spring", "Summer", "Fall"))

for (i in unique(dataset3$Season)){
  name <- paste("Season",i,sep ="")
  dataset3[name] <- as.numeric(dataset3$Season == i)
}

grep('Season',colnames(dataset3))

dat_transformed <- dataset3[,c(1:2, 201:207,11:12, 
                               3:5, 8, 13, 23, 25:26, 111:133, 170, 172, 174, 176, # continuous
                               6, 18, 29, 33:34, 36, 38:39, 41, 43:44, 60:100, 134:146, 148:150, 2317:2320, # binary
                               2297:2331, # categorical transformed
                               10, 17, 19, 24, 28, 110, 152:161, 184:200, 147, 2314:2316, # categorical original
                               208:2296)] # OTUs

## CONTINUOUS: u3talteru : u3lk_tria
## BINARY: u3csex : u3tmhp_eradik, SeasonSummer : SeasonFall
## CATEGORICAL: u3tberufb1 : u3tsexual5
## CATEGORICAL (original): u3tphact : u3tsexual, u3c039as5 : u3c039es6 (cancer codes), lca_sf14 : lca_icd5_sf14 (cancer codes), u3tmdbi (drug burden index score)

# save(dat_transformed, file = '/Volumes/GoogleDrive/My\ Drive/DOCTORATE/Thesis/Microbiome\ -\ NOV18/dat_transformed_NOV18.RData')

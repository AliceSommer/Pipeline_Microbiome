
library(sas7bdat)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(igraph)

###############################################################################

# set working directory
setwd('/Users/alicesommer/Desktop/Bureau/DOCTORATE')

# load data formated 
load('data_pipeline_microbiome/dat_transformed_NOV18.RData')

# download long-term exposure data 
dat_pollution = read.sas7bdat('data_pipeline_microbiome/pv_14117g_sommer_gc_20180806.sas7bdat')
head(dat_pollution)

# processing
sum(dat_transformed$ff4_prid %in% dat_pollution$ff4_labid)

# merge 
data_pollution = merge(dat_transformed, dat_pollution, 
                       by.x = 'ff4_prid', by.y = 'ff4_labid', all.y = TRUE)

hist(data_pollution$GC_PM25_14, main = "", xlab = "Long-term PM2.5", breaks = 40)
abline(v=10.3, col="blue")
text(10.7, 160, "10%",
     cex = .8)
abline(v=13, col="blue")
text(13.5, 160, "90%",
     cex = .8)

quantile(data_pollution$GC_PM25_14, prob = seq(0, 1, .10)) 

######################
# Exclusion criteria #
######################

# antiobiotics intake: u3tmabio_j01, u3tmabio_son
dat_excl_temp1 <- data_pollution[(data_pollution$u3tmabio_j01 %in% 0) & (data_pollution$u3tmabio_son %in% 0), ]
# cancer of disgestive organ
dat_excl_temp2 <- dat_excl_temp1[!(dat_excl_temp1$u3c039as6 %in% c('C15','C16','C17','C18','C25','C26') | dat_excl_temp1$lca_icd1_sf14 %in% c('C15','C16','C17','C18','C20','C21','C25')), ]

# prepare the before matching data
data <- dat_excl_temp2
# create the exposure variable
data$W <- NA
data[data$GC_PM25_14 >= 13,]$W <- 0
data[data$GC_PM25_14 <= 10.3,]$W <- 1
data <- subset(data,W!="NA")
dim(data)
table(data$W)

#########################
# Before matching plots #
#########################
source("http://stat.duke.edu/courses/Spring14/sta320.01/CausalInference.R")

colnames_covs = c("u3talteru", "u3tbmi", "u3csex", 
                  "u3tcigsmk2", "u3tcigsmk3", "u3tcigsmk1",
                  "u3tphys", "u3talkkon", "u3tedyrs",
                  #"SeasonSummer", "SeasonSpring", "SeasonWinter", "SeasonFall",
                  "u3tsysmm", "u3tdiamm", "u3lk_chola", "u3lk_hdla", "u3lk_ldla" , "u3lk_tria",
                  "u3tmi", "u3tap", "u3tschl", "u3tca")

balance.plot <- function(d, cov.names = NULL, d2= NULL, main="", left.mar=1.5, right.mar=.1, bottom.mar=.8, top.mar=NULL, xlim=NULL, ...) {
  #d is a vector of differences (probably standardized)
  
  k = length(d)
  
  if(is.null(xlim)) {
    m = max(abs(d))
    xlim=c(-m,m)
  }
  
  
  #if covariate names not given
  if (is.null(cov.names)) {
    cov.names = rep(NA, k)
    for (j in 1:k) cov.names[j] = paste("x",sep="", j)
  }
  
  #setting the top margin, depending on whether there is a title or not
  if (is.null(top.mar)) {
    if (main=="") top.mar=.1
    else top.mar=.8
  }
  cur.mai = par()$mai 
  par(mai=c(bottom.mar, left.mar, top.mar, right.mar))
  
  #creating the plot
  plot(d, k:1, xlim=xlim, pch=19, xlab="Standardized Difference in Covariate Means", axes=FALSE, ylab="", main=main, ...)
  axis(1)
  axis(2, at=k:1, lab = cov.names, las=1, ...)
  abline(v=0)
  abline(v=0.1, lty=3)
  abline(v=-0.1, lty=3)
  
  if (!is.null(d2)) points(d2, k:1, pch=19, col="red")
  
  #reseting margins
  par(mai = cur.mai)
}
cov.balance(data[,colnames_covs], data$W)

### Age ###
g_age <- ggplot(data, aes(x = u3talteru, fill = factor(W)))  +
  geom_density(alpha = .8) + xlab("Age") +
  scale_fill_manual(name = "PM2.5", breaks = c(0,1),
                    labels=c("High","Low"), values = c('gray','green4')) +
  xlim(c(20,95)) +  
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in"))

### BMI ###
colnames(data)[grep('bmi',colnames(data))]
g_bmi <- ggplot(data, aes(x = u3tbmi, fill = factor(W)))  +
  geom_density(alpha = .8) + xlab("BMI (kg/m2)") +
  scale_fill_manual(name = "PM2.5", breaks = c(0,1),
                    labels=c("High","Low"), values = c('gray','green4')) +
  xlim(c(15,60)) +   
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in"))

### Sex ###
g_sex <- ggplot(data, aes(x = factor(W), fill = factor(u3csex))) +
  geom_bar(position = "fill") +
  scale_fill_manual(name = "Sex", breaks = c(0,1),
                    labels=c("Men","Women"), values = c('darkgray','lightgray')) +
  scale_x_discrete(name = "PM2.5", breaks = c(0,1), labels = c("High","Low")) +
  #geom_text(stat = 'count', position = position_fill(vjust = .5),
  #          aes(label = ..count..), angle = 30, size = 3) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in"))

### Smoking ###
colnames(data)[grep('smk', colnames(data))]
table(data$u3tcigsmk)

g_smoke <- ggplot(data, aes(x = factor(W), fill = factor(u3tcigsmk))) + 
  geom_bar(position = "fill")  + 
  scale_fill_brewer(name = "Smoking", palette="RdYlGn", labels = c("Smoker","Ex-Smoker","Never-Smoker")) +
  scale_x_discrete(name = "PM2.5", breaks = c(0,1), labels = c("High","Low")) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) + guides(fill=guide_legend(nrow=3,byrow=TRUE))

### Physical activity ###
g_phys <- ggplot(data, aes(x = factor(W), fill = factor(u3tphys))) + 
  geom_bar(position = "fill")  +
  scale_fill_manual(name = "Physical Acivity", breaks = c(0,1),
                    labels=c("Inactiv","Activ"), values = c('darkgray','lightgray')) +
  scale_x_discrete(name = "PM2.5", breaks = c(0,1), labels = c("High","Low")) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in"))

### Alcohol ###
g_alcohol <- ggplot(data, aes(x = u3talkkon, fill = factor(W)))  +
  geom_density(alpha = .8) + xlab("Alcohol (g/Tag)") +
  scale_fill_manual(name = "PM2.5", breaks = c(0,1),
                    labels=c("High","Low"), values = c('gray','green4')) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in"))

### Education ####
grep('beruf', colnames(data))
g_educ <- ggplot(data, aes(x = factor(W), fill = factor(u3tberufb))) + 
  geom_bar(position = "fill")  + 
  scale_fill_brewer(name = "Education", palette="RdBu", labels = c("non","Professional", "Technical", "Engeneering", "University")) +
  scale_x_discrete(name = "PM2.5", breaks = c(0,1), labels = c("High","Low")) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) + guides(fill=guide_legend(nrow=3,byrow=TRUE))

g_educ_years <- ggplot(data, aes(x = u3tedyrs, fill = factor(W)))  +
  geom_density(alpha = .8) + xlab("Years of education") +
  scale_fill_manual(name = "PM2.5", breaks = c(0,1),
                    labels=c("High","Low"), values = c('gray','green4')) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in"))


### Season ###
data$Season <- factor(data$Season, levels=c("Spring","Summer","Fall","Winter"))

g_season <- ggplot(data, aes(x = factor(W), fill = factor(Season))) + 
  geom_bar(position = "fill")  + 
  scale_fill_brewer(name = "Season", palette="YlGnBu") +
  scale_x_discrete(name = "PM2.5", breaks = c(0,1), labels = c("High","Low")) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) + guides(fill=guide_legend(nrow=2,byrow=TRUE))

g <- grid.arrange(g_age, g_bmi, g_sex, g_alcohol, g_phys, g_smoke, g_educ, g_season,
                  top="Before Matching", ncol = 2)

### Lab variables ###

grep('u3tsysmm', colnames(data)); grep('u3lk_chola',colnames(data)); grep('u3lk_tria',colnames(data)) ## +2 cos no need of WHR
dat_melt_cont_1 = melt(data, id.vars = "W", measure.vars = c(18:19,43:46))

# lvls = levels(as.factor(dat_melt_cont_1$variable))
# # nacounts <- by(dat_melt_cont_1, as.factor(dat_melt_cont_1$variable), function(x) sum(is.na(x$value)))
# # levels(dat_melt_cont_1$variable) = paste(lvls," (NA=",as.integer(nacounts),")",sep="")

levels(dat_melt_cont_1$variable) <- c('sys. BP', 'dia. BP', 'cholesterol', 'HDL chol.', 'LDL chol.', 'triglyceride')

g_lab <- ggplot(dat_melt_cont_1, aes(x=value)) +
  geom_density(aes(group=factor(W), fill=factor(W)),
               alpha = .8) + facet_wrap(~variable, scales = "free") +
  scale_fill_manual(name = "PM2.5", breaks = c(0,1),
                    labels=c("High","Low"), values = c('gray','green4')) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) 

### Disease variables ###

grep('u3tmi', colnames(data)); grep('u3tca',colnames(data))

dat_melt_bin2 = melt(data, id.vars = "W", measure.vars = c(49,52,53,55,56))

levels(dat_melt_bin2$variable) <- c('diabetes', 'myocard. inf.', 'angina p.', 'stroke', 'cancer')

g_disease <- ggplot(dat_melt_bin2, aes(x = factor(W), fill = factor(value))) +
  geom_bar(position = "fill") + facet_wrap(~variable, nrow = 1) +
  scale_fill_manual(name = "Disease", breaks = c(0,1,NA),
                    labels=c("No","Yes","NA"), values = c('darkgray','lightgray','purple'))  +
  scale_x_discrete(name = "PM2.5", breaks = c(0,1), labels = c("High","Low")) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in"),
        axis.text.x = element_text(angle = 30, hjust = 1)) 


### Medication variables ###
grep('u3tmeddia',colnames(data)); grep('u3tmppi',colnames(data))

dat_melt_bin3 = melt(data, id.vars = "W", measure.vars = 58:111)

g_medic <- ggplot(dat_melt_bin3, aes(x = factor(W), fill = factor(value))) +
  geom_bar(position = "fill") + facet_wrap(~variable, nrow = 10) +
  scale_fill_manual(name = "Medication", breaks = c(0,1,NA),
                    labels=c("No","Yes","NA"), values = c('darkgray','lightgray','purple'))  +
  scale_x_discrete(name = "PM2.5", breaks = c(0,1), labels = c("High","Low")) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in"),
        axis.text.x = element_text(angle = 30, hjust = 1)) 


#############
# Matching #
############

setwd("/Volumes/GoogleDrive/My\ Drive/Bureau/Stephane_matching") ### put code elsewhere ### ADD functions to github !
source("matchingOzone_functions_2.R")

data_match = data

data_match$is_treated = as.logical(data_match$W)
data_match$pair_nb = NA

# Optional weights for each covariate when computing the distances
# WARNING: the order of the items in scaling needs to be the same as the order of the covariates (i.e. columns)
scaling =  rep(list(1),ncol(data_match))
names(scaling) = colnames(data_match)
# set the thresholds for each covariate, default is Inf (i.e. no matching)
thresholds = rep(list(Inf),ncol(data_match))
names(thresholds) = colnames(data_match)
# set particular values

thresholds$u3talkkon = 25
thresholds$u3tbmi = 4
thresholds$u3tcigsmk = 0
# thresholds$u3tberufb = 1 # education
# thresholds$u3tedyrs = 3
# thresholds$u3tbmiwho = 0
# #thresholds$u3twhrat = 0 ### waist-hip-ratio
thresholds$u3talteru = 5
thresholds$u3tdiabet = 0
# thresholds$u3tsysmm = 25
# thresholds$u3tdiamm = 20
# thresholds$u3lk_ldla = 50
# thresholds$u3lk_hdla = 25
thresholds$u3csex = 0
thresholds$u3tphys = 0 ### Aktif, Inaktif
# thresholds$u3twhom1 = 0 ### Hypertonieklassifiaktion kategorisch; 1 = Normoton (<140/90 mmHg)
# thresholds$u3tmi = 0 
# thresholds$u3tap = 0 
# thresholds$u3tca = 0 
# thresholds$u3tschl = 0 
# thresholds$Season = 0 

# 
# ## MEDICINE ## 
# thresholds$u3tmbbl = 0 ### Beta-Bockern
# thresholds$u3tmace = 0 ### ACE-Hemmern
# thresholds$u3tmdiu = 0 ### Diuretika
# thresholds$u3tmata = 0 ### Angiotensin-Antagonisten
# thresholds$u3tantihy = 0 ### Antihypertensiva
# thresholds$u3tmhypol = 0 ### Lipidsenkern incl. pflanzliche Stoffe (same variable than u3tlipi)
# thresholds$u3tmstati = 0 ### Statinen
# thresholds$u3tmantco = 0 ### Antikoagulantien
# thresholds$u3tmzop = 0 ### Zopiclon/Zolpidem (Nervensystem Psyche Schlaf)
# thresholds$u3tmibu_b = 0 ### Einnahme von Ibuprofen bzw. Dexibuprofen bei Bedarf, d.h. mindestens eine Dosis innerhalb der letzten 24 Stunden 
# thresholds$u3tmabio_j01 = 0 ### Einnahme von systemischen Antibiotika mit dem ATC= J01
# thresholds$u3tmaadep = 0 ## Einnahme von Andere Antidepressiva 
# thresholds$u3tmbenzo = 0 ## Einnahme von Benzodiazepine als Anxiolytika (ATC= N05BA)
# thresholds$u3tmparacet_b = 0 ### Einnahme von Paracetamol bei Bedarf, d.h. mindestens eine Dosis innerhalb der letzten 24 Stunden 
# # (Präparate mit dem ATC-Code N02AA59, N02AA69, N02AX62, N02BE01, N02BE51, N02BE61 oder R05XA01 und Einnahmemodus: bei Bedarf)

#thresholds$u3tmjodid = 0 # Einnahme von Jodid (ATC= H03AA51 Kombinationen, H03C Iodtherapie)
#thresholds$u3tmmirt = 0 # Einnahme von Mirtazapin
#thresholds$u3tmlith = 0 # Einnahme von Lithium
#thresholds$u3tmlacs = 0 ### very precise medicine
#thresholds$u3tmlc07ab = 0 ### very precise medicine
#thresholds$u3tmlr01b = 0 ### very precise medicine
#thresholds$u3tmlr06a = 0 ### very precise medicine
#thresholds$u3tmace = 5 ### very precise medicine
#thresholds$u3tmcablo = 0 ### Calcium-Antagonisten
#thresholds$u3tmfibra = 0 ### Fibraten
#thresholds$u3tmlr01ad = 0 ### Rhinologika: Corticosteroide zur topischen Anwendung 
#thresholds$u3tsexual7 = 0 ### Andere Sexualhormone (Androgene, Antiandrogene, Enzym-Inhibitoren, Gonadotropin-Releasing-Hormon-Analago und Anabolika)
#thresholds$u3tmppi = 0 ### PPI 
#thresholds$u3tmlacs = 0 # Corticosteroidhaltige Mittel lokal wirksam (ausgenommen D07) --- (ATC = A01AC, A07EA, C05AA, S01BA, S02B, S03C)
#thresholds$u3tmln02a = 0 ### very precise medicine
#thresholds$u3tmln05 = 0 ### very precise medicine
#thresholds$u3tmlr02a = 0  # Hals- und Rachentherapeutika --- (ATC = R02A)


relevant_fields = colnames(data_match)[which(unlist(thresholds)<Inf)]
relevant_fields = c(relevant_fields, "is_treated")

matched_df = data.frame()
total_nb_match = 0
count = 0

start_time = Sys.time()
pb = txtProgressBar(min = 0, max = length(unique(data_match$ff4_prid)), initial = 0, char = "=", style = 3)
count = 0

N = nrow(data_match)
#--------- explore treated units ---------#
treated_units = subset(data_match,is_treated)
control_units = subset(data_match,!is_treated)
N_treated = nrow(treated_units)
# if (N_treated==0){
#   next
# }
N_control = nrow(control_units)
# if (N_control==0){
#   next
# }
cat("Number of treated units:", N_treated,"\nNumber of control units:", N_control,"\n")
#--------------------------------------------------------------------------------------------------------------#
# Compue the discrepancies
discrepancies = discrepancyMatrix(treated_units, control_units, thresholds, scaling)
N_possible_matches = sum(rowSums(discrepancies<Inf)>0)
cat("Number of prospective matched treated units =", N_possible_matches,"\n")
# if (N_possible_matches==0){
#   next
# }
#------------------ Force pair matching via bipartite maximal weighted matching -----------------#
adj = (discrepancies<Inf)
edges_mat = which(adj,arr.ind = TRUE)
weights = 1/(1+sapply(1:nrow(edges_mat),function(i)discrepancies[edges_mat[i,1],edges_mat[i,2]]))
edges_mat[,"col"] = edges_mat[,"col"] + N_treated
edges_vector = c(t(edges_mat))
#-----------------------------------------------------------------------------
# Build the graph from the list of edges
BG = make_bipartite_graph(c(rep(TRUE,N_treated),rep(FALSE,N_control)), edges = edges_vector)
MBM = maximum.bipartite.matching(BG,weights = weights)

# List the dates of the matched pairs
pairs_list = list()
N_matched = 0
for (i in 1:N_treated){
  if (!is.na(MBM$matching[i])){
    N_matched = N_matched + 1
    pairs_list[[N_matched]] = c(i,MBM$matching[i]-N_treated)
  }
}
# Quick sanity check for matched pairs
for (i in 1:N_matched){
  
  total_nb_match = total_nb_match + 1
  # save pair number
  treated_units[pairs_list[[i]][1],"pair_nb"] = total_nb_match
  control_units[pairs_list[[i]][2],"pair_nb"] = total_nb_match
  
  matched_df = rbind(matched_df,treated_units[pairs_list[[i]][1],])
  matched_df = rbind(matched_df,control_units[pairs_list[[i]][2],])
  cat("\n-------------------- Matched pair", total_nb_match,"--------------------\n")
  print(treated_units[pairs_list[[i]][1],relevant_fields])
  print(control_units[pairs_list[[i]][2],relevant_fields])
} 
count = count + 1
setTxtProgressBar(pb,count)
print(Sys.time()-start_time)

table(matched_df$W)

###################
# Diagnostic plot #
###################

### Age ###-
g_age_after <- ggplot(matched_df, aes(x = u3talteru, fill = factor(W)))  +
  geom_density(alpha = .8) + xlab("Age") +
  scale_fill_manual(name = "PM2.5", breaks = c(0,1),
                    labels=c("High","Low"), values = c('gray','green4')) +
  xlim(c(20,95)) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in"))

### BMI ###
colnames(matched_df)[grep('bmi',colnames(matched_df))]
g_bmi_after <- ggplot(matched_df, aes(x = u3tbmi, fill = factor(W)))  +
  geom_density(alpha = .8) + xlab("BMI (kg/m2)") +
  xlim(c(15,60)) +  
  scale_fill_manual(name = "PM2.5", breaks = c(0,1),
                    labels=c("High","Low"), values = c('gray','green4')) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in"))

### Sex ###
g_sex_after <- ggplot(matched_df, aes(x = factor(W), fill = factor(u3csex))) +
  geom_bar(position = "fill") +
  scale_fill_manual(name = "Sex", breaks = c(0,1),
                    labels=c("Men","Women"), values = c('darkgray','lightgray')) +
  scale_x_discrete(name = "PM2.5", breaks = c(0,1), labels = c("High","Low")) +
  #geom_text(stat = 'count', position = position_fill(vjust = .5),
  #          aes(label = ..count..), angle = 30, size = 3) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in"))

### Smoking ###
colnames(matched_df)[grep('smk', colnames(matched_df))]
table(matched_df$u3tcigsmk)

g_smoke_after <- ggplot(matched_df, aes(x = factor(W), fill = factor(u3tcigsmk))) + 
  geom_bar(position = "fill")  + 
  scale_fill_brewer(name = "Smoking", palette="RdYlGn", labels = c("Smoker","Ex-Smoker","Never-Smoker")) +
  scale_x_discrete(name = "PM2.5", breaks = c(0,1), labels = c("High","Low")) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) + guides(fill=guide_legend(nrow=3,byrow=TRUE))

### Physical activity ###
g_phys_after <- ggplot(matched_df, aes(x = factor(W), fill = factor(u3tphys))) + 
  geom_bar(position = "fill")  +
  scale_fill_manual(name = "Physical Acivity", breaks = c(0,1),
                    labels=c("Inactiv","Activ"), values = c('darkgray','lightgray')) +
  scale_x_discrete(name = "PM2.5", breaks = c(0,1), labels = c("High","Low")) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in"))

### Alcohol ###
g_alcohol_after <- ggplot(matched_df, aes(x = u3talkkon, fill = factor(W)))  +
  geom_density(alpha = .8) + xlab("Alcohol (g/Tag)") +
  scale_fill_manual(name = "PM2.5", breaks = c(0,1),
                    labels=c("High","Low"), values = c('gray','green4')) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in"))

### Education ####
grep('beruf', colnames(matched_df))
g_educ_after <- ggplot(matched_df, aes(x = factor(W), fill = factor(u3tberufb))) + 
  geom_bar(position = "fill")  + 
  scale_fill_brewer(name = "Education", palette="RdBu", labels = c("non","Professional", "Technical", "Engineering", "University")) +
  scale_x_discrete(name = "PM2.5", breaks = c(0,1), labels = c("High","Low")) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) + guides(fill=guide_legend(nrow=3,byrow=TRUE))

g_educ_years_after <- ggplot(matched_df, aes(x = u3tedyrs, fill = factor(W)))  +
  geom_density(alpha = .8) + xlab("Years of education") +
  scale_fill_manual(name = "PM2.5", breaks = c(0,1),
                    labels=c("High","Low"), values = c('gray','green4')) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in"))


### Season ###
matched_df$Season <- factor(matched_df$Season, levels=c("Spring","Summer","Fall","Winter"))

g_season_after <- ggplot(matched_df, aes(x = factor(W), fill = factor(Season))) + 
  geom_bar(position = "fill")  + 
  scale_fill_brewer(name = "Season", palette="YlGnBu") +
  scale_x_discrete(name = "PM2.5", breaks = c(0,1), labels = c("High","Low")) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) + guides(fill=guide_legend(nrow=2,byrow=TRUE))

g_after <- grid.arrange(arrangeGrob(g_age, g_bmi, g_sex, g_alcohol, g_phys, g_smoke, g_educ, g_season,
                                    top="Before Matching", ncol = 2),
                        arrangeGrob(g_age_after, g_bmi_after, g_sex_after, g_alcohol_after, 
                                    g_phys_after, g_smoke_after, g_educ_after, g_season_after,
                                    top="After Matching", ncol = 2), ncol = 2)

# ggsave(file = 'plots_pipeline_microbiome/PM25_balance.jpeg',
#        g_after,
#        dpi=300,
#        width = 300,
#        height = 350,
#        units = "mm")

### Lab variables ###

grep('u3tsysmm', colnames(matched_df)); grep('u3lk_chola',colnames(matched_df)); grep('u3lk_tria',colnames(matched_df)) ## +2 cos no need of WHR
dat_melt_cont_1_after = melt(matched_df, id.vars = "W", measure.vars = c(18:19,43:46))

# nacounts <- by(dat_melt_cont_1_after, as.factor(dat_melt_cont_1_after$variable), function(x) sum(is.na(x$value)))
# lvls = levels(as.factor(dat_melt_cont_1_after$variable))
# levels(dat_melt_cont_1_after$variable) = paste(lvls," (NA=",as.integer(nacounts),")",sep="")

levels(dat_melt_cont_1_after$variable) <- c('sys. BP', 'dia. BP', 'cholesterol', 'HDL chol.', 'LDL chol.', 'triglyceride')

g_lab_after <- ggplot(dat_melt_cont_1_after, aes(x=value)) +
  geom_density(aes(group=factor(W), fill=factor(W)),
               alpha = .8) + facet_wrap(~variable, scales = "free") +
  scale_fill_manual(name = "PM2.5", breaks = c(0,1),
                    labels=c("High","Low"), values = c('gray','green4')) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) 

### Disease variables ###

grep('u3tmi', colnames(matched_df)); grep('u3tca',colnames(matched_df))

dat_melt_bin2_after = melt(matched_df, id.vars = "W", measure.vars = c(49,52,53,55,56))

levels(dat_melt_bin2_after$variable) <- c('diabetes', 'myocard. inf.', 'angina p.', 'stroke', 'cancer')

g_disease_after <- ggplot(dat_melt_bin2_after, aes(x = factor(W), fill = factor(value))) +
  geom_bar(position = "fill") + facet_wrap(~variable, nrow = 1) +
  scale_fill_manual(name = "Disease", breaks = c(0,1,NA),
                    labels=c("No","Yes","NA"), values = c('darkgray','lightgray','purple'))  +
  scale_x_discrete(name = "PM2.5", breaks = c(0,1), labels = c("High","Low")) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in"),
        axis.text.x = element_text(angle = 30, hjust = 1)) 

g_after_other <- grid.arrange(arrangeGrob(g_lab, g_disease,
                                          top="Before Matching", nrow = 2),
                              arrangeGrob(g_lab_after, g_disease_after,
                                          top="After Matching", nrow = 2), ncol = 2)

# ggsave(file = 'plots_pipeline_microbiome/PM25_balance_other.jpeg',
#        g_after_other,
#        dpi=300,
#        width = 300,
#        height = 200,
#        units = "mm")

### Medication variables ###
grep('u3tmeddia',colnames(matched_df)); grep('u3tmppi',colnames(matched_df))

dat_melt_bin3_after = melt(matched_df, id.vars = "W", measure.vars = 58:111)

g_medic_after <- ggplot(dat_melt_bin3_after, aes(x = factor(W), fill = factor(value))) +
  geom_bar(position = "fill") + facet_wrap(~variable, nrow = 10) +
  scale_fill_manual(name = "Medication", breaks = c(0,1,NA),
                    labels=c("No","Yes","NA"), values = c('darkgray','lightgray','purple'))  +
  scale_x_discrete(name = "PM2.5", breaks = c(0,1), labels = c("High","Low")) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in"),
        axis.text.x = element_text(angle = 30, hjust = 1)) 

g_after_medic <- grid.arrange(arrangeGrob(g_medic,
                                          top="Before Matching"),
                              arrangeGrob(g_medic_after,
                                          top="After Matching"), ncol = 2)

# ggsave(file = 'plots_pipeline_microbiome/PM25_balance_medic.jpeg',
#        g_after_medic,
#        dpi=300,
#        width = 300,
#        height = 300,
#        units = "mm")


################
# SAVE DATASET #
################
save(matched_df, file = 'data_pipeline_microbiome/dat_matched_PM25_bis.RData')

#######################
#### RANDOMIZATION ####
#######################

n_total = length(matched_df$W)
n_treated  = sum(matched_df$W)
2^n_treated

table(matched_df$W)
head(matched_df$W)
head(matched_df$pair_nb)

# generate a matrix with some possible randomizations (W_sim)
n_col = 10^6
W_sim = matrix(NA, ncol=n_col, nrow=n_total)

for(t in 1:n_col){
  
  W_sim_to_fill = NULL
  flip_coin <- rbinom(n=n_treated,prob=.5,size=1)
  W_sim_to_fill[seq(from = 1, to = n_total - 1, by = 2)] <- flip_coin
  W_sim_to_fill[seq(from = 2, to = n_total, by = 2)] <- 1 - flip_coin
  
  W_sim[,t] = W_sim_to_fill
}

W_unique <- unique(W_sim, MARGIN = 2)
dim(W_unique)

#####################
## SAVE W_balanced ##
#####################

# reorder the data by ff4_prid
W_paired <- W_unique[order(matched_df$ff4_prid),]

# save(W_paired, file = "data_pipeline_microbiome/W_paired_PM25.Rdata")

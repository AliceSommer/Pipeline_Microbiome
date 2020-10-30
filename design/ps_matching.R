
library(sas7bdat)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(Matching)

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
# source("http://stat.duke.edu/courses/Spring14/sta320.01/CausalInference.R")

colnames_covs = c("u3talteru", "u3tbmi", "u3csex", 
                  "u3tcigsmk2", "u3tcigsmk3", "u3tcigsmk1",
                  "u3tphys", "u3talkkon", "u3tedyrs",
                  #"SeasonSummer", "SeasonSpring", "SeasonWinter", "SeasonFall",
                  "u3tsysmm", "u3tdiamm", "u3lk_chola", "u3lk_hdla", "u3lk_ldla" , "u3lk_tria",
                  "u3tmi", "u3tap", "u3tschl", "u3tca")

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

lvls = levels(as.factor(dat_melt_cont_1$variable))
nacounts <- by(dat_melt_cont_1, as.factor(dat_melt_cont_1$variable), function(x) sum(is.na(x$value)))
levels(dat_melt_cont_1$variable) = paste(lvls," (NA=",as.integer(nacounts),")",sep="")

g_lab <- ggplot(dat_melt_cont_1, aes(x=value)) +
  geom_density(aes(group=factor(W), fill=factor(W)),
               alpha = .8) + facet_wrap(~variable, scales = "free") +
  scale_fill_manual(name = "PM2.5", breaks = c(0,1),
                    labels=c("High","Low"), values = c('gray','green4')) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) 

### Disease variables ###

grep('u3tmi', colnames(data)); grep('u3tca',colnames(data))

dat_melt_bin2 = melt(data, id.vars = "W", measure.vars = c(49,52,53,55,56))

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

ggplot(dat_melt_bin3, aes(x = factor(W), fill = factor(value))) +
  geom_bar(position = "fill") + facet_wrap(~variable, nrow = 10) +
  scale_fill_manual(name = "Medication", breaks = c(0,1,NA),
                    labels=c("No","Yes","NA"), values = c('darkgray','lightgray','purple'))  +
  scale_x_discrete(name = "PM2.5", breaks = c(0,1), labels = c("High","Low")) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in"),
        axis.text.x = element_text(angle = 30, hjust = 1)) 

########
## 2. ## Love plot (before matching)
######## 

# remove the subject for whom we don't have alk data
data <- data[!is.na(data$u3talkkon),]

# create a standardized difference in means function (stdif)
stdif <- function(X, W){
  mean_diff = mean(X[W==1], na.rm = TRUE) - mean(X[W==0], na.rm = TRUE)
  SE = sqrt(var(X[W==1], na.rm = TRUE)/sum(W==1, na.rm = TRUE) + var(X[W==0], na.rm = TRUE)/sum(W==0, na.rm = TRUE))
  return(mean_diff/SE)
}

# subset the data to apply the function 
# on the covariates only (dat_covs)
dat_covs <-  data[,colnames_covs]

# apply the function 
# save in stdif_before vector
stdif_before <- sapply(dat_covs, function(x) stdif(x, data$W))
stdif_before

# creating the Love plot
par(mfrow=c(1,1))
plot(stdif_before, 1:length(colnames_covs), axes=FALSE, ylab = "", 
     xlab="T-test statistic", xlim = c(-14,14))
axis(1)
axis(2, at=1:length(colnames_covs), lab = colnames_covs, las=1)
abline(v=0)
abline(v=1.96, lty=3)
abline(v=-1.96, lty=3)

########
## 3. ## Propensity score model
########

# create a propensity score model (ps_model)
?glm
ps_model <- glm(W ~ u3lk_tria + u3tcigsmk3 + u3talkkon + u3talteru, data = data, family="binomial") 

# retrieve the fitted values (ps)
?predict
ps <- predict(ps_model, type="response") 

length(ps)
length(data$W)
########
## 4. ## Overlap assessment and trimming
######## 

# create a plot showing the propensity score distributions of the treated and control units
plot(density(ps[data$W==0], from = min(ps[data$W ==0]), 
             to = max(ps[data$W ==0])), 
     xlim=range(ps), xlab = "Propensity score", col="blue", main="")
lines(density(ps[data$W==1], from = min(ps[data$W ==1]), 
              to = max(ps[data$W ==1])), 
      col="red")
legend("topleft", c("W = 1", "W = 0"), col=c("red", "blue"), lty=c(1,2))

# campare the range 
range(ps[data$W==1])
range(ps[data$W==0])

# eliminate non-comparable cases and create new dataset (dat_elim)
dat_elim = data[ps >= min(ps[data$W==1]) & ps <= max(ps[data$W==0]),] 

dim(dat_elim)
dim(data)

# refit a propensity score model with the trimmed data (ps_model_2)
ps_model_2 <- glm(W ~ u3lk_tria + u3tcigsmk3 + u3talteru, data = dat_elim, family="binomial") 

# retrieve the fitted values (ps_2)
ps_2 <- predict(ps_model_2, type="response") 

# campare the range 
range(ps_2[dat_elim$W==1])
range(ps_2[dat_elim$W==0])

########
## 5. ## Matching
######## 
set.seed(16)
# order data with highest propensity scores first (dat_match)
dat_match = dat_elim[order(ps_2, decreasing = TRUE),]

# order propensity scores with highest propensity scores first (ps_match)
ps_match = ps_2[order(ps_2, decreasing=TRUE)]

head(cbind(ps_match, dat_match))

# basic matching on linear propensity score only, 1:1, without replacement, caliper = 1 (matches)
matches = Match(Tr = dat_match$W, X = ps_match, replace=FALSE, caliper = 1)
# caliper = 1 means that all matches not equal to or within 1 sd 
# of the propensity score are dropped. 

# See ?Match for possible options

# subset data after matching (matched_data)
# hint: matches$index.treated are indices for treated units to keep after matching
matched_data = dat_match[c(matches$index.treated, matches$index.control),]
dim(matched_data)
table(matched_data$W)

# love plot ##

# dataset with only covariates after matching (dat_covs_after)
dat_covs_after <- matched_data[,colnames_covs]

# apply the standardized difference in means function 
# save in stdif_after vector
stdif_after <- sapply(dat_covs_after, function(x) stdif(x, matched_data$W))
stdif_after

# create the Love plot after matching by keeping the before matching information
par(mfrow=c(1,1))
plot(stdif_before, 1:length(colnames_covs), axes=FALSE, ylab = "",
     xlab="T-test statistic", xlim = c(-5,5))
axis(1)
axis(2, at=1:length(colnames_covs), lab = colnames_covs, las=1)
abline(v=0)
abline(v=1.96, lty=3)
abline(v=-1.96, lty=3)
points(stdif_after, 1:length(colnames_covs), pch = 17, col = "red")

###################
# Diagnostic plot #
###################

### Age ###-
g_age_after <- ggplot(matched_data, aes(x = u3talteru, fill = factor(W)))  +
  geom_density(alpha = .8) + xlab("Age") +
  scale_fill_manual(name = "PM2.5", breaks = c(0,1),
                    labels=c("High","Low"), values = c('gray','green4')) +
  xlim(c(20,95)) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in"))

### BMI ###
colnames(matched_data)[grep('bmi',colnames(matched_data))]
g_bmi_after <- ggplot(matched_data, aes(x = u3tbmi, fill = factor(W)))  +
  geom_density(alpha = .8) + xlab("BMI (kg/m2)") +
  xlim(c(15,60)) +  
  scale_fill_manual(name = "PM2.5", breaks = c(0,1),
                    labels=c("High","Low"), values = c('gray','green4')) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in"))

### Sex ###
g_sex_after <- ggplot(matched_data, aes(x = factor(W), fill = factor(u3csex))) +
  geom_bar(position = "fill") +
  scale_fill_manual(name = "Sex", breaks = c(0,1),
                    labels=c("Men","Women"), values = c('darkgray','lightgray')) +
  scale_x_discrete(name = "PM2.5", breaks = c(0,1), labels = c("High","Low")) +
  #geom_text(stat = 'count', position = position_fill(vjust = .5),
  #          aes(label = ..count..), angle = 30, size = 3) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in"))

### Smoking ###
colnames(matched_data)[grep('smk', colnames(matched_data))]
table(matched_data$u3tcigsmk)

g_smoke_after <- ggplot(matched_data, aes(x = factor(W), fill = factor(u3tcigsmk))) + 
  geom_bar(position = "fill")  + 
  scale_fill_brewer(name = "Smoking", palette="RdYlGn", labels = c("Smoker","Ex-Smoker","Never-Smoker")) +
  scale_x_discrete(name = "PM2.5", breaks = c(0,1), labels = c("High","Low")) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) + guides(fill=guide_legend(nrow=3,byrow=TRUE))

### Physical activity ###
g_phys_after <- ggplot(matched_data, aes(x = factor(W), fill = factor(u3tphys))) + 
  geom_bar(position = "fill")  +
  scale_fill_manual(name = "Physical Acivity", breaks = c(0,1),
                    labels=c("Inactiv","Activ"), values = c('darkgray','lightgray')) +
  scale_x_discrete(name = "PM2.5", breaks = c(0,1), labels = c("High","Low")) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in"))

### Alcohol ###
g_alcohol_after <- ggplot(matched_data, aes(x = u3talkkon, fill = factor(W)))  +
  geom_density(alpha = .8) + xlab("Alcohol (g/Tag)") +
  scale_fill_manual(name = "PM2.5", breaks = c(0,1),
                    labels=c("High","Low"), values = c('gray','green4')) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in"))

### Education ####
grep('beruf', colnames(matched_data))
g_educ_after <- ggplot(matched_data, aes(x = factor(W), fill = factor(u3tberufb))) + 
  geom_bar(position = "fill")  + 
  scale_fill_brewer(name = "Education", palette="RdBu", labels = c("non","Professional", "Technical", "Engineering", "University")) +
  scale_x_discrete(name = "PM2.5", breaks = c(0,1), labels = c("High","Low")) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) + guides(fill=guide_legend(nrow=3,byrow=TRUE))

g_educ_years_after <- ggplot(matched_data, aes(x = u3tedyrs, fill = factor(W)))  +
  geom_density(alpha = .8) + xlab("Years of education") +
  scale_fill_manual(name = "PM2.5", breaks = c(0,1),
                    labels=c("High","Low"), values = c('gray','green4')) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in"))


### Season ###
matched_data$Season <- factor(matched_data$Season, levels=c("Spring","Summer","Fall","Winter"))

g_season_after <- ggplot(matched_data, aes(x = factor(W), fill = factor(Season))) + 
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

grep('u3tsysmm', colnames(matched_data)); grep('u3lk_chola',colnames(matched_data)); grep('u3lk_tria',colnames(matched_data)) ## +2 cos no need of WHR
dat_melt_cont_1_after = melt(matched_data, id.vars = "W", measure.vars = c(18:19,43:46))

lvls = levels(as.factor(dat_melt_cont_1_after$variable))
nacounts <- by(dat_melt_cont_1_after, as.factor(dat_melt_cont_1_after$variable), function(x) sum(is.na(x$value)))
levels(dat_melt_cont_1_after$variable) = paste(lvls," (NA=",as.integer(nacounts),")",sep="")

g_lab_after <- ggplot(dat_melt_cont_1_after, aes(x=value)) +
  geom_density(aes(group=factor(W), fill=factor(W)),
               alpha = .8) + facet_wrap(~variable, scales = "free") +
  scale_fill_manual(name = "PM2.5", breaks = c(0,1),
                    labels=c("High","Low"), values = c('gray','green4')) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) 

### Disease variables ###

grep('u3tmi', colnames(matched_data)); grep('u3tca',colnames(matched_data))

dat_melt_bin2_after = melt(matched_data, id.vars = "W", measure.vars = c(49,52,53,55,56))

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

# ggsave(file = 'plots_pipeline_microbiome/M25_balance_other.jpeg',
#        g_after_other,
#        dpi=300,
#        width = 300,
#        height = 200,
#        units = "mm")

### Medication variables ###
grep('u3tmeddia',colnames(matched_data)); grep('u3tmppi',colnames(matched_data))

dat_melt_bin3_after = melt(matched_data, id.vars = "W", measure.vars = 58:111)

ggplot(dat_melt_bin3_after, aes(x = factor(W), fill = factor(value))) +
  geom_bar(position = "fill") + facet_wrap(~variable, nrow = 10) +
  scale_fill_manual(name = "Medication", breaks = c(0,1,NA),
                    labels=c("No","Yes","NA"), values = c('darkgray','lightgray','purple'))  +
  scale_x_discrete(name = "PM2.5", breaks = c(0,1), labels = c("High","Low")) + 
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in"),
        axis.text.x = element_text(angle = 30, hjust = 1)) 

################
# SAVE DATASET #
################
save(matched_data, file = 'data_pipeline_microbiome/dat_matched_PM25_ps.RData')

#######################
#### RANDOMIZATION ####
#######################

n_total = length(matched_data$W)
n_treated  = sum(matched_data$W)
choose(n_total,n_treated)

table(matched_data$W)
head(matched_data$W)

# generate a matrix with some possible randomizations (W_sim)
n_col = 10^6
W_sim = matrix(NA, ncol=n_col, nrow=n_total)

for(t in 1:n_col){
  W_sim[,t] = sample(matched_data$W)
}

# keep only the randomizations for which the balance criteria holds (col_to_keep)
# add the column name in col_to_keep vector
col_to_keep <- NULL

stdif_after_rand <- sapply(dat_covs_after, function(x) stdif(x, W_sim[,1]))
stdif_after_rand 
abs(stdif_after_rand) < 2.2
sum(abs(stdif_after_rand) < 2.2)

for(col in 1:10^6){
  stdif_after_rand <- sapply(dat_covs_after, function(x) stdif(x, W_sim[,col]))
  if (sum(abs(stdif_after_rand) < 2.2) == dim(dat_covs_after)[2]) {
    col_to_keep <- c(col_to_keep,col)
  }
}

10^6 - length(col_to_keep)
head(col_to_keep,20)
tail(col_to_keep,20)

# create a W_balance matrix with only the col_to_keep
W_balanced <- W_sim[,col_to_keep]
dim(W_balanced)

W_unique <- unique(W_balanced, MARGIN = 2)
dim(W_unique)

#####################
## SAVE W_balanced ##
#####################

# reorder the data by ff4_prid
W_ps <- W_unique[order(matched_data$ff4_prid),]

# save(W_ps, file = "data_pipeline_microbiome/W_ps_PM25.Rdata")

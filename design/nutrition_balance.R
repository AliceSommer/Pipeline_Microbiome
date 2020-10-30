library(ggplot2)
library(reshape2)

hi <- read.csv('/Volumes/GoogleDrive/My\ Drive/DOCTORATE/Thesis/KORA\ DATA/Microbiome_data/data_extension_Jan19_nutrition/K14117g_Sommer_tra20181206.csv')
head(hi)

dim(hi)
dat_work <- hi[!is.na(hi$u3v_ena),]

dim(dat_work)

# load('/Users/alicesommer/Desktop/Bureau/DOCTORATE/data_pipeline_microbiome/dat_matched_smoke_bis.RData')
load('/Users/alicesommer/Desktop/Bureau/DOCTORATE/data_pipeline_microbiome/dat_matched_PM25_bis.RData')

head(matched_df$ff4_prid)
head(dat_work$ff4_prid)

table(matched_df$W[matched_df$ff4_prid %in% dat_work$ff4_prid])

# nut_vars <- c('u3v_e07', 'u3v_zb', 'u3v_zk', 'u3v_ze', 'u3v_zf')
# nut_vars <- c('u3v_e02','u3v_e04','u3v_e07', 'u3v_e08')
nut_vars <- paste0(c(rep('u3v_e0', 9), rep('u3v_e', 2)), 1:11)

summary(dat_work[dat_work$ff4_prid %in% matched_df$ff4_prid, nut_vars])

dat_plot <- dat_work[dat_work$ff4_prid %in% matched_df$ff4_prid, c(nut_vars,"ff4_prid")]

dat_plot_mer <- merge(dat_plot, matched_df[,c('ff4_prid','W')], by = "ff4_prid")

dat_melt <- melt(dat_plot_mer, id.vars = "W", measure.vars = nut_vars)

levels(dat_melt$variable) <- c('Potatoes/Roots','Vegetables', "Legumes",'Fruits/Nuts', 'Dairy products',
                               'Cereal products','Meat', 'Fish', 'Egg products', 
                               'Fat', 'Sugar')

# nut_smoke <- ggplot(dat_melt, aes(x=value)) +
#   geom_density(aes(group=factor(W), fill=factor(W)),
#                alpha = .8) + facet_wrap(~variable, scales = "free", ncol = 3)  +
#   scale_fill_manual(name = "Smoking", breaks = c(0,1),
#                     labels=c("Yes (n = 176)","No (n = 220)"), values = c('gray','green4')) + 
#   theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) 

nut_PM <- ggplot(dat_melt, aes(x=value)) +
  geom_density(aes(group=factor(W), fill=factor(W)),
               alpha = .8) + facet_wrap(~variable, scales = "free", ncol = 3) +
  scale_fill_manual(name = "PM2.5", breaks = c(0,1),
                    labels=c("High (n = 65)","Low (n = 78)"), values = c('gray','green4')) +
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in"))

ggsave(file = '/Users/alicesommer/Desktop/Bureau/DOCTORATE/plots_pipeline_microbiome/nut_PM.jpeg',
       nut_PM,
       dpi=300,
       width = 250,
       height = 300,
       units = "mm")

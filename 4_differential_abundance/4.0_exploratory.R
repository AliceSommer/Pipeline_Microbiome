library(compositions)
library(ggplot2)
library(RColorBrewer)
library(reshape2)

## filtered data at 5% prevalence threshold
# ps_work <- readRDS("/Users/alicesommer/Desktop/diff_means_cluster_PM/ps_diff_mean_7_PM.rds")
ps_work <- readRDS("/Users/alicesommer/Desktop/DACOMP_cluster/ps_filt_7_PM.rds") 

x2 = merge_samples(ps_work, "W")
x3 = transform_sample_counts(x2, function(x) x/sum(x))

unname(tax_table(x3)[,"Phylum"])
unname(otu_table(x3)) 

phylcol <- c(brewer.pal(9, "Set1")[1:7], brewer.pal(9, "Greens")[4], brewer.pal(9, "Set1")[8:9])

stack_plot <- plot_bar(x3, fill="Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  # scale_x_discrete(name = "PM2.5", breaks = c(0,1), labels = c("High","Low")) +
  scale_x_discrete(name = "Smoking", breaks = c(0,1), labels = c("Yes","No")) +
  scale_color_manual(values = phylcol) + 
  scale_fill_manual(values = phylcol)

ggsave(stack_plot, file = "/Users/alicesommer/Desktop/Bureau/DOCTORATE/plots_pipeline_microbiome/stack_smoke.png",
       dpi=300,
       width = 85,
       height = 120,
       units = "mm")

## ratio comparison ##

# ps_clr <- transform_sample_counts(ps_work, function(x){x <- x + 1; clr(x)})
# ps_comp <- transform_sample_counts(ps_work, function(x){x <- x + 1; x/sum(x)})
ps_ra <- transform_sample_counts(ps_work, function(x){x/sum(x)})

X <- data.frame(otu_table(ps_ra))
names(X) <- unname(tax_table(ps_ra)[,"Phylum"])

sample_data(ps_ra)$fir_bac <- log(X[,"Firmicutes"]/X[,"Bacteroidetes"])
sample_data(ps_ra)$fir <- X[,"Firmicutes"]
sample_data(ps_ra)$bac <- X[,"Bacteroidetes"]

which(sample_data(ps_ra)$fir_bac > 400)
sample_data(ps_ra)$fir[143]
sample_data(ps_ra)$bac[143]
sample_data(ps_ra)$pro[143]

g_ratio <- ggplot(sample_data(ps_ra), aes(color = factor(W), y = fir_bac)) +
  geom_boxplot(alpha = .5) + ylab('log(Firmicutes/Bacteroidetes)') +
  scale_x_discrete(name = "") +
  scale_colour_manual(values = c("darkgreen","blue4"), limits=c("1","0"), name ="Long-term PM2.5", labels = c("Low","High")) +
  # scale_colour_manual(values = c("darkgreen","blue4"), limits=c("1","0"), name ="Smoking", labels = c("No","Yes")) +
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) +
  guides(color=guide_legend(nrow=2,byrow=TRUE))

ggsave(g_ratio, file = "/Users/alicesommer/Desktop/Bureau/DOCTORATE/plots_pipeline_microbiome/ratio_F_B_PM.png",
       dpi=300,
       width = 85,
       height = 120,
       units = "mm")


# ggplot(sample_data(ps_ra)[!sample_data(ps_ra)$fir_bac > 5,], aes(color = factor(W), y = fir_bac)) +
#   geom_boxplot(alpha = .5) + ylab('Firmicutes/Bacteroidetes') +
#   scale_x_discrete(name = "") +
#   # scale_colour_manual(values = c("darkgreen","blue4"), limits=c("1","0"), name ="Long-term PM2.5", labels = c("Low","High")) +
#   scale_colour_manual(values = c("darkgreen","blue4"), limits=c("1","0"), name ="Smoking", labels = c("No","Yes")) +
#   theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) +
#   guides(color=guide_legend(nrow=2,byrow=TRUE))

# ggplot(sample_data(ps_ra), aes(color = factor(W), y = fir)) +
#   geom_boxplot(alpha = .5) + ylab('Firmicutes') +
#   scale_x_discrete(name = "") +
#   # scale_colour_manual(values = c("darkgreen","blue4"), limits=c("1","0"), name ="Long-term PM2.5", labels = c("Low","High")) +
#   scale_colour_manual(values = c("darkgreen","blue4"), limits=c("1","0"), name ="Smoking", labels = c("No","Yes")) +
#   theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) +
#   guides(color=guide_legend(nrow=2,byrow=TRUE))
# 
# ggplot(sample_data(ps_ra), aes(color = factor(W), y = bac)) +
#   geom_boxplot(alpha = .5) + ylab('Bacteroidetes') +
#   scale_x_discrete(name = "") +
#   # scale_colour_manual(values = c("darkgreen","blue4"), limits=c("1","0"), name ="Long-term PM2.5", labels = c("Low","High")) +
#   scale_colour_manual(values = c("darkgreen","blue4"), limits=c("1","0"), name ="Smoking", labels = c("No","Yes")) +
#   theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) +
#   guides(color=guide_legend(nrow=2,byrow=TRUE))

## Family level
ps_work <- readRDS("/Users/alicesommer/Desktop/diff_means_cluster/ps_diff_mean_4.rds")

ps_ra_Fam <- transform_sample_counts(ps_work, function(x){x/sum(x)})

X_Fam <- data.frame(otu_table(ps_ra_Fam))
names(X_Fam) <- unname(tax_table(ps_ra_Fam)[,"Family"])

sample_data(ps_ra_Fam)$Rum_Bac <- X_Fam[,"Ruminococcaceae"]/X[,"Bacteroidetes"]

ggplot(sample_data(ps_ra_Fam), aes(color = factor(W), y = Rum_Bac)) +
  geom_boxplot(alpha = .5) + ylab('F_Ruminococcaceae/P_Bacteroidetes') +
  scale_x_discrete(name = "") +
  # scale_colour_manual(values = c("darkgreen","blue4"), limits=c("1","0"), name ="Long-term PM2.5", labels = c("Low","High")) +
  scale_colour_manual(values = c("darkgreen","blue4"), limits=c("1","0"), name ="Smoking", labels = c("No","Yes")) +
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) +
  guides(color=guide_legend(nrow=2,byrow=TRUE))

ggplot(sample_data(ps_ra_Fam)[!sample_data(ps_ra_Fam)$Rum_Bac > 100,], aes(color = factor(W), y = Rum_Bac)) +
  geom_boxplot(alpha = .5) + ylab('F_Ruminococcaceae/P_Bacteroidetes') +
  scale_x_discrete(name = "") +
  # scale_colour_manual(values = c("darkgreen","blue4"), limits=c("1","0"), name ="Long-term PM2.5", labels = c("Low","High")) +
  scale_colour_manual(values = c("darkgreen","blue4"), limits=c("1","0"), name ="Smoking", labels = c("No","Yes")) +
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) +
  guides(color=guide_legend(nrow=2,byrow=TRUE))

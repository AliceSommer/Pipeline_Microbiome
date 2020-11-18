
library(VGAM)
library(ggplot2)
library(reshape2)

# setwd('/Users/alicesommer/Desktop/Bureau/DOCTORATE/Pipeline_Microbiome/3_mean_diff_test/Tarray_diff_mean_smoke')
setwd('/Users/alicesommer/Desktop/Bureau/DOCTORATE/Pipeline_Microbiome/3_mean_diff_test/Tarray_diff_mean_PM')

# # smoking
# p_vec <- c(7409, 479, 271,  81, 48, 31, 16)
# PM
p_vec <- c(4370, 414, 252, 74, 44, 29, 15)

Tarray_mat <- matrix(NA, ncol = 7, nrow = 10001)
p_asymp <- matrix(NA, ncol = 7, nrow = 1)

for (d in 1:7){
  # open and concatenate data
  # load(paste0("Tarray_diff_mean", d, ".RData"))
  load(paste0("Tarray_diff_mean", d, "_PM.RData"))
  Tarray_mat[,d] <- Tarray
  
  # randonmized p-values
  clr_x_stat <- Tarray[1]
  p_value <- mean(Tarray >= clr_x_stat)
  
  # asymptotic p-value
  p_clrx <- 1-exp(-1/sqrt(pi)*exp(-(clr_x_stat-(2*log(p_vec[d])-log(log(p_vec[d]))))/2))
  p_asymp[,d] <- p_clrx
  
  # print p-values
  print(paste(Tarray[1], ":", p_value, "-", p_vec[d], p_clrx))
  
  rm("Tarray")
}

# rename columns
colnames(Tarray_mat) <- paste(c('ASV', 'Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum'), '- p =', p_vec)

# randomized p-value for each level
p_vals <- apply(Tarray_mat, 2, function(x) mean(x >= x[1]))

# constant for each level
csts <- sapply(p_vec, function(x) 2*log(x)-log(log(x)))

# matrix of Tarray minus the constant 
m_csts <- matrix(rep(csts,10001), ncol=7, nrow=10001, byrow=T)
Tarray_cst <- Tarray_mat - m_csts

## PLOT ##
# melt data 
Tarray_cst_melt <- melt(Tarray_cst)
colnames(Tarray_cst_melt) <- c('','variable','value')
# save observed test data in data.frame
dat_text_lab <- data.frame(variable = colnames(Tarray_mat))
dat_text_lab$obs_stat <- as.numeric(Tarray_cst[1,])
# save pvalues rand
dat_pvals <- data.frame(variable = colnames(Tarray_mat),
                        p_val = paste('p-val. (rand.) =', round(p_vals,4)))
# save pvalues asymp
dat_p_asymp <- data.frame(variable = colnames(Tarray_mat),
                        p_val = paste('p-val. (asy.) =', as.numeric(round(p_asymp,4))))
# parameters of the Gumbel
beta <- 2
mu <- -log(pi)
# plot
plot_gumbel <- ggplot(data = Tarray_cst_melt, aes(x = value)) + 
  geom_histogram(aes(y = ..density..) , binwidth=.2) +  
  facet_wrap(~variable, ncol = 2) +
  stat_function(fun = VGAM::dgumbel, args = list(location = mu, scale = beta), color = 'blue') +
  geom_vline(data = dat_text_lab, mapping = aes(xintercept = obs_stat), 
             linetype = "dashed", colour = "red", size = .3) +
  geom_text(data = dat_pvals, mapping = aes(x = 10, y = .2, label = p_val), 
          colour = "black", size = 3) + 
  geom_text(data = dat_p_asymp, mapping = aes(x = 13, y = .025, label = p_val), 
            colour = "blue", size = 2.5) + 
  xlab('M − 2log p + log log p')

# ggsave(file = '/Users/alicesommer/Desktop/Bureau/DOCTORATE/plots_pipeline_microbiome/gumbel_diff_means_PM.jpeg', 
#        plot_gumbel,
#        dpi=300,
#        width = 120,
#        height = 180,
#        units = "mm")

# sanity check
test_dist <- rgumbel(90000000, location = mu, scale = beta)

for(i in 1:7){
  # randomized p-value
  print(paste('random:',mean(Tarray_mat[,i] - csts[i] >= Tarray_mat[1,i] - csts[i]),
              # p-value based on generated Gumbel with assumed parameters mu and beta
              'asymp:',mean(test_dist >= Tarray_mat[1,i] - csts[i])))
}

# #######################
# # Multiple comparison #
# #######################
# 
# nrep = dim(Tarray_mat)[1]
# ntest = dim(Tarray_mat)[2]
# 
# hyp_matrix <- Tarray_mat
# hyp_p_value <- matrix(NA, ncol = ntest, nrow = nrep)
# 
# verbose = T
# 
# # based on value (hyp_obs) of each row
# for (r in 1:nrep){
#   if(verbose)
#     if(r%% ceiling(nrep/100) == 1)
#       cat(paste0('Testing rep : ',r,'/',nrep,' \n\r'))
#   # calc. hypothetical p_value on each column of the matrix 
#   hyp_p_value[r,] <- apply(hyp_matrix, 2, function(x) mean(x >= x[r]))
# }
# 
# # for each rep. take the min. p_value
# min_p_nrep <- apply(hyp_p_value, 1, function(x) min(x, na.rm = TRUE))
# head(min_p_nrep,20)
# 
# # calculate the proportion of min_p_nrep that is sm/eq. p_value (for obs.)
# p_value_adj <- sapply(p_vals, function(x) mean(min_p_nrep <= x))
# p_value_adj

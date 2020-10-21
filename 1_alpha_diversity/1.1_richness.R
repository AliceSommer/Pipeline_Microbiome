library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(breakaway)
library(dplyr)
library(gridExtra)

###############################################################################

# set working directory
setwd('/Users/alicesommer/Desktop/Bureau/DOCTORATE/data_pipeline_microbiome')

# load microbiome data
ASV_table <- readRDS('dada2output/seqtab2020.rds')
taxon_assign <- readRDS('dada2output/taxa2020.rds')

# load sample/matched_data
# load('dat_matched_PM25_bis.RData')
load('dat_matched_smoke_bis.RData')

# load W matrix for randomization test
# load("W_paired_PM25.Rdata")
load("W_paired_smoke_bis.Rdata")

# matched AP data
sample_df <- matched_df[order(matched_df$ff4_prid),]
sample_df$W <- as.factor(sample_df$W)
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

### 1. ESTIMATE THE TOTAL RICHNESS FOR EACH SAMPLE ###
rich <- sample_richness(ps_prune)
ba <- breakaway(ps_prune)
ba[[1]]
# plot(ba, ps_prune, color = "W")

sample_data(ps_prune)[,"breakaway_W"] <- summary(ba)$estimate
sample_data(ps_prune)[,"ba_error_W"] <- summary(ba)$error
sample_data(ps_prune)[,"lower"] <- summary(ba)$lower
sample_data(ps_prune)[,"upper"] <- summary(ba)$upper
# upper is NA when SE so small that confidence interval is the estimate...

# add sequencing depth as variable
sample_data(ps_prune)[,"sample_counts"] <- sample_sums(ps_prune)
sample_data(ps_prune)[,"richness"] <- summary(rich)$estimate

sample_data(ps_prune)[,"dim"] <- as.factor(1:dim(sample_df)[1])

ggplot(sample_data(ps_prune), aes(x = dim)) + 
  geom_point(aes(y = breakaway_W), colour = "red") + 
  # geom_errorbar(aes(ymin=lower, ymax=upper)) +
  geom_point(aes(y = richness), colour = "blue", alpha = .5) +
  xlab('Samples') + ylab("BA estimate")

hist(sample_data(ps_prune)$ba_error_W, breaks = 50)
head(sample_data(ps_prune)[sample_data(ps_prune)$breakaway_W > 400, c("ba_error_W","sample_counts","breakaway_W","pair_nb")])

summary(sample_data(ps_prune)$sample_counts)

# remove the outlier in seq. depth ?
ps_prune_out <- subset_samples(ps_prune, pair_nb != 62)
  
###############
## BREAKAWAY ##
###############

g_PM <- ggplot(sample_data(ps_prune_out), aes(color = factor(W), y = breakaway_W)) +
  geom_boxplot(alpha = .5) + ylab('breakaway richness measure') +
  scale_x_discrete(name = "") +
  # scale_colour_manual(values = c("darkgreen","blue4"), limits=c("1","0"), name ="Long-term PM2.5", labels = c("Low","High")) +
  scale_colour_manual(values = c("darkgreen","blue4"), limits=c("1","0"), name ="Smoking", labels = c("No","Yes")) +
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) +
  guides(color=guide_legend(nrow=2,byrow=TRUE)) + 
  annotate(geom="text",x=.7, y=89, label="stat. = 3.9576") +
  annotate(geom="text",x=.7, y=100, label="p-value = 0.1518") 

# PM: p-value = 0.0008; test-statistic = 19.0775
# Smoke: p-value = 0.1518; test-statistic = 3.9576

# ggsave(g_PM, file = "/Users/alicesommer/Desktop/Bureau/DOCTORATE/plots_pipeline_microbiome/box_break_smoke.png",
#        dpi=300,
#        width = 85,
#        height = 120,
#        units = "mm")

g_PM_sex <- ggplot(sample_data(ps_prune), aes(x = factor(u3csex), y = breakaway_W, color = factor(W))) +
  geom_boxplot(alpha = .5) + ylab('breakaway richness measure') +
  scale_x_discrete(name = "Sex", limits=c("0","1"), labels = c("Female", "Male")) +
  scale_colour_manual(values = c("gray","blue4"), limits=c("1","0"), name ="Long-term PM2.5", labels = c("Low","High")) +
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) +
  guides(color=guide_legend(nrow=2,byrow=TRUE))

g_arrange <- grid.arrange(g_PM,g_PM_sex, nrow = 1)

### 2. USE BETTA FUNCTION ###
x <- cbind(1, 
           # sample_data(ps_prune)$u3tcigsmk,
           # sample_data(ps_prune)$u3csex
          sample_data(ps_prune_out)$W)
head(sample_data(ps_prune_out)$breakaway_W)
head(sample_data(ps_prune_out)$ba_error_W)

reg <- betta(sample_data(ps_prune_out)$breakaway_W,
             sample_data(ps_prune_out)$ba_error_W, X = x)
reg$table
estim_obs <- reg$table[2,1]

### 3. PERFORM A RANDOMIZATION TEST ###
dim(W_paired_smoke)
W_paired_smoke <- W_paired_smoke[-which(sample_data(ps_prune)$pair_nb == 62), ]
# dim(W_paired)

# set the number of randomizations
nrep <- ncol(W_paired_smoke)/100
# nrep <- ncol(W_paired)/100

# create a matrix where the t_rand will be saved
t_array <- NULL

for(i in 1:nrep){
  print(i)
  x = cbind(1, 
            # sample_data(ps_prune)$u3tcigsmk, 
            # sample_data(ps_prune)$u3csex
            W_paired_smoke[,i])
            # W_paired[,i]) 
            
  
  reg = betta(sample_data(ps_prune_out)$breakaway_W,
              sample_data(ps_prune_out)$ba_error_W, X = x)
  
  # fill t_array
  t_array[i] = reg$table[2,1] 
}

## calculate p_value
p_value <- mean(t_array >= estim_obs, na.rm = TRUE)
p_value

# 1. Open png file
png("/Users/alicesommer/Desktop/Bureau/DOCTORATE/plots_pipeline_microbiome/beta_break_smoke.png", width = 350, height = 350)
# 2. Create the plot
hist(t_array, breaks = 30, main = "", xlab = "breakaway beta (Smoking)")
abline(v = estim_obs, col = 'red', lwd = 2, lty = 2)
# 3. Close the file
dev.off()

# ##############################
# ### betta (tweeked solver) ###
# ##############################
# betta_tweek <- function(chats, ses, X){
#   if (isTRUE(is.na(X))) {
#     X <- matrix(rep(1, length(chats)), ncol = 1)
#   }
#   consider <- !(is.na(chats) | is.na(ses) | apply(is.na(X),
#                                                   1, sum))
#   chats_effective <- chats[consider]
#   ses_effective <- ses[consider]
#   X_effective <- as.matrix(X[consider, ])
#   n <- dim(X_effective)[1]
#   p <- dim(X_effective)[2]
#   likelihood <- function(input) {
#     ssq_u <- input[1]
#     beta <- input[2:length(input)]
#     W <- diag(1/(ssq_u + ses_effective^2))
#     -0.5 * (sum(log(ssq_u + ses_effective^2) + (chats_effective -
#                                                   X_effective %*% beta)^2/(ssq_u + ses_effective^2)) +
#               log(det(t(X_effective) %*% W %*% X_effective)))
#   }
#   if (any(is.null(initial_est))) {
#     initial_est <- c(var(chats_effective), solve(t(X_effective) %*%
#                                                    X_effective) %*% t(X_effective) %*% chats_effective)
#   }
#   output <- try(optim(initial_est, likelihood, hessian = FALSE,
#                       control = list(fnscale = -1), lower = c(0, rep(-Inf,
#                                                                      p)), method = "L-BFGS-B"), silent = TRUE)
#   i = 0
#   while ("try-error" %in% class(output) & i < 200) {
#     i <- i + 1
#     perturb <- rnorm(n = length(initial_est), mean = c(0, 0), sd = 0.001 * i * abs(initial_est))
#     initial_est_perturbed <- pmax(c(0, rep(-Inf, p)), initial_est + perturb)
#     output <- try(optim(initial_est_perturbed, likelihood,
#                         hessian = FALSE, control = list(fnscale = -1), lower = c(0,
#                                                                                  rep(-Inf, p)), method = "L-BFGS-B"), silent = TRUE)
#   }
#   if ("try-error" %in% class(output)) {
#     stop(paste("The starting value and 200 perturbations were not",
#                "enough to find a maximum likelihood solution.",
#                "Please try again with a new choice of `initial_est`."))
#   }
#   
#   ssq_u <- output$par[1]
#   beta <- output$par[2:length(output$par)]
#   
#   W <- diag(1/(ssq_u + ses_effective^2))
#   vars <- 1/diag(t(X_effective) %*% W %*% X_effective)
#   global <- t(beta) %*% (t(X_effective) %*% W %*% X_effective) %*%
#     beta
#   Q <- sum((chats_effective - X_effective %*% beta)^2/ses_effective^2)
#   R <- diag(ses_effective^2)
#   G <- diag(ssq_u, n)
#   mytable <- list()
#   mytable$table <- cbind(Estimates = beta, `Standard Errors` = sqrt(vars),
#                          `p-values` = round(2 * (1 - pnorm(abs(beta/sqrt(vars)))), 3))
#   mytable$cov <- solve(t(X_effective) %*% W %*% X_effective)
#   mytable$ssq_u <- ssq_u
#   mytable$homogeneity <- c(Q, 1 - pchisq(Q, n - p))
#   mytable$global <- c(global, 1 - pchisq(global, p - 1))
#   us <- c(ssq_u * W %*% (chats_effective - X_effective %*%
#                            beta))
#   blups <- rep(NA, length(chats))
#   blups[consider] <- c(X_effective %*% beta + us)
#   mytable$blups <- blups
#   var_matrix <- matrix(NA, nrow = (n + p), ncol = (n + p))
#   var_matrix[1:p, 1:p] <- t(X_effective) %*% solve(R, tol = 1e-20) %*% X_effective
#   var_matrix[(p + 1):(n + p), (p + 1):(n + p)] <- solve(R, tol = 1e-20) +
#     MASS::ginv(G)
#   var_matrix[1:p, (p + 1):(n + p)] <- t(X_effective) %*% solve(R, tol = 1e-20)
#   var_matrix[(p + 1):(n + p), 1:p] <- solve(R, tol = 1e-20) %*% X_effective
#   var_matrix_inv <- MASS::ginv(var_matrix)
#   blupvars <- rep(NA, length(chats))
#   blupvars[consider] <- (cbind(X_effective, diag(1, n)) %*% var_matrix_inv %*% t(cbind(X_effective, diag(1, n)))) %>% diag %>% sqrt %>% c
#   mytable$blupses <- blupvars
#   logLhat <- -0.5 * (n * log(2 * pi) + sum(log(ssq_u + ses_effective^2) +
#                                              (chats_effective - X_effective %*% beta)^2/(ssq_u + ses_effective^2)))
#   mytable$loglikelihood <- logLhat
#   mytable$aic <- -2 * logLhat + 2 * (1 + p)
#   mytable$aicc <- mytable$aic + (2 * (1 + p)^2 + 2 * (1 + p))/(n -
#                                                                  (1 + p) - 1)
#   mytable$r_squared_wls <- 1 - sum((chats_effective - X_effective %*%
#                                       beta)^2)/(sum(chats_effective^2) - n * (mean(chats_effective))^2)
#   
#   return(mytable)
# }
# 
# chats = sample_data(ps_prune)$breakaway_W
# ses = sample_data(ps_prune)$ba_error_W
# X = x
# initial_est = NULL
# 
# betta_tweek(chats,ses,X)$table
# 
# ### PLOT PAIRED DATA ###
# 
# sample_data(ps_prune)$alpha_ci_low <- summary(ba)$lower
# sample_data(ps_prune)$alpha_ci_high <- summary(ba)$upper
# 
# dat_plot = sample_data(ps_prune)
# 
# # order plot
# dat_plot$pair_nb <- factor(dat_plot$pair_nb, levels = dat_plot$pair_nb[dat_plot$W == 1][order(dat_plot$breakaway_W[dat_plot$W == 1])])
# dat_plot$pair_nb  # notice the changed order of factor levels
# 
# ggplot(dat_plot, 
#        aes(color = factor(W), y = breakaway_W, x = pair_nb)) +
#   geom_point() + 
#   geom_errorbar(aes(ymin=alpha_ci_low, ymax=alpha_ci_high), width=.1) +
#   ylab('Breakaway richness estimate') +
#   scale_x_discrete(name = "") +
#   scale_colour_manual(values = c("gray","blue4"), limits=c("1","0"), 
#                       name ="Long-term PM2.5", 
#                       labels = c("Low","High")) +
#   theme(axis.text.x=element_text(angle =- 70, vjust = 0.5)) + ylim(c(30,300))
# 
# g_ba <- ggplot(dat_plot, aes(color = factor(W), y = breakaway_W)) +
#   geom_boxplot(alpha = .5) + ylab('Breakaway richness estimate') +
#   scale_x_discrete(name = "") +
#   scale_colour_manual(values = c("gray","blue4"), limits=c("1","0"), 
#                       name ="Long-term PM2.5", 
#                       labels = c("Low (<= 11)","High (>= 12)")) +
#   theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) +
#   guides(color=guide_legend(nrow=2,byrow=TRUE))

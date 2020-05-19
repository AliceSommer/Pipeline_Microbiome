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
load('dat_matched_PM25_bis.RData')

# load W matrix for randomization test
load("W_paired_PM25.Rdata")

# load original dataset (try other Ws)
# df <- read.csv('/Volumes/GoogleDrive/My\ Drive/DOCTORATE/Thesis/KORA\ DATA/Microbiome_data/KORA_microbiome_variables.csv')

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

### 1. ESTIMATE THE TOTAL DIVERSITY FOR EACH SAMPLE ###
rich <- sample_richness(ps_prune)
ba <- breakaway(ps_prune)
# ba <- breakaway_nof1(ps_prune)
ba[[1]]
# plot(ba, ps_prune, color = "W") 

sample_data(ps_prune)[,"breakaway"] <- summary(ba)$estimate

g_PM <- ggplot(sample_data(ps_prune), aes(color = factor(W), y = breakaway)) +
  geom_boxplot(alpha = .5) + ylab('breakaway richness measure') +
  scale_x_discrete(name = "") +
  scale_colour_manual(values = c("gray","blue4"), limits=c("1","0"), name ="Long-term PM2.5", labels = c("Low","High")) +
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) +
  guides(color=guide_legend(nrow=2,byrow=TRUE))

g_PM_sex <- ggplot(sample_data(ps_prune), aes(x = factor(u3csex), y = breakaway, color = factor(W))) +
  geom_boxplot(alpha = .5) + ylab('breakaway richness measure') +
  scale_x_discrete(name = "Sex", limits=c("0","1"), labels = c("Female", "Male")) +
  scale_colour_manual(values = c("gray","blue4"), limits=c("1","0"), name ="Long-term PM2.5", labels = c("Low","High")) +
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) +
  guides(color=guide_legend(nrow=2,byrow=TRUE))

g_arrange <- grid.arrange(g_PM,g_PM_sex, nrow = 1)

x <- cbind(1, sample_data(ps_prune)$W, sample_data(ps_prune)$u3csex, sample_data(ps_prune)$u3tcigsmk1)

head(summary(ba)$estimate)
head(summary(ba)$error)
  
### 2. USE BETTA FUNCTION ###
reg <- betta(summary(rich)$estimate,
             summary(ba)$error, X = x)
reg$table
estim_obs <- reg$table[2,1]

### 3. PERFORM A RANDOMIZATION TEST ###
dim(W_paired)

# set the number of randomizations
nrep <- ncol(W_paired)/100

# create a matrix where the t_rand will be saved
t_array <- NULL

for(i in 1:nrep){
  print(i)
  x = cbind(1, W_paired[,i], 
             sample_data(ps_prune)$u3csex, 
             sample_data(ps_prune)$u3tcigsmk1)
  
  reg = betta(summary(rich)$estimate,
               summary(ba)$error, X = x)
  
  # fill t_array
  t_array[i] = reg$table[2,1] 
}

## calculate p_value
p_value <- mean(t_array >= estim_obs)
p_value
hist(t_array, breaks = 30)

##############################
### betta (tweeked solver) ###
##############################
data_check <- data.frame(summary(ba)$estimate, summary(ba)$error, W = sample_data(ps)$W)
head(data_check)

sum(duplicated(data_check))

chats = data_check[,1]
ses = data_check[,2]
X = cbind(1, data_check[,3])
initial_est = NULL

  if (isTRUE(is.na(X))) {
    X <- matrix(rep(1, length(chats)), ncol = 1)
  }
  consider <- !(is.na(chats) | is.na(ses) | apply(is.na(X),
                                                  1, sum))
  chats_effective <- chats[consider]
  ses_effective <- ses[consider]
  X_effective <- as.matrix(X[consider, ])
  n <- dim(X_effective)[1]
  p <- dim(X_effective)[2]
  likelihood <- function(input) {
    ssq_u <- input[1]
    beta <- input[2:length(input)]
    W <- diag(1/(ssq_u + ses_effective^2))
    -0.5 * (sum(log(ssq_u + ses_effective^2) + (chats_effective -
                                                  X_effective %*% beta)^2/(ssq_u + ses_effective^2)) +
              log(det(t(X_effective) %*% W %*% X_effective)))
  }
  if (any(is.null(initial_est))) {
    initial_est <- c(var(chats_effective), solve(t(X_effective) %*%
                                                   X_effective) %*% t(X_effective) %*% chats_effective)
  }
  output <- try(optim(initial_est, likelihood, hessian = FALSE,
                      control = list(fnscale = -1), lower = c(0, rep(-Inf,
                                                                     p)), method = "L-BFGS-B"), silent = TRUE)
  i = 0
  while ("try-error" %in% class(output) & i < 200) {
    i <- i + 1
    perturb <- rnorm(n = length(initial_est), mean = c(0, 0), sd = 0.001 * i * abs(initial_est))
    initial_est_perturbed <- pmax(c(0, rep(-Inf, p)), initial_est + perturb)
    output <- try(optim(initial_est_perturbed, likelihood,
                        hessian = FALSE, control = list(fnscale = -1), lower = c(0,
                                                                                 rep(-Inf, p)), method = "L-BFGS-B"), silent = TRUE)
  }
  if ("try-error" %in% class(output)) {
    stop(paste("The starting value and 200 perturbations were not",
               "enough to find a maximum likelihood solution.",
               "Please try again with a new choice of `initial_est`."))
  }

  ssq_u <- output$par[1]
  beta <- output$par[2:length(output$par)]

  W <- diag(1/(ssq_u + ses_effective^2))
  vars <- 1/diag(t(X_effective) %*% W %*% X_effective)
  global <- t(beta) %*% (t(X_effective) %*% W %*% X_effective) %*%
    beta
  Q <- sum((chats_effective - X_effective %*% beta)^2/ses_effective^2)
  R <- diag(ses_effective^2)
  G <- diag(ssq_u, n)
  mytable <- list()
  mytable$table <- cbind(Estimates = beta, `Standard Errors` = sqrt(vars),
                         `p-values` = round(2 * (1 - pnorm(abs(beta/sqrt(vars)))), 3))
  mytable$cov <- solve(t(X_effective) %*% W %*% X_effective)
  mytable$ssq_u <- ssq_u
  mytable$homogeneity <- c(Q, 1 - pchisq(Q, n - p))
  mytable$global <- c(global, 1 - pchisq(global, p - 1))
  us <- c(ssq_u * W %*% (chats_effective - X_effective %*%
                           beta))
  blups <- rep(NA, length(chats))
  blups[consider] <- c(X_effective %*% beta + us)
  mytable$blups <- blups
  var_matrix <- matrix(NA, nrow = (n + p), ncol = (n + p))
  var_matrix[1:p, 1:p] <- t(X_effective) %*% solve(R, tol = 1e-20) %*% X_effective
  var_matrix[(p + 1):(n + p), (p + 1):(n + p)] <- solve(R, tol = 1e-20) +
    MASS::ginv(G)
  var_matrix[1:p, (p + 1):(n + p)] <- t(X_effective) %*% solve(R, tol = 1e-20)
  var_matrix[(p + 1):(n + p), 1:p] <- solve(R, tol = 1e-20) %*% X_effective
  var_matrix_inv <- MASS::ginv(var_matrix)
  blupvars <- rep(NA, length(chats))
  blupvars[consider] <- (cbind(X_effective, diag(1, n)) %*% var_matrix_inv %*% t(cbind(X_effective, diag(1, n)))) %>% diag %>% sqrt %>% c
  mytable$blupses <- blupvars
  logLhat <- -0.5 * (n * log(2 * pi) + sum(log(ssq_u + ses_effective^2) +
                                             (chats_effective - X_effective %*% beta)^2/(ssq_u + ses_effective^2)))
  mytable$loglikelihood <- logLhat
  mytable$aic <- -2 * logLhat + 2 * (1 + p)
  mytable$aicc <- mytable$aic + (2 * (1 + p)^2 + 2 * (1 + p))/(n -
                                                                 (1 + p) - 1)
  mytable$r_squared_wls <- 1 - sum((chats_effective - X_effective %*%
                                      beta)^2)/(sum(chats_effective^2) - n * (mean(chats_effective))^2)

mytable

sample_data(ps_prune)$ba_estimate <- summary(ba)$estimate
sample_data(ps_prune)$alpha_ci_low <- summary(ba)$lower
sample_data(ps_prune)$alpha_ci_high <- summary(ba)$upper

dat_plot = sample_data(ps_prune)

# order plot
dat_plot$pair_nb <- factor(dat_plot$pair_nb, levels = dat_plot$pair_nb[dat_plot$W == 1][order(dat_plot$ba_estimate[dat_plot$W == 1])])
dat_plot$pair_nb  # notice the changed order of factor levels

ggplot(dat_plot, 
       aes(color = factor(W), y = ba_estimate, x = pair_nb)) +
  geom_point() + 
  geom_errorbar(aes(ymin=alpha_ci_low, ymax=alpha_ci_high), width=.1) +
  ylab('Breakaway richness estimate') +
  scale_x_discrete(name = "") +
  scale_colour_manual(values = c("gray","blue4"), limits=c("1","0"), 
                      name ="Long-term PM2.5", 
                      labels = c("Low (<= 11)","High (>= 12)")) +
  theme(axis.text.x=element_text(angle =- 70, vjust = 0.5)) + ylim(c(30,300))

g_ba <- ggplot(dat_plot, aes(color = factor(W), y = ba_estimate)) +
  geom_boxplot(alpha = .5) + ylab('Breakaway richness estimate') +
  scale_x_discrete(name = "") +
  scale_colour_manual(values = c("gray","blue4"), limits=c("1","0"), 
                      name ="Long-term PM2.5", 
                      labels = c("Low (<= 11)","High (>= 12)")) +
  theme(legend.position = "top", legend.key.size =  unit(0.1, "in")) +
  guides(color=guide_legend(nrow=2,byrow=TRUE))

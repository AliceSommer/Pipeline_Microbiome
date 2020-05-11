library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(breakaway)
library(dplyr)

###############################################################################

# set working directory
setwd('/Users/alicesommer/Desktop/Bureau/DOCTORATE/data_pipeline_microbiome')

# load microbiome data
ASV_table <- readRDS('dada2output/seqtab2020.rds')
taxon_assign <- readRDS('dada2output/taxa2020.rds')

# load sample/matched_data
load('dat_matched_PM25.RData')

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

### 1. ESTIMATE THE TOTAL DIVERSITY FOR EACH SAMPLE ###
ba <- breakaway(ps_prune)
# ba <- breakaway_nof1(ps)
ba[[1]]
plot(ba, ps_prune, color = "W") 

x = cbind(1, sample_data(ps)$W)

## retrieve estimates and ses of breakaway
break_estimates <- NULL
break_ses <- NULL

for (i in 1:dim(sample_data(ps))[1]) {
  est = ba[[i]]
  break_estimates[i] <- est$estimate
  break_ses[i] <- est$error
}

head(break_estimates)
head(break_ses)

### 2. USE BETTA FUNCTION ###
summary(lm(break_estimates ~ sample_data(ps)$W))

reg <- betta(break_estimates, break_ses, X = x)
reg

head(cbind(break_estimates, break_ses, sample_data(ps)$W), u3tbmi, u3talteru)

data_check <- data.frame(break_estimates, break_ses, W = sample_data(ps)$W)
data_X_W <- data_check[!(data_check[,2] < 0.43),]
head(data_X_W)
dim(data_X_W)

# there are two samples that are not unique !

# run regression when removes
reg_2 <- betta(data_X_W[,1], data_X_W[,2], 
               # X = cbind(1, data_check[,3], data_check[,4], data_check[,5]))
               X = cbind(1, data_X_W[,3]))
reg_2


### betta (tweeked solver) ###

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

# # try the objective bayes method
# frq_table <- build_frequency_count_tables(otu_table(ps))
# 
# dim(table(frq_table[[4]]))
# 
# length_freq <- NULL
# sample_id <- NULL
# 
# for(i in 1:486){
#   sample_id <- c(sample_id,i)
#   length_freq <- c(length_freq,length(frq_table[[i]]$Freq))
# }
# 
# breakaway(frq_table[[1]])
# 
# # retrieve estimates and ses of ob_nb
# ob_nb_estimates <- NULL
# ob_nb_ses <- NULL
# 
# for (i in 1:dim(sample_data(ps))[1]) {
#   ob_nb <- objective_bayes_negbin(frq_table[[i]], output = F, plot = F, answers = T)
#   ob_nb_estimates[i] <- ob_nb$est
#   ob_nb_ses[i] <- ob_nb$semeanest
# }
# 
# head(ob_nb_estimates)
# head(ob_nb_ses)
# 
# x = cbind(1, sample_data(ps)$W)
# 
# reg <- betta(ob_nb_estimates, ob_nb_ses, X = x)
# reg

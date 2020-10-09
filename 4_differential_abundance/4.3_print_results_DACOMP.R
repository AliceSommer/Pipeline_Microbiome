library(phyloseq)

setwd('/Users/alicesommer/Desktop/DACOMP_cluster/')
rank_names <- c( "ASV" , "Species" , "Genus"  , "Family" , "Order"  , "Class"  , "Phylum")

###########
# SMOKING #
###########

load('dacomp_results/dacomp_results_smoke.RData')

for (d in 1:7){
  condition <- dacomp_results$p_value_adj[dacomp_results$rank == rank_names[d]] <= 0.2
  print(paste(d, ":", sum(condition, na.rm = TRUE)))
  print(dacomp_results$p_value_adj[dacomp_results$rank == rank_names[d]][which(condition)])
  print(dacomp_results$raw_test_stat[dacomp_results$rank == rank_names[d]][which(condition)])
  ps <- readRDS(paste0("ps_filt_", d, ".rds"))
  print(unname(tax_table(ps)[which(condition),]))
}

#############
# POLLUTION #
#############

load('dacomp_results/dacomp_results_PM.RData')

for (d in 1:7){
  condition <- dacomp_results$p_value_adj[dacomp_results$rank == rank_names[d]] <= 0.2
  print(paste(d, ":", sum(condition, na.rm = TRUE)))
  print(dacomp_results$p_value_adj[dacomp_results$rank == rank_names[d]][which(condition)])
  print(dacomp_results$raw_test_stat[dacomp_results$rank == rank_names[d]][which(condition)])
  ps <- readRDS(paste0("ps_filt_", d, "_PM.rds"))
  if(sum(condition > 0, na.rm = TRUE)){print(unname(tax_table(ps)[which(condition),]))}
}

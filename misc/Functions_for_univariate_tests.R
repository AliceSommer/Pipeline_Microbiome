#Libraries used
library(Hmisc)
library(coin)
library(dacomp)
dt_inst = installed.packages()
if("doRNG" %in% rownames(dt_inst))
  library(doRNG)
if("parallel" %in% rownames(dt_inst))
  library(parallel)
if("doParallel" %in% rownames(dt_inst))
  library(doParallel)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Definitions for early stopping of permutations:
# Note that this is used only in test found under UNIVARIATE.TESTS.MULTI.STRATA.FULL.PERM.
EARLY_STOP = F # note that I didn't implement a wald sequential per se, just a stopping rule,
# A wald sequential with multiplicity has a problem in maintaining good power.
# Specifically, you need to adjust alpha0 and beta0 in the wald seq test for multiplicity.
# Since weve dropped permutations tests as our recommended approach, I am not investing more effort
# in this. If you want to play with Wald Sequential for this framework, feel free to contact me.

#Will stop at this number of permutations:
STOPPING_RULE_PERM = 1000
# If P-val will be above:
STOPPING_RULE_PV = 0.2
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exact stats for permutations:
# If pack states is called, will calculate the percent of permutations with P-values under the values found in DS.FDR.Grid
DS.FDR.Grid.Constant = 0.00001
DS.FDR.Grid = c(DS.FDR.Grid.Constant*c(1:50,seq(60,200,10),seq(250,0.1/DS.FDR.Grid.Constant,by=50)))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Definitions for the different types of tests.

# All tests available
UNIVARIATE.TESTS.MULTI.STRATA.TYPES = c('KW',                     # A KW test with block permutations. May be used with mean value adjustment per block
                                        'Spearman',               # A Spearman test with block permutations. May be used with mean value adjustment per block
                                        'Spearman-Fisher',        # Asymptotic Spearman tests in each strata, aggregating over the P-value of different blocks. Aggregation is done by Fisher combination over right sided tests and left sided tests seperately, with the PV taken to be the bonferroni PV of both Fisher combinations.
                                        'KW-Fisher',              # Asymptotic KW tests in each strata, aggregating over two sided P-values using Fisher combination scores
                                        'Hoeffding-Fisher',       # Asymptotic Hoeffding tests in each strata, aggregating over two sided P-values using Fisher combination scores
                                        'KW-Strata-Asymp',        # Asymptotic KC tests in each strata, aggregating over tests statistics by summation, as implemented in the 'coin' package
                                        'Spearman-Strata-Asymp',  # Asymptotic Spearman tests in each strata,  aggregating over tests statistics by summation, as implemented in the 'coin' package
                                        'Wilcoxon-Fisher',        # Asymptotic Wilcoxon tests in each strata,  aggregating over the P-value of different blocks. Aggregation is done by Fisher combination over right sided tests and left sided tests seperately, with the PV taken to be the bonferroni PV of both Fisher combinations.
                                        "Wilcoxon-Strata-Asymp",  # Asymptotic Wilcoxon tests in each strata,  aggregating over tests statistics by summation, as implemented in the 'coin' package
                                        "Wilcoxon-Paired-Asymp",  # Asymptotic Wilcoxon tests for paired data, as implemented in the coin package.
                                        'C_Wilcoxon',             # A C variant of Wilcoxon, as implemented in dacomp.test
                                        'C_Spearman',             # A C variant of Spearman, as implemented in dacomp.test
                                        'Wilcoxon-Paired'         # A C level implementation of the Wilcoxon paired test, with P-value computed by permutations.
)

#These tests use block permutations.
UNIVARIATE.TESTS.MULTI.STRATA.FULL.PERM   = c('KW','Spearman','C_Wilcoxon','C_Spearman','Wilcoxon-Paired')
UNIVARIATE.TESTS.MULTI.STRATA.FULL.PERM.C = c('C_Wilcoxon','C_Spearman','Wilcoxon-Paired')

#These tests are for K-Samples, 2-Samples or paired.
UNIVARIATE.TESTS.MULTI.STRATA.TYPES.KS = c('KW','C_Wilcoxon','KW-Strata-Asymp','KW-Fisher','Wilcoxon-Fisher',"Wilcoxon-Strata-Asymp","Wilcoxon-Paired-Asymp",'Wilcoxon-Paired')

#These tests are for numeric or ordered traits
UNIVARIATE.TESTS.MULTI.STRATA.TYPES.ORDER = c('Spearman','C_Spearman','Spearman-Fisher','Spearman-Strata-Asymp','Hoeffding-Fisher')

# These tests will filter small blocks, as defined in minimal block size
UNIVARIATE.TESTS.MULTI.STRATA.TYPES.ASYMP.FILTERING.REQUIRED = c('Spearman-Fisher','KW-Fisher','Hoeffding-Fisher','KW-Strata-Asymp','Spearman-Strata-Asymp','Wilcoxon-Fisher',"Wilcoxon-Strata-Asymp","Wilcoxon-Paired-Asymp")

# These tests support categorical traits with only two values.
UNIVARIATE.TESTS.MULTI.STRATA.TYPES.KS.ONLY.TWO.LEVELS = c('Wilcoxon-Fisher',"Wilcoxon-Strata-Asymp","Wilcoxon-Paired-Asymp",'C_Wilcoxon','Wilcoxon-Paired')

# These tests aggregate over strata by P-values combination methods
UNIVARIATE.TESTS.MULTI.STRATA.TYPES.AGGREGATE.OVER.BLOCKS = c("Spearman-Fisher","KW-Fisher","Hoeffding-Fisher","Wilcoxon-Fisher")

# Tests for paired designs.
UNIVARIATE.TESTS.MULTI.STRATA.TYPES.PAIRED.TEST = c("Wilcoxon-Paired-Asymp",'Wilcoxon-Paired')

# Auxilliary function for computing the percent of permutation P-values which are under set thresholds. Used in DSFDR multiplicity corrections
Univariate.Test.Multi.Strata.DSFDR.aux.function=function(stats_values){
  ret = rep(NA,length(DS.FDR.Grid))
  temp_pv_ranked = rank(c(stats_values,Inf),
                        na.last = NA,
                        ties.method = 'min')
  temp_pv_ranked =  (length(temp_pv_ranked)- temp_pv_ranked)/length(temp_pv_ranked)
  for(i in 1:length(DS.FDR.Grid)){
    ret[i] = mean(temp_pv_ranked <= DS.FDR.Grid[i])
  }
  return(ret)
}

# Auxilliary function for packing the results of a univariate tests.
Univariate.Test.Multi.Strata.Pack.Results=function(stats_values,Description,
                                                   pack_stats = F,
                                                   is_PV_given = F,
                                                   given_PV = NA){
  ret = list()
  if(is_PV_given){
    temp_pv = given_PV
  }else{
    temp_pv = mean(stats_values >= stats_values[1],na.rm = T)
    ret$Number_of_perms_under_value = Univariate.Test.Multi.Strata.DSFDR.aux.function(stats_values)
  }
  
  ret$P.value = temp_pv
  ret$Description = Description
  
  if(pack_stats){
    ret$stats_values = stats_values
  }
  return(ret)
}

#' Function for generating permutations of trait value Y by the blocks given in Z.
#'
#' @param Z stratification variable, values of \code{Y} will only be permuted in blocks of this variable
#' @param Y variable to generate permutations of
#' @param nr.perm number of permutations to generate
#'
#' @return A matrix with \code{nr.perm} columns, giving block permutations for the values of \code{Y}, by the strata defined by the values of \code{Z}
#' @export
#'
#' @examples
Univariate.Test.Generate.Permutations.Object=function(Z,Y,nr.perm = 200){
  
  # We identify the different blocks
  Z.values = unique(Z)
  Z.ind.list = list()
  ind = 1:length(Z)
  
  for(Zi in 1:length(Z.values)){
    current_block = which(Z == Z.values[Zi])
    if(length(current_block)>1){
      Z.ind.list[[length(Z.ind.list) + 1]] = current_block
    }
  }
  
  # and generate a permutation matrix by over indices (perm_matrix), the matrix (Y_perm_matrix) will contain the actual permuted values
  perm_matrix = matrix(NA,nrow = length(Z),ncol = nr.perm+1)
  Y_perm_matrix = matrix(NA,nrow = length(Z),ncol = nr.perm+1)
  for(j in 1:(nr.perm+1)){
    perm_matrix[,j] = 1:length(Z)
    for(Zi in 1:length(Z.ind.list)){
      perm_matrix[Z.ind.list[[Zi]],j] = sample(ind[Z.ind.list[[Zi]]])
    }
    Y_perm_matrix[,j] = Y[perm_matrix[,j]]
  }
  
  return(Y_perm_matrix)
}



#' A mega-wrapper function for computing non parametric tests over strata.
#'
#' Run the test defined by \code{Method} to test for association between ranks of \code{X} and \code{Y}, while adjusting for the strata given by \code{Z}. The function supports several testing methods, for several possible designs, this includes an ordered Y, categorical Y, or paired design (Y has two values, and Z defined pairs of observations).
#'
#' @details This function provides a mega-wrapper for many different rank tests. \code{X} is assumed to be normalizaed micorbiome counts (after DACOMP or DACOMP-ratio normalization). \code{Y} is either numeric for Hoeffding or Spearman type tests, or a character vector for K-sample, 2-Sample and paired design tests.
#' The method begins by ranking X, with ties in observation replaced by average ranks. The flag \code{do.block.mean.normalization} sets whther the mean value in each block of \code{X} (as defined by \code{Z}) is subtracted from observations after ranking. This is relevant for tests which compute the test statistic over all strata combined ('Spearman' or 'KW' tests, found in \code{UNIVARIATE.TESTS.MULTI.STRATA.FULL.PERM}) but also for the description string returned (see value for details).
#' The list of tests available is defined by \code{UNIVARIATE.TESTS.MULTI.STRATA.TYPES} and also listed below.
#' Tests with the word 'Fisher' in name (such as 'Spearman-Fisher','Wilcoxon-Fisher','Hoeffding-Fisher','KW-Fisher') combine over tests in block using the Fisher combination score. For 'Wilcoxon-Fisher' and 'Spearman-Fisher', two Fisher combination scores are computed, for right and left sided tests across block, and the final P-value is the bonferroni combination of the two P-values.
#' Tests with the word 'Strata' in name (such as 'KW-Strata-Asymp', 'Spearman-Strata-Asymp','Wilcoxon-Strata-Asymp')aggregate over test statistics by summation over the different strata, as implemented in the \code{\link{coin}} package
#'
#' The different tests available (values for 'Method') are:
#' \itemize{
#' \item{'KW'}{A KW test with block permutations. May be used with mean value adjustment per block}
#' \item{'Spearman'}{A Spearman test with block permutations. May be used with mean value adjustment per block}
#' \item{'Spearman-Fisher'}{Asymptotic Spearman tests in each strata, aggregating over the P-value of different blocks. Aggregation is done by Fisher combination over right sided tests and left sided tests seperately, with the PV taken to be the bonferroni PV of both Fisher combinations.}
#' \item{'KW-Fisher'}{Asymptotic KW tests in each strata, aggregating over two sided P-values using Fisher combination scores}
#' \item{'Hoeffding-Fisher'}{Asymptotic Hoeffding tests in each strata, aggregating over two sided P-values using Fisher combination scores}
#' \item{'KW-Strata-Asymp'}{Asymptotic KC tests in each strata, aggregating over tests statistics by summation, as implemented in the 'coin' package}
#' \item{'Spearman-Strata-Asymp'}{ Asymptotic Spearman tests in each strata,  aggregating over tests statistics by summation, as implemented in the 'coin' package}
#' \item{'Wilcoxon-Fisher'}{Asymptotic Wilcoxon tests in each strata,  aggregating over the P-value of different blocks. Aggregation is done by Fisher combination over right sided tests and left sided tests seperately, with the PV taken to be the bonferroni PV of both Fisher combinations.}
#' \item{'Wilcoxon-Strata-Asymp'}{ Asymptotic Wilcoxon tests in each strata,  aggregating over tests statistics by summation, as implemented in the 'coin' package}
#' \item{"Wilcoxon-Paired-Asymp"}{Asymptotic Wilcoxon tests for paired data, as implemented in the coin package.}
#' \item{"C_Wilcoxon"}{ A C variant of Wilcoxon, as implemented in dacomp.test}
#' \item{"C_Spearman"}{ A C variant of Spearman, as implemented in dacomp.test}
#'
#' Description, P-value, permutation stats and permutation statistics.
#' }
#' @param X A numeric variable. For microbiome testing using DACOMP, it is assumed X is given after normalization for compositionality (DACOMP or DACOMP-ratio).
#' @param Y If \code{Method} is in \code{UNIVARIATE.TESTS.MULTI.STRATA.TYPES.KS} (KW, Wilcoxon or paired tests) a character vector defining groups. If \code{Method} is in \code{UNIVARIATE.TESTS.MULTI.STRATA.TYPES.ORDER } (Spearman or Hoeffding tests) a vector of numeric values.
#' @param Z A character vector, defining the strata for testing.
#' @param Method The type of test to perform, must be one of the methods available under \code{UNIVARIATE.TESTS.MULTI.STRATA.TYPES}. See details for explanation on the different methods.
#' @param nr.perm Number of permutations to perform. This is relevant only for tests which perform actual permutations by blocks, and
#' @param do.block.mean.normalization Whether the mean value of ranks of \code{X} in each stratum should be subtracted before testing. This is relevant for two elements: 1) For the permutations tests in \code{UNIVARIATE.TESTS.MULTI.STRATA.FULL.PERM}, this changes the computation method (are the values of ranks of \code{X}, with or without block adjustment); and 2) The description field will change by the block adjustment (even for tests aggregating over blocks).
#' @param given_permutation A table defining the permutations of \code{Y}, as generated by \code{Univariate.Test.Generate.Permutations.Object}. If generated manually, the table must contain a number of columns equal to \code{nr.perm}, and rows by the length of \code{Y}. Each column must contain the values of \code{Y}, with values permuted by the strata defined by \code{Z}
#' @param pack_stats For tests defined under \code{UNIVARIATE.TESTS.MULTI.STRATA.FULL.PERM}, turning this flag to \code{T} will return a vector of test statistics, the first one for the original data and the rest computed for permuted values of \code{Y} (by blocks of \code{Z}). Relevant only for tests in \code{UNIVARIATE.TESTS.MULTI.STRATA.FULL.PERM}.
#' @param Minimum_Block_Size Tests that aggregate over different strata either by summing test statistics in strata (have 'strata' in name), or by combining P-values (have 'Fisher' in name) may lose power if aggregating over small strata. Tests in \code{UNIVARIATE.TESTS.MULTI.STRATA.TYPES.ASYMP.FILTERING.REQUIRED} will disregard strata having fewer parameters than defined by this function.
#' @param print.permutations.progress Whether the progress of permutations should be printed to screen.
#'
#' @return
#' \itemize{
#' \item{P.value}{P-value for the test}
#' \item{Description }{For correlation tests (Spearman), will give the spearman correlation coeffcient between X and Y. For K-sample, 2-Sample and paired tests, will provide a string describing the ordering of mean rank of X, across levels of Y. Note that the description changes by the value of \code{do.block.mean.normalization}, even if the test performed if of the 'Strata' or 'Fisher' types. }
#' \item{Number_of_perms_under_value}{ percent of permutation P-values with values under \code{DS.FDR.Grid}}
#' \item{stats_values}{If \code{pack_stats} is set to true, this will hold a vector of test statistics, the first one for the original data and the rest computed for permuted values of \code{Y} (by blocks of \code{Z}). Relevant only for tests in \code{UNIVARIATE.TESTS.MULTI.STRATA.FULL.PERM}.}
#' }
#' @export
#'
#' @examples
#'
#'
#'
Univariate.Test.Multi.Strata = function(X,Y,Z,Method = 'KW',nr.perm = 200,
                                        do.block.mean.normalization = T,
                                        given_permutation = NULL,
                                        pack_stats = F,
                                        Minimum_Block_Size = 2,
                                        print.permutations.progress = F){
  Description = 'NONE'
  
  if(max(X)-min(X)<= 10^(-8)){
    return(Univariate.Test.Multi.Strata.Pack.Results(c(1,1),Description,pack_stats = F,is_PV_given = T,given_PV = NA))
  }
  
  Z.values = unique(Z)
  Z.ind.list = list()
  X_rank = X
  if(!(Method %in% UNIVARIATE.TESTS.MULTI.STRATA.TYPES)){
    stop('Univariate method not known')
  }
  
  is_full_perm = (Method %in% UNIVARIATE.TESTS.MULTI.STRATA.FULL.PERM)
  is_full_perm_c = (Method %in% UNIVARIATE.TESTS.MULTI.STRATA.FULL.PERM.C)
  
  for(Zi in 1:length(Z.values)){
    current_block = which(Z == Z.values[Zi])
    if(length(current_block)>1){
      Z.ind.list[[length(Z.ind.list) + 1]] = current_block
      if(do.block.mean.normalization)
        X_rank[current_block] = X_rank[current_block] - mean(X_rank[current_block]) #might be worthwhile to try with median, but need to check for zeros...
    }else{
      X_rank[current_block] = 0
    }
    
    if(sum(is.na(X_rank[current_block]))>0){
      print(paste0('Error:  block id :', Zi,' with NAs:',sum(is.na(X_rank[current_block]))))
    }
  }
  
  X_rank = rank(X_rank,ties.method = 'average')
  
  Y_perm = Y
  
  if(Method %in% UNIVARIATE.TESTS.MULTI.STRATA.TYPES.KS){
    group.values = unique(Y)
    group.sizes = rep(NA,length(group.values))
    group.sums = rep(NA,length(group.values))
    for(g in 1:length(group.values)){
      group.sizes[g] = sum(Y==group.values[g])
    }
  }else if (Method %in% UNIVARIATE.TESTS.MULTI.STRATA.TYPES.ORDER){
    X_rank = X_rank - mean(X_rank)
    Y_perm = rank(Y_perm,ties.method = 'average')
    Y_perm = Y_perm - mean(Y_perm)
  }
  
  #two families of tests:
  if(is_full_perm & !is_full_perm_c){
    # need to call function
    return(Univariate.Test.Multi.Strata.Wrapper.for.Perm.Based.Test(nr.perm,
                                                                    Method,
                                                                    group.values,
                                                                    group.sums,
                                                                    group.sizes,
                                                                    X_rank,
                                                                    Y_perm,
                                                                    given_permutation,
                                                                    EARLY_STOP,
                                                                    print.permutations.progress,
                                                                    Z.ind.list,
                                                                    pack_stats))
  }else if(is_full_perm_c){
    # calling the C level tests from DACOMP
    if(Method == 'Wilcoxon-Paired' & !is.null(given_permutation) ){
      warning('Pregenerated permutation matrix given for Wilcoxon-Paired, but is not used.')
    }
    if(!(Method == 'Wilcoxon-Paired')){ #currently we are not supporting pregenerated permutations for Wilcoxon-Paired. It is too costly to generate the permutation matrix and iterate over it.
      if(is.null(given_permutation) | Method == 'C_Spearman'){ #currently we are not supporting pregenerated normalization for spearman, we perform normalization on Y_perm before this point
        given_permutation = Univariate.Test.Generate.Permutations.Object(Z,Y_perm,nr.perm = nr.perm)
      }
      perm_object = cbind(Y_perm,given_permutation)
    }
    if(Method %in% c('Wilcoxon-Paired')){
      ind_g1 = (Y_perm == Y_perm[1])
      ind_g2 = (Y_perm != Y_perm[1])
      X1 = X_rank[ind_g1]
      X2 = X_rank[ind_g2]
      X1 = X1[order(Z[ind_g1])]
      X2 = X2[order(Z[ind_g1])]
      if(all(X1-X2==0)){
        return(Univariate.Test.Multi.Strata.Pack.Results(stats_values = c(1,1),Description = "All pair differences are zero, cannot test",pack_stats = pack_stats,is_PV_given = T,given_PV = NA))
      }
      res = rep(NA,nr.perm)
      
      
      differences = X1-X2
      second_is_bigger = 1*(differences > 0)
      differences[differences==0] = NA
      rank_differences = rank(abs(differences), ties.method = 'average',na.last = "keep")
      
      res[1] = dacomp:::rcpp_Compute_Wilcoxon_Signed_Rank_Stat(rank_differences,second_is_bigger) #(dacomp:::SignedRankWilcoxon.statistic(X1,X2))^2
      
      for(b in 2:nr.perm){
        second_is_bigger = rbinom(length(X1),size = 1,prob = 0.5)  
        res[b] = dacomp:::rcpp_Compute_Wilcoxon_Signed_Rank_Stat(rank_differences,second_is_bigger) #(dacomp:::SignedRankWilcoxon.statistic(p_X1,p_X2))^2
        
      }
      
      stat_mean = sum(rank_differences,na.rm = T)/2
      stat_var = sum(rank_differences^2,na.rm = T)/4
      res = (res - stat_mean)/sqrt(stat_var)
      res = res^2
      Description = description_for_KS(data.frame(X_rank = X_rank,Y_perm = Y_perm))
    }
    if(Method %in% c('C_Wilcoxon')){
      if(is.factor(Y_perm))
        perm_object = 1*(perm_object == as.numeric(Y_perm[1]))
      else
        perm_object = 1*(perm_object == (Y_perm[1]))
      
      test_as_param = dacomp::DACOMP.TEST.NAME.WILCOXON
      res = dacomp:::Compute.resample.test(X_matrix = t(t(X_rank)),Y_matrix = perm_object,statistic = test_as_param)
      
      Description = description_for_KS(data.frame(X_rank = X_rank,Y_perm = Y_perm))
      
    }
    if(Method == 'C_Spearman'){
      res = dacomp:::Compute.resample.test(X_matrix = t(t(X)),Y_matrix = perm_object,statistic = dacomp::DACOMP.TEST.NAME.SPEARMAN)
      Description = cor(X_rank,Y_perm,method = 'spearman')
    }
    return(Univariate.Test.Multi.Strata.Pack.Results(stats_values = res,Description = Description,pack_stats = pack_stats,is_PV_given = F))
    
  }else{
    # we are left with asymptotic tests
    return(Univariate.Test.Multi.Strata.Wrapper.for.Asymp.Based.Test(Method,
                                                                     X_rank,
                                                                     Y_perm,
                                                                     Z.ind.list,
                                                                     Z,
                                                                     Minimum_Block_Size))
    
  }
  
}

check.dt_filtering_for_Coin_Wilcoxn_Paired = function(dt_filtering){
  
  if(any(!(table(dt_filtering$Z) %in% c(0,2)))){
    stop('Strata defined by Z are not only pairs')
  }
  
  ind_g1 = which(dt_filtering$Y_perm == unique(dt_filtering$Y_perm)[1])
  ind_g2 = which(dt_filtering$Y_perm == unique(dt_filtering$Y_perm)[2])
  
  xval_by_sorted_Z_g1 = (dt_filtering$X_rank[ind_g1])[order(dt_filtering$Z[ind_g1])]
  xval_by_sorted_Z_g2 = (dt_filtering$X_rank[ind_g2])[order(dt_filtering$Z[ind_g2])]
  
  if(!all( xval_by_sorted_Z_g1 - xval_by_sorted_Z_g2 == 0)){
    return(T)
  }else{
    return(F)
  }
}


# An internal function for performing tests based on permutations
Univariate.Test.Multi.Strata.Wrapper.for.Perm.Based.Test=function(nr.perm,
                                                                  Method,
                                                                  group.values,
                                                                  group.sums,
                                                                  group.sizes,
                                                                  X_rank,
                                                                  Y_perm,
                                                                  given_permutation,
                                                                  EARLY_STOP,
                                                                  print.permutations.progress,
                                                                  Z.ind.list,
                                                                  pack_stats){
  
  stats_values = rep(NA,nr.perm + 1) # will store stat for original data and permutation stats
  #handle statistics with full permutations
  
  #over permutations
  for( b in 1:(nr.perm+1)){
    
    if(print.permutations.progress)
      print(paste0('Permutation : ',b))
    
    #compute statistics. Note the for b==1, the Y_perm values haven't been permuted yet, so we are computing it for the original data
    
    if(Method == 'KW'){
      #Compute the KW stat
      for(g in 1:length(group.values)){
        group.sums[g] = sum(X_rank[Y_perm == group.values[g]])
      }
      overall_mean = sum(group.sizes*group.sums)/length(Y_perm)
      
      stats_values[b]  = sum(group.sizes *(group.sums - overall_mean)^2)/max(sum((X_rank - overall_mean)^2),1)
      #For the null permutation, we create a description ordering the mean ranks of groups from lowest to highest
      if(b==1){
        Description = paste(paste0(group.values,' (',round(group.sums/group.sizes),')')[order(group.sums/group.sizes)],collapse = '<=')
      }
    }else if (Method == 'Spearman'){
      #Compute the spearman test statistic with description
      if(b==1){
        Description = as.character(cor(X_rank,Y_perm,method = 'spearman'))
      }
      stats_values[b] = (sum(X_rank*Y_perm))^2
    }
    
    #permute. If permutations are given, take them from object
    if(!is.null(given_permutation)){
      Y_perm = given_permutation[,b]
    }else{
      # If permutations not given, we permute in each block.
      for(Zi in 1:length(Z.ind.list)){
        Y_perm[Z.ind.list[[Zi]]] = sample(Y_perm[Z.ind.list[[Zi]]])
      }
    }
    
    #If early stopping is activated, check if we need to stop.
    if(EARLY_STOP){
      temp_pv = mean(stats_values >= stats_values[1],na.rm = T)
      if((b == STOPPING_RULE_PERM & temp_pv >= STOPPING_RULE_PV )){
        return(Univariate.Test.Multi.Strata.Pack.Results(stats_values,Description,pack_stats = pack_stats))
      }
    }
  }
  
  #done all permutations on a full perm test, return results
  return(Univariate.Test.Multi.Strata.Pack.Results(stats_values,Description,pack_stats = pack_stats))
}


# A wrapper for Asymptotic tests.
Univariate.Test.Multi.Strata.Wrapper.for.Asymp.Based.Test = function(
  Method,
  X_rank,
  Y_perm,
  Z.ind.list,
  Z,
  Minimum_Block_Size){
  #check the user hasn't required Minimum block size of 2 in paired tests. If the user requested that, override - obviously it is a mistake
  if(Method %in% UNIVARIATE.TESTS.MULTI.STRATA.TYPES.PAIRED.TEST & Minimum_Block_Size >2){
    warning('Minimum block size requested was higher than 2 in a paired test, setting it to 2')
    Minimum_Block_Size = 2
  }
  
  #we create a data table, for processing strata and iterating over them.
  dt_filtering = data.frame(X_rank=X_rank,Y_perm=Y_perm,Z=Z)
  
  #remove small blocks , but only if needed
  if(Method %in% UNIVARIATE.TESTS.MULTI.STRATA.TYPES.ASYMP.FILTERING.REQUIRED){
    dt_filtering$Z =(as.factor(dt_filtering$Z))
    Z_to_keep = names(which(table(dt_filtering$Z) >= Minimum_Block_Size ))
    dt_filtering = dt_filtering[dt_filtering$Z %in% Z_to_keep,,drop=F]
    if(nrow(dt_filtering) == 0){
      stop('No block larger than minimal block size required')
    }
    # for KS tests that require filtering, check each strata has at least two trait values in Y (otherwise - nothing to sum over in test in stratum)
    if(Method %in% UNIVARIATE.TESTS.MULTI.STRATA.TYPES.KS){
      Z_to_keep = names(which(apply(table(dt_filtering$Y_perm,dt_filtering$Z) >0,2,sum) >= 2))
      dt_filtering = dt_filtering[dt_filtering$Z %in% Z_to_keep,,drop = F]
      if(nrow(dt_filtering) == 0){
        stop('No blocks with Y values in at least two levels')
      }
    }
  }
  
  #check that for 2-Samples tests Y has exactly two values.
  if(Method %in% UNIVARIATE.TESTS.MULTI.STRATA.TYPES.KS.ONLY.TWO.LEVELS){
    if(length(unique(dt_filtering$Y_perm))>2)
      stop('Wilcoxon-type test (or other two sample test) with more than two levels in trait variable')
  }
  
  # Generate description field. For order tests, this is the spearman correlation coffeficient.
  # For KS tests, we generate a string ordering group labels of Y by their average rank over X values.
  # Note that the description is affected by block mean adjustment.
  # Even if "Strata" or "Fisher" type tests are computed in blocks, the description is computed while disregarding blocks (hence, it makes sense to always have mean block adjustment turned on).
  if(Method %in% UNIVARIATE.TESTS.MULTI.STRATA.TYPES.ORDER){
    Description = as.character(cor(dt_filtering$X_rank,dt_filtering$Y_perm,method = 'spearman'))
  }else{
    Description = description_for_KS(dt_filtering)
  }
  
  # we have to run over blocks, run tests, and collect results
  # The tests handled here are #"Spearman-Fisher","KW-Fisher","Hoeffding-Fisher","Wilcoxon-Fisher"
  if(Method %in% UNIVARIATE.TESTS.MULTI.STRATA.TYPES.AGGREGATE.OVER.BLOCKS){
    
    unique_Z = unique(dt_filtering$Z)
    p.val.vec = rep(NA,length(unique_Z))
    
    #iterate over strata
    for(Zi in 1:length(unique_Z)){
      X1 = dt_filtering$X_rank[dt_filtering$Z == unique_Z[Zi]]
      Y1 = dt_filtering$Y_perm[dt_filtering$Z == unique_Z[Zi]]
      dt_temp = data.frame(X1 = X1,Y1=Y1)
      # run test in strata:
      if(Method == 'KW-Fisher'){
        
        test_res = coin::kruskal_test(X1~Y1,data = dt_temp,
                                      distribution = "asymptotic")
        p.val.vec[Zi] = coin::pvalue(test_res)
        
      }else if(Method == 'Wilcoxon-Fisher'){
        test_res = coin::wilcox_test(X1~Y1,data = dt_temp,
                                     distribution = "asymptotic",
                                     alternative = "greater")
        
        p.val.vec[Zi] = coin::pvalue(test_res)
        
      }else if(Method == 'Spearman-Fisher'){
        
        test_res = coin::spearman_test(X1~Y1,distribution = "asymptotic",
                                       alternative = "greater")
        p.val.vec[Zi] = coin::pvalue(test_res)
      }else if(Method == 'Hoeffding-Fisher'){
        
        Hoeffding_stat = Hmisc::hoeffd(X1,Y1)
        p.val.vec[Zi] = Hoeffding_stat$P[1,2]
        
      }
    } # done iterating over strata
    
    #combine strata
    p.val.vec = p.val.vec[(!is.na(p.val.vec))] #remove strata which were NA (for example, if all X values were identical)
    All_Removed = F
    if(length(p.val.vec)==0)
      All_Removed = T
    if(!All_Removed){
      if(Method %in% c('Hoeffding-Fisher','KW-Fisher')){
        #P-values for two-sided tests are aggregated by Fisher combination scores
        test_statistic = -2 * sum(log(p.val.vec))
        pv.Fisher = 1-pchisq(test_statistic,df = 2*length(p.val.vec),lower.tail = T)
      }else if(Method %in% c('Spearman-Fisher','Wilcoxon-Fisher')){
        #For 'Spearman-Fisher','Wilcoxon-Fisher', we can test directional hypotheses so it makes sense to Fisher-combine right sided tests and left sided tests seperatly, and then Bonferroni combine both Fisher P-values
        test_statistic_G = -2 * sum(log(p.val.vec))
        test_statistic_L = -2 * sum(log(1-p.val.vec))
        pv.Fisher.G = 1-pchisq(test_statistic_G,df = 2*length(p.val.vec),lower.tail = T)
        pv.Fisher.L = 1-pchisq(test_statistic_L,df = 2*length(p.val.vec),lower.tail = T)
        pv.Fisher = min(2*min(pv.Fisher.G,pv.Fisher.L),1)
      }
      return(Univariate.Test.Multi.Strata.Pack.Results(c(-1,-1),Description,pack_stats = F,is_PV_given = T,given_PV = pv.Fisher))
    }else{
      # we got here if all strata were discarded, and than we can simply put NA as the P-values and avoid testing this hypothesis
      return(Univariate.Test.Multi.Strata.Pack.Results(c(-1,-1),'All strata filtered',pack_stats = F,is_PV_given = T,given_PV = NA))
    }
    
    # Done over methods that require aggregation over strata
    
  }else{
    dt_filtering$Z = as.factor(as.character(dt_filtering$Z))
    # P-value can be evaluated by a single function call.
    # The tests handled here are "KW-Strata-Asymp" "Spearman-Strata-Asymp" "KW-Asymp" "Spearman-Asymp"
    if(Method == "Spearman-Strata-Asymp"){
      test_res = coin::spearman_test(X_rank~Y_perm|Z,data = dt_filtering, distribution = "asymptotic")
      pv = coin::pvalue(test_res)
    }else if(Method == "KW-Strata-Asymp"){
      test_res = coin::kruskal_test(X_rank~Y_perm|Z,data = dt_filtering, distribution = "asymptotic")
      pv = coin::pvalue(test_res)
    }else if(Method == "Wilcoxon-Strata-Asymp"){
      test_res = coin::wilcox_test(X_rank~Y_perm|Z,data = dt_filtering, distribution = "asymptotic")
      pv = coin::pvalue(test_res)
    }else if(Method == "Wilcoxon-Paired-Asymp"){
      check_data_for_coin_signed_wilcoxon = check.dt_filtering_for_Coin_Wilcoxn_Paired(dt_filtering)
      
      if(check_data_for_coin_signed_wilcoxon){
        test_res = coin::wilcoxsign_test(X_rank~Y_perm|Z,data = dt_filtering, distribution = "asymptotic");
        pv = coin::pvalue(test_res)                    
      }else{
        #cannot test wilcoxon - all pair differences are zero
        pv = NA
      }
    }
    return(Univariate.Test.Multi.Strata.Pack.Results(c(-1,-1),Description,pack_stats = F,is_PV_given = T,given_PV = pv))
  }
  
}

#Aux function for providing a description, for KW-Strata-Asymp, KW-Fisher and Wilcox-Fisher
description_for_KS = function(dt_for_kSamples){
  aggregated_means = aggregate(x = dt_for_kSamples$X_rank,by=list(Y_perm = dt_for_kSamples$Y_perm),mean)
  aggregated_means = aggregated_means[order(aggregated_means$x),]
  return(paste0(aggregated_means$Y_perm,collapse = '<='))
}
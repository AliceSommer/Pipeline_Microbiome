
#' Function for DACOMP tests with stratification
#'
#' Run non-parametric tests of association, with strata.
#' This function may be used for computing P-values for marginal tests (for reference selection, \code{taxa_to_normalize_by = NULL}), or for running tests with DACOMP-ratio or DACOMP normalization (\code{taxa_to_normalize_by} set to the references of taxa indices)
#'
#' @details The function performs univariate non-parameteric tests for differentiall abundance. If no reference taxa are selected (\code{taxa_to_normalize_by = NULL}), counts are normalized by the total number of counts in sample. This form of analysis us used for selecting non-rare taxa that are stable with respect to change in \code{Y} across different taxa, and form a reference set using them. If reference taxa are selected, the function performs DACOMP or DACOMP-ratio normalization, as in the function \code{\link{dacomp.test}}.
#' The function calls \code{\link{Univariate.Test.Multi.Strata}} over normalized taxa counts for testing. The types of tests available are defined in \code{\link{Univariate.Test.Multi.Strata}}.
#' The parameter \code{Method} sets the type of test. The parameters \code{do.block.mean.normalization, nr.perm, Minimum_Block_Size} are used as parameters when calling \code{\link{Univariate.Test.Multi.Strata}}.
#' See \code{\link{Univariate.Test.Multi.Strata}} for a full description on the types of tests available.
#'
#' @param X A matrix of 16S counts, un normalized
#' @param Y A vector of Trait values
#' @param Z A chacter vector, specifying the stratification.
#' @param Method test type, as defined in \code{\link{Univariate.Test.Multi.Strata}}
#' @param do.block.mean.normalization Should block stratification be performed? see \code{\link{Univariate.Test.Multi.Strata}} for details.
#' @param nr.perm Number of permutations for tests based on permuatations. See \code{\link{Univariate.Test.Multi.Strata}} for details on which tests are asymptotoc and which are permutation based.
#' @param Minimum_Block_Size Minimal size of strata. Some tests filter out small strata that may hinder statistical power. See \code{\link{Univariate.Test.Multi.Strata}} for a list of tests which use filtration.
#' @param taxa_to_normalize_by The set of taxa to normalize by. The default value of \code{NULL} performs TSS (Total Sum in Sample) normalization. Inputing a numeric vector sets a set of reference taxa (by indicies) to normalize by. You can also input the result of \code{\link{dacomp.select_references}} or \code{\link{dacomp.select_references.by.split}}
#' @param normalize_by_DACOMP_ratio Should the DACOMP-ratio normalization be used? (default:\code{T}). DACOMP-ratio normalization divids counts for a taxon by the total number of counts in the taxon and the reference. If set to \code{F}, DACOMP normalization with rarefaction will be used: the number of counts in a taxon used for analysis will be a hypergeometric sample from the total number of reads available under the taxon normalized and the reference taxa. The sampling depth is set by the maximum possible depth with no samples removed.
#' @param run.in.parallel Should computation be run in parallel.
#' @param select.taxa.for.DSFDR A logical (boolean) vector specifying for which taxa should the DSFDR adjusted P-values be computed. Taxa whose corresponding entries in this vector are set to F are excluded from testing, and will not have Pvalues are DSFDR adjusted Pvalues computed. When using this function to perform marginal differential abundance test,e.g., for selecting taxa to serve as a reference set, leave this vector at it's default value to make sure marginal tests are computed for all taxa.
#' @return A list of class \code{dacomp.strata.result.object}. The list will contain \code{P.values} - vector of P-values. If \code{taxa_to_normalize_by} is set to a value other than \code{NULL}, reference taxa would not be tested for differentiall abundance. If \code{Method} was set to a test computing P-values by permutations, you can use the returned object will also contain \code{DSFDR.AdjustedPvalues} a vector of DSFDR adjusted P-values. Note that taxa whose corresponding entries in the parameter argument \code{select.taxa.for.DSFDR} is set to \code{F} will not be tested for differential abundance, and will have P-values set to \code{NA}.
#' @export
#'
#' @examples
#' 
#' 
dacomp.test.with.strata =  function( X,
                                Y,
                                Z,
                                Method,
                                do.block.mean.normalization=T,
                                nr.perm=1000,
                                Minimum_Block_Size=4,
                                taxa_to_normalize_by = NULL,
                                normalize_by_DACOMP_ratio = T,
                                run.in.parallel = F,
                                select.taxa.for.DSFDR = rep(T,ncol(X)) ){
  
  
  Ref.Select.Start.Time = Sys.time()
  
  given_perm_object = NULL
  if(Method %in% c('KW','Spearman','C_Wilcoxon')){
    given_perm_object = Univariate.Test.Generate.Permutations.Object(Z,Y,nr.perm = nr.perm)
  }
  
  
  Skip_references = T 
  if(is.null(taxa_to_normalize_by)){
    Skip_references = F #if we do marginal testing, there is no "skipping the reference", we have to test all taxa
    taxa_to_normalize_by = 1:ncol(X)
  }else{
    #we check if we got a reference object. If so, we use only the selected refrences.
    if(class(taxa_to_normalize_by) == dacomp:::CLASS.LABEL.REFERENCE_SELECTION_OBJECT)
      taxa_to_normalize_by = taxa_to_normalize_by$selected_references
  }
  
  # counts total number of reads in reference
  total_counts = apply(X[,taxa_to_normalize_by],1,sum)

  single_pv_computer = function(j){
    
    if((Skip_references & (j %in% taxa_to_normalize_by)) | !select.taxa.for.DSFDR[j])
      return(NA)
    
    X_j = X[,j]
    
    #compute total number of counts in j and reference across taxa, note that j might also be in the reference
    if(j %in% taxa_to_normalize_by){
      total_counts_j = total_counts
    }else{
      total_counts_j = total_counts + X_j
    }
    
    #normalize, either by DACOMP or DACOMP ratio
    if(normalize_by_DACOMP_ratio){
      keep_Counts_marginal = X_j/total_counts_j  
    }else{
      lambda_j = min(total_counts_j)
      keep_Counts_marginal = rhyper(length(X_j),
                                    m = X_j,
                                    n = total_counts_j - X_j,
                                    k = lambda_j)
    }
    
    #test, with all relevant parameters
    temp = Univariate.Test.Multi.Strata(X = keep_Counts_marginal,
                                        Y = Y,
                                        Z = Z,
                                        Method = Method,
                                        do.block.mean.normalization = do.block.mean.normalization,
                                        nr.perm = nr.perm,
                                        Minimum_Block_Size = Minimum_Block_Size,given_permutation = given_perm_object)
    return(temp)
  }
  
  export_vars = c('EARLY_STOP',
  'STOPPING_RULE_PERM',
  'STOPPING_RULE_PV',
  'DS.FDR.Grid.Constant',
  'DS.FDR.Grid',
  'UNIVARIATE.TESTS.MULTI.STRATA.TYPES',
  'UNIVARIATE.TESTS.MULTI.STRATA.FULL.PERM',
  'UNIVARIATE.TESTS.MULTI.STRATA.FULL.PERM.C',
  'UNIVARIATE.TESTS.MULTI.STRATA.TYPES.KS',
  'UNIVARIATE.TESTS.MULTI.STRATA.TYPES.ORDER',
  'UNIVARIATE.TESTS.MULTI.STRATA.TYPES.ASYMP.FILTERING.REQUIRED',
  'UNIVARIATE.TESTS.MULTI.STRATA.TYPES.KS.ONLY.TWO.LEVELS',
  'UNIVARIATE.TESTS.MULTI.STRATA.TYPES.AGGREGATE.OVER.BLOCKS',
  'UNIVARIATE.TESTS.MULTI.STRATA.TYPES.PAIRED.TEST',
  'Univariate.Test.Multi.Strata.DSFDR.aux.function',
  'Univariate.Test.Multi.Strata.Pack.Results',
  'Univariate.Test.Generate.Permutations.Object',
  'Univariate.Test.Multi.Strata',
  'Univariate.Test.Multi.Strata.Wrapper.for.Perm.Based.Test',
  'Univariate.Test.Multi.Strata.Wrapper.for.Asymp.Based.Test',
  'check.dt_filtering_for_Coin_Wilcoxn_Paired',
  'description_for_KS')
  
  #run tests on all taxa by calling the above function
  if(run.in.parallel){
    cl <- makeCluster(parallel::detectCores() - 1)
    registerDoParallel(cl)  
    pvals_parallelized_result <- foreach(j=1:ncol(X), .options.RNG=1234,.export = export_vars) %dorng% {
      # source('Functions_for_univariate_tests.R')
      # source('dacomp_testing_and_reference_selection_by_split.R')
      return(single_pv_computer(j))
    }
    stopCluster(cl)  
  }else{
    pvals_parallelized_result = list()
    for(j in 1:ncol(X)){
      pvals_parallelized_result[[j]] = single_pv_computer(j)
    }
  }
  
  P_vals_vector = unlist(lapply(pvals_parallelized_result,
                                FUN = function(x){
                                if(is.null(names(x)))
                                  return(NA)
                                if(!('P.value' %in% names(x)))
                                  return(NA)
                                return(x$P.value)
                                }
                                ))
  
  ret = list()
  class(ret) = 'dacomp.strata.result.object'
  P.values = as.numeric(P_vals_vector) #this will drop annoying atributes coming from doRNG
  ret$P.values = P.values
  
  #This section computes DS-FDR adjusted P-values.
  DSFDR_filter = select.taxa.for.DSFDR #This is what the user requested, all == TRUE by default
  if(Skip_references){ #we dont test references
    DSFDR_filter[taxa_to_normalize_by] = F
  }
  for(j in 1:ncol(X)){ #we make sure to have the DSFDR statistics for what we want to test, and that a test was run properly, e.g. not everything was rarefied to zero
    if(!('Number_of_perms_under_value' %in% names(pvals_parallelized_result[[j]])) | is.na(P_vals_vector[j]))
      DSFDR_filter[j] = F
  }
  
  # If we have a non zero number of taxa to test:
  if(sum(DSFDR_filter)>0){
    Q.by.grid.value = rep(NA,length(DS.FDR.Grid)) #these will hold the extimates for E,V,Q (as by the DS-FDR paper) for the different possible critical P-values
    E.by.grid.value = rep(0,length(DS.FDR.Grid))
    V.by.grid.value = rep(1,length(DS.FDR.Grid))
    DS.FDR.AdjustedPvalues = rep(1,ncol(X)) # At the end of the algorihm, this will hold the DSFDR adjusted P-values by hypothesis
    
    #iterate over possible critical P-values
    for(k in 1:length(DS.FDR.Grid)){
      
      # compute the expected number of rejections, under the null, for a cutoff value.
      # This is computed by summing the probability for rejection, per hypothesis, as given by the 
      # proportion of permutation P-values that are under (or including) the critical P-value.
      for(j in 1:ncol(X)){
        if(DSFDR_filter[j])
          E.by.grid.value[k] = E.by.grid.value[k] + pvals_parallelized_result[[j]]$Number_of_perms_under_value[k]
      }
      
      # The number of rejections for a critical P-value is given the number of PVs smaller or equal than the threshold
      V.by.grid.value[k] = max(sum(P_vals_vector[DSFDR_filter]<=DS.FDR.Grid[k]),1)
      current_Q = E.by.grid.value[k]/V.by.grid.value[k] #this is the current estimated Qhat = E[V]/R for a given critical PV
      
      # The next code section updates updated the DSFDR adjusted P-values for each hypothesis,
      # by checking if the hypothesis was rejected at a lower Q than was ever recorded.
      temp_to_pmin_with_AdjustedPvalues = rep(1,ncol(X)) #we will use this to hold the current Q values. For non-rejections they will be 1, and for rejections we will put Qhat instead.
      
      #we check which of the hypotheses were rejected. rejected_for_current_PV_crit will have TRUE for the hypotheses that were rejected for the current PV crit 
      rejected_for_current_PV_crit = rep(F,ncol(X)) 
      ind_of_rejected_for_current_PV_crit = which(P_vals_vector[DSFDR_filter]<=DS.FDR.Grid[k])
      if(length(ind_of_rejected_for_current_PV_crit) >0){
        rejected_for_current_PV_crit[(which(DSFDR_filter == T))[ind_of_rejected_for_current_PV_crit]] = T
      }
      
      # for the rejected, we put Qhat instead of 1.
      if(sum(rejected_for_current_PV_crit)>0)
        temp_to_pmin_with_AdjustedPvalues[rejected_for_current_PV_crit] = current_Q
      
      # we take the element-wise (hypothesis-wise) minimum between the two vectors.
      DS.FDR.AdjustedPvalues = pmin(DS.FDR.AdjustedPvalues,temp_to_pmin_with_AdjustedPvalues)
    }
    
    # These are the Qhats by possible PVcrit.
    Q.by.grid.value = E.by.grid.value / V.by.grid.value
    names(Q.by.grid.value) = paste0('PV.crit = ', DS.FDR.Grid)
    
    
    ret$DSFDR.AdjustedPvalues = DS.FDR.AdjustedPvalues
    #ret$Q.by.grid.value = Q.by.grid.value # I am not exporting them now. No need to confuse the user with something that is not for the main pipeline. 

  }
  
  return(ret)
}


#' Delecting reference taxa, based on P-values from independent data
#'
#' The function selects reference taxa based on the P-values given by \code{P_vals_marginal}. Taxa will be selected into the reference set when by one, from the highest to lowest P-value,
#' until all taxa with P-values above and including \code{REFERENCE_SELECTION_PVAL_THRESHOLD} are selected, or the number of counts in the reference set of taxa is at least \code{MAXIMAL_TA} in all samples.
#' #' In any case, reference taxa will be selected until at least \code{MINIMAL_TA} counts are available under the reference taxa in each sample.
#' Taxa will be selected into the reference set, only if they are available in at least \code{MINIMAL_NR_SUBJECTS_FOR_REFERENCE_TEST} samples in the test data. For that, the counts matrix for the test data need to be supplied in \code{Counts_Mat_testing}.
#' 
#' @param P_vals_marginal A vector of P-values from marginal tests. Can also be the result of \code{dacomp.test.with.strata}.
#' @param Counts_Mat_testing The counts matrix for test data. Used for diregarding rare taxa, and for verifying the abundance of the test data.
#' @param REFERENCE_SELECTION_PVAL_THRESHOLD Taxa with marginal test P-values above and including this threshold will be allowed to enter the reference set of taxa
#' @param MINIMAL_NR_SUBJECTS_FOR_REFERENCE_TEST The minimal number of samples required to include a taxon, for the taxon to be considered a valid choice for the reference set of taxa.
#' @param MINIMAL_TA The selected reference set of taxa will include at least this number of counts in each sample.
#' @param MAXIMAL_TA The selected reference set of taxa will not include more than this number of counts in all samples.
#'
#' @return The function returns an object of type "dacomp.reference.selection.object", similar to the object returned from \code{dacomp.select_references}.
#' The main differene is that the median_SD field will hold the value \code{-log(REFERENCE_SELECTION_PVAL_THRESHOLD)}
#' 
#' @export
#'
#' @examples
dacomp.select_references.by.split = function(P_vals_marginal,
                                              Counts_Mat_testing,
                                              REFERENCE_SELECTION_PVAL_THRESHOLD = 0.5,
                                              MINIMAL_NR_SUBJECTS_FOR_REFERENCE_TEST = 5,
                                              MINIMAL_TA= 10,
                                              MAXIMAL_TA = 10000
                                                  ){
  
  if(class(P_vals_marginal) == 'dacomp.strata.result.object'){
    P_vals_marginal = P_vals_marginal$P.values
  }
  test_mat_prevalence = apply(1*(Counts_Mat_testing>0),2,sum)
  
  ind.nan = which(is.nan(P_vals_marginal) | is.na(P_vals_marginal) | test_mat_prevalence < MINIMAL_NR_SUBJECTS_FOR_REFERENCE_TEST)
  ind.not.nan = which(!is.nan(P_vals_marginal)& !is.na(P_vals_marginal) & test_mat_prevalence >= MINIMAL_NR_SUBJECTS_FOR_REFERENCE_TEST)
  
  if(length(ind.nan)>0)
    P_vals_marginal[ind.nan] = 1
  
  current_scores = -log(P_vals_marginal)
  current_scores[ind.nan] = max(current_scores) +1
  current_ref_select_for_variable = parallel_reference_select(X = Counts_Mat_testing,
                                                              median_SD_threshold = -log(REFERENCE_SELECTION_PVAL_THRESHOLD),
                                                              minimal_TA = MINIMAL_TA,
                                                              maximal_TA = MAXIMAL_TA,
                                                              Pseudo_Count_used = 1,
                                                              Precomputed_scores = current_scores,
                                                              select_from = ind.not.nan)
  
  Ref.Select.End.Time = Sys.time()
  
  return(current_ref_select_for_variable)
  
}

parallel_reference_select = function(X, median_SD_threshold, 
                                     minimal_TA = 10,
                                     maximal_TA = 200,
                                     Pseudo_Count_used = 1,
                                     verbose = F,
                                     select_from = NULL,
                                     Previous_Reference_Selection_Object = NULL,
                                     Nr.Cores = parallel::detectCores()-1,
                                     Precomputed_scores = NULL){
  #check inputs
  input_check_result = dacomp:::check.input.select_references(X, median_SD_threshold, minimal_TA,maximal_TA, Pseudo_Count_used, verbose, select_from, Previous_Reference_Selection_Object)
  if(!input_check_result)
    stop('Input check failed on dacomp.select_references')
  if(is.null(Previous_Reference_Selection_Object) & is.null(Precomputed_scores)){
    m = dim(X)[2]
    
    #compute the SD_{j,k}, see paper for definition
    ratio_matrix = matrix(NA, ncol = m, nrow = m)
    
    
    library(doParallel)
    cl <- makeCluster(Nr.Cores)
    registerDoParallel(cl)
    
    SDs_for_fixed_i = function(i){
      ret = rep(NA,m)
      X_i = X[,i]
      for( j in (i+1):m ){
        X_j = X[,j]
        r = log(((X_i+Pseudo_Count_used) / (X_j+Pseudo_Count_used)))
        ret[j] = sd(r)
      }
      return(ret)
    }
    ratio_matrix = foreach(i=1:(m-1), .options.RNG=1234,.combine = 'cbind') %dorng% {
      SDs_for_fixed_i(i)
    }
    stopCluster(cl)
    ratio_matrix = cbind(ratio_matrix,rep(NA,nrow(ratio_matrix)))
    for( i in 1:(m-1) ){
      for( j in (i+1):m ){
        #print(paste0("i ",i, "j ",j))
        ratio_matrix[i,j] = ratio_matrix[j,i]
      }
    }
    
  }else{
    if(is.null(Precomputed_scores)){
      m = ncol(Previous_Reference_Selection_Object$ratio_matrix)
      ratio_matrix = Previous_Reference_Selection_Object$ratio_matrix  
    }else{
      m = dim(X)[2]
    }
    
  }
  
  
  
  # for the default, null, all taxa can serve as reference
  if(is.null(select_from)){
    select_from = 1:m
  }
  
  
  #compute the different median SD scores
  if(is.null(Precomputed_scores))
    scores         = apply(ratio_matrix,2,function(x){median(x,na.rm = T)}) 
  else
    scores = Precomputed_scores         
  #compute prevalance of taxa
  prevalence_mat = 1*(X > 0)
  mean_prevalance = apply(prevalence_mat,2,mean)
  
  #filter out columns that the user asked to not be considered
  filter_cols = 1:m
  filter_cols = filter_cols[which(filter_cols %in% select_from)]
  scores = scores[filter_cols]
  sorted_columns = order(scores)
  original_ind = (filter_cols)[sorted_columns]
  sorted_scores = scores[sorted_columns]
  
  #sort the prevalence matrix, by the order of medianSD scores.  
  sorted_prevalence_mat = prevalence_mat[ , original_ind]
  sorted_X_mat = X[,original_ind]
  
  # compute the cummulative prevalence and abundance of samples, for taxa from the lowest medianSD score up to any column on the matrix.
  cummulative_sorted_prevalence_mat = sorted_prevalence_mat
  cummulative_sorted_X_mat = sorted_X_mat
  for(i in 1:nrow(cummulative_sorted_prevalence_mat)){
    cummulative_sorted_prevalence_mat[i,] = cummax(cummulative_sorted_prevalence_mat[i,])
    cummulative_sorted_X_mat[i,] = cumsum(cummulative_sorted_X_mat[i,])
  }
  
  # minimum abundance across taxa, for each possible set of references, by including one additional taxon at a time in the reference set:
  mean_prevalence_over_the_sorted = as.numeric(apply(cummulative_sorted_prevalence_mat,2,mean))
  min_abundance_over_the_sorted = as.numeric(apply(cummulative_sorted_X_mat,2,min))
  
  #find the possible cut point by required abundance (these will be the fall back, if there are no point to find by requested threshold:
  
  possible_cut_points = which(min_abundance_over_the_sorted>=minimal_TA & min_abundance_over_the_sorted<=maximal_TA)
  # if no point are found, go the the firt place (effectivly increase the threshold) the the minimal number of counts is found in each subject
  if(length(possible_cut_points) == 0)
    possible_cut_points = min(which(min_abundance_over_the_sorted>=minimal_TA))
  
  # if not possible, find a cut point which is possible, with a lower number of counts
  if(length(possible_cut_points) == 1 & is.infinite(possible_cut_points[1]))
    possible_cut_points = min(which(min_abundance_over_the_sorted>=1))
  
  # out of all possible points that meet demands for abundance, we first consider those that also meet the demand for median SD:
  scores_possible_cut_points = sorted_scores[possible_cut_points] 
  
  threshold_to_use = median_SD_threshold
  
  possible_cut_points_above_threshold = which(scores_possible_cut_points >= threshold_to_use)
  if(length(possible_cut_points_above_threshold)>0) # if there are points which also meet the medianSD demand, we take the smallest one
    cut_point_selected = possible_cut_points[min(possible_cut_points_above_threshold)]
  else
    cut_point_selected = max(possible_cut_points) #we fall back to requiring only a minimal number of counts
  
  #the selected references
  selected_references = original_ind[1:cut_point_selected]
  
  #the minimum number of counts observed for reference taxa in a subject
  selected_MinAbundance = min_abundance_over_the_sorted[cut_point_selected] 
  
  #return results:
  ret = list()
  ret$selected_references = selected_references
  ret$mean_prevalence_over_the_sorted = mean_prevalence_over_the_sorted
  ret$min_abundance_over_the_sorted = min_abundance_over_the_sorted
  if(is.null(Precomputed_scores))
    ret$ratio_matrix = ratio_matrix
  ret$scores = scores
  ret$selected_MinAbundance = selected_MinAbundance
  ret$median_SD_threshold = median_SD_threshold
  ret$minimal_TA = minimal_TA
  ret$maximal_TA = maximal_TA
  class(ret) = dacomp:::CLASS.LABEL.REFERENCE_SELECTION_OBJECT
  return(ret)
}

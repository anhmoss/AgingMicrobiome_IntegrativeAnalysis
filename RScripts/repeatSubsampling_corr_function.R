#repeat subsampling with Correlation stat test
#input: number of iterations/permutations, desired subsample size, lognorm file, metadata, correlation method
#output: list with three elements: subsample size, mean of significant taxa sum, standard deviation of significant taxa sum


repeatSubsampling_corr_function = function(maxiteration, subsamplesize, lognormFile, metadata, corrMethod) {
  
  fullIteration_results = NULL
  
  for(i in 1:maxiteration) {
    
    results =  subsampling_corr_function(lognormFile, metadata, subsamplesize, corrMethod)
    fullIteration_results[i] = results[length(results)]
    
  }
  
  subsample_mean = mean(as.numeric(fullIteration_results))
  subsample_sd = sd(as.numeric(fullIteration_results))
  
  repeated_results = list("SampleSize" = subsamplesize,
                          "SubSample_Mean"= subsample_mean, 
                          "SubSample_SD" = subsample_sd)
  return(repeated_results)
  
}

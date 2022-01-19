#power analysis: subsamplings and repeating x 50

#____________________________splitting the steps____________

# stat test: LM
subsampling_LM_function = function(lognormFile, metadata, subsamplesize) {
  file = lognormFile
  file$Age = metadata 
  
  lognorm_subsample = slice_sample(file, n=subsamplesize) 
  
  lognorm_subsample_taxa = lognorm_subsample[,-(ncol(lognorm_subsample))] 
  lognorm_subsample_taxa= lognorm_subsample_taxa[, colSums(lognorm_subsample_taxa)!=0]
  # lognorm_subsample_taxa= lognorm_subsample_taxa[, grep("g__", colnames(lognorm_subsample_taxa))]
  
  meta_subsample = lognorm_subsample$Age 
  
  statsResults_FDR = stat_simpleLM_function(lognorm_subsample_taxa, meta_subsample)
  FDR_pvalues = statsResults_FDR[,3]
  totalSignifTaxa=sum(as.numeric(FDR_pvalues) < 0.05)
  FDR_and_sigifTaxaSum = c(FDR_pvalues, "totalSignifTaxa" =totalSignifTaxa)
  return(FDR_and_sigifTaxaSum)
}


repeatSubsampling_LM_function = function(maxiteration, subsamplesize, lognormFile, metadata) {
  
  fullIteration_results = NULL
  
  for(i in 1:maxiteration) {
    
    results =  subsampling_LM_function(lognormFile, metadata, subsamplesize)
    fullIteration_results[i] = results[length(results)]
    
  }

  subsample_mean = mean(as.numeric(fullIteration_results))
  subsample_sd = sd(as.numeric(fullIteration_results))

  repeated_results = list("SampleSize" = subsamplesize,
                        "SubSample_Mean"= subsample_mean, 
                        "SubSample_SD" = subsample_sd)
  return(repeated_results)
  
}

# stat test: kendall corr

subsampling_corr_function = function(lognormFile, metadata, subsamplesize, corrMethod) {
  
  file = lognormFile
  file$Age = metadata 
  
  lognorm_subsample = slice_sample(file, n=subsamplesize) 
  
  lognorm_subsample_taxa = lognorm_subsample[,-(ncol(lognorm_subsample))] 
  lognorm_subsample_taxa= lognorm_subsample_taxa[, colSums(lognorm_subsample_taxa)!=0]
  # lognorm_subsample_taxa= lognorm_subsample_taxa[, grep("g__", colnames(lognorm_subsample_taxa))]
  
  meta_subsample = lognorm_subsample$Age 
  
  statsResults_FDR = stat_correlation_function(lognorm_subsample_taxa, meta_subsample, corrMethod)
  FDR_pvalues = statsResults_FDR[,3]
  totalSignifTaxa=sum(as.numeric(FDR_pvalues) < 0.05)
  FDR_and_sigifTaxaSum = c(FDR_pvalues, "totalSignifTaxa" =totalSignifTaxa)
  return(FDR_and_sigifTaxaSum)
  
}


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


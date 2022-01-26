#random subsetting of dataset, then runs stat test
# stat test: correlation 
#input: lognorm counts, metadata, and desired sub-sampling size
#output: vector of FDR pvals, with last value the total number of FDR pvals < 0.05

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

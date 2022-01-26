#input: counts, meta
#output: matrix of wilcox.test test results: wilcox.test t-value (for direction), and pval

stat_wilcox_function = function(counts, meta) {
  
  statResults = matrix(nrow=ncol(counts), ncol=4)
  simple_wilcoxon = NULL

  
  #run wilcox.test and obtain stat and pval
  for(i in 1:ncol(counts))
  {
    simple_wilcoxon = wilcox.test(counts[,i] ~ meta, conf.int=TRUE)
    
    
    if(class(simple_wilcoxon)=="try-error"){
      statResults[i,]  =NA
    }
    
    statResults[i,1] = simple_wilcoxon$estimate
    statResults[i,2] = simple_wilcoxon$p.value
    
    if(statResults[i,1] < 0){
      statResults[i,4] = "negative"
    } else {statResults[i,4] = "positive"}
    
  }
  statResults[,3] = p.adjust(statResults[,2], method="BH")
  colnames(statResults) = c("stats", "pval", "pvalFDR", "enrichment")
  row.names(statResults) = colnames(counts)
  return(statResults)
  
} 


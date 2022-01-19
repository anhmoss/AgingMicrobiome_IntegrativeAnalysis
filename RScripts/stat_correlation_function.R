## obtaining lm results pval and estimate from corr test: OTU vs age 

stat_correlation_function= function(counts, age, corrMethod) {
  
  statResults = matrix(nrow=ncol(counts), ncol=4)
  simple_corrTest = NULL
  # taxadf = data.frame(counts)
  #run lm and obtain stat and pval
  for(i in 1:ncol(counts))
  {
    simple_corrTest = cor.test(age, counts[,i], method=corrMethod)
    statResults[i,1] = simple_corrTest$estimate
    statResults[i,2] = simple_corrTest$p.value
    
    #account for pval direction based on stat  
    if(statResults[i,1] < 0){
      statResults[i,4] = "negative"
    } else {statResults[i,4] = "positive"}
    
  }
  
  statResults[,3] = p.adjust(statResults[,2], method="BH")
  colnames(statResults) = c("stats", "pval", "pvalFDR", "enrichment")
  row.names(statResults) = colnames(counts)
  return(statResults)
  
} 
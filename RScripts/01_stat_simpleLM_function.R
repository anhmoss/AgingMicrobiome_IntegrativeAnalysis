#input: counts, age
#output: matrix of LM test results: lm t-value (for direction), and pval

stat_simpleLM_function = function(counts, age) {
  
  statResults = matrix(nrow=ncol(counts), ncol=4)
  simple_lm = NULL
  # taxadf = as.data.frame(counts)
  
  
  #run lm and obtain stat and pval
  for(i in 1:ncol(counts))
  {
    simple_lm = lm(counts[,i] ~ age)
    
    
    if(class(simple_lm)=="try-error"){
      statResults[i,]  =NA
    }
    
    statResults[i,1] = summary(simple_lm)$coefficients[2,3]
    statResults[i,2] = summary(simple_lm)$coefficients[2,4]

  if(statResults[i,1] < 0){
    statResults[i,4] = "negative"
  } else {statResults[i,4] = "positive"}
    
 }
  statResults[,3] = p.adjust(statResults[,2], method="BH")
  colnames(statResults) = c("stats", "pval", "pvalFDR", "enrichment")
  row.names(statResults) = colnames(counts)
  return(statResults)
  
} 


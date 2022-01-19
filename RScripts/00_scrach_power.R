#power with LM models 


#input: counts: lognorm counts table with only taxa, age, i=number of samples/rows to subsample from indicated cohort: qiime_genus
# power_subsample_function = function(counts, i) {
#   
#   sliced_counts = slice_sample(counts, n=i)
#   qiime_genus_n = sliced_counts[1:725]
#   simple_lm = NULL
#   lm_pvals=NULL
#   for(i in 1:ncol(qiime_genus_n))
#   {
#     simple_lm = lm(qiime_genus_n[,i] ~ sliced_counts$Age)
#     lm_pvals[i] = summary(simple_lm)$coefficients[2,4]
#     
#   }
#   
#   fdr_pvals = p.adjust(lm_pvals, method = "BH")
#   total_signifTaxa = sum(fdr_pvals < 0.05)
#   return(total_signifTaxa)
# }

#input: x= number of times to repeat the calculation; v = subsample size 

repeat_samples_function = function(x,counts,i) {
  repeat_subsamples = NULL
  samplesize = rep(i, length(x))
  
  for(i in 1:x)
  {
    repeat_subsamples[i] = power_subsample_function(counts,i)
  }
  
  output_df = data.frame(samplesize, repeat_subsamples)
  colnames(output_df) = c("SampleSize", "SignifTaxa")
  return(output_df)
}

mySeq= c(20, 40, 60)

myseq_test = list()

for(i in mySeq){
  myseq_test[[i]] = repeat_samples_function(50, lognorm_agp, mySeq)
}


##my test 12/20/21
manyiterations_matrix = matrix(nrow=length(testoutput_matrix), ncol=6)

testoutput_matrix = head(stats_fdrs_list[[11]][,3])
sum(as.numeric(testoutput_matrix) < 0.05)
testvalue = 5
testoutput_matrix  = c(testoutput_matrix, "sum" = testvalue)
manyiterations_matrix[,6] = testoutput_matrix

for(i in 1:5){
  
manyiterations_matrix[,i] = testoutput_matrix }
rownames(manyiterations_matrix) = names(testoutput_matrix)
names(testoutput_matrix)

dim(testoutput_matrix)

testsubset = slice_sample(lognorm_agp, n=10)
dim(testsubset)
testsubset$Age
testsubset_b = testsubset[,-(ncol(testsubset))]
colnames(testsubset_b)
dim(testsubset_b)

powerRarefaction_function = function(maxiteration, subsamplesize, 
                                     lognormFile, metadata, corrMethod) {
  
  file = lognormFile
  file$Age = metadata 
  lognorm_subsample = slice_sample(file, n=subsamplesize) 
  lognorm_subsample_taxa = lognorm_subsample[,-(ncol(lognorm_subsample))] 
  # lognorm_subsample_taxa= lognorm_subsample_taxa[, colSums(lognorm_subsample_taxa)!=0]
  # lognorm_subsample_taxa= lognorm_subsample_taxa[, grep("g__", colnames(lognorm_subsample_taxa))]

  meta_subsample = lognorm_subsample$Age 
  
  fullIteration_matrix = matrix(nrow = ncol(lognorm_subsample_taxa)+1, ncol=maxiteration) 
  
  for(i in 1:maxiteration) {
    
    statsResults_FDR = stat_correlation_function(lognorm_subsample_taxa, meta_subsample, corrMethod)
    FDR_pvalues = statsResults_FDR[,3]
    totalSignifTaxa=sum(as.numeric(FDR_pvalues) < 0.05)
    FDR_and_sigifTaxaSum = c(FDR_pvalues, "totalSignifTaxa" =totalSignifTaxa)
    fullIteration_matrix[,i] = FDR_and_sigifTaxaSum
    
  }
  
  return(fullIteration_matrix)
}

testmorgan = lognorm_morgan
testmorgan_taxa = testmorgan[, grep("g__", colnames(testmorgan))]
testmorgan_taxa = testmorgan_taxa[, colSums(testmorgan_taxa) != 0]
testmorgan_age = testmorgan$Age

testmorgan_powerfunction = powerRarefaction_function(maxiteration = 50,
                                                     subsamplesize = nrow(lognorm_file),
                                                     lognormFile = lognorm_file, 
                                                     metadata = age_meta,
                                                     corrMethod="kendall")

dim(testmorgan_powerfunction)
names(testmorgan_powerfunction) #need to add taxa names
testmorgan_powerfunction[165,]
testmorgan_powerfunction[nrow(testmorgan_powerfunction),]


##testing from main code body: 
#coutning signif taxa after fdr...this is corr test results...
#corrected to include only genus level taxa/genera and no cols ==0
totalSignifTaxa = matrix(nrow = 11, ncol = 2)

for(i in 1:11){
  totalSignifTaxa[i,1] = coh_l[i]
  totalSignifTaxa[i,2] = sum(as.numeric(stats_fdrs_list[[i]][,3]) < 0.05)
}

##lm
totalSignifTaxa_lm = matrix(nrow = 11, ncol = 2)

for(i in 1:11){
  totalSignifTaxa_lm[i,1] = coh_l[i]
  totalSignifTaxa_lm[i,2] = sum(as.numeric(lm_results_list[[i]][,3]) < 0.05)
}

##power analysis results...with full sample sizes to test..
test_power_totalsigniftaxa_kencorr = matrix(nrow = 11, ncol = 2)

for(i in 1:11){
  test_power_totalsigniftaxa_kencorr[i,1]=coh_l[i]
  test_power_totalsigniftaxa_kencorr[i,2]= test_powerfunction_list[[i]][length(test_powerfunction_list[[i]])]
}

# table(as.numeric(stats_fdrs_list[[11]][,3]) < 0.05) #agp: 39 signif taxa..75 if 'as numeric'
# stats_fdrs_list[stats_fdrs_list[[11]][,3] < 0.05]
# stats_fdrs_list[[11]][,3] < 0.05
# 
# testing_statoutput_df_agp = as.data.frame(stats_fdrs_list[[11]])
# typeof(testing_statoutput_df_agp$pvalFDR)
# sum(as.numeric(testing_statoutput_df_agp$pvalFDR) < 0.05) #39 vs 75
# testing_statoutput_df_agp$pvalFDR[as.numeric(testing_statoutput_df_agp$pvalFDR) < 0.05]

#simple plot comparing pvals from kendall corr and lm
plot(as.numeric(stats_fdrs_list[[1]][,2]), as.numeric(lm_results_list[[1]][,2]),
     main="Non-parametric vs Parametric P-values (Unadjusted)",
     xlab="Non-parametric (Kendall correlation) P-values",
     ylab="Parametric (linear regression) P-values")
all.equal(rownames(stats_fdrs_list[[1]]), rownames(lm_results_list[[1]]))
length(stats_fdrs_list[[1]][,2])
length(lm_results_list[[1]][,2])

for(i in 1:length(vector))
  
  
  ### troubleshooting the powerRarefaction funtion...got it now though
  #needed to also exclude zero columns after the random selection/subsetting
  
  ##
  
  
  testpower_baxter_334 = powerRarefaction_function(maxiteration = 2,
                                               subsamplesize = 40,
                                               lognormFile = testbaxter,
                                               metadata = lognorm_baxter$Age,
                                               corrMethod = "kendall")

all.equal(testpower_baxter_334, testpower_baxter_333)
power_testlist = list()
power_testlist = list(testpower_baxter_333, testpower_baxter_334)
names(power_testlist) = c(40,41)

power_testlist[[2]][nrow(power_testlist[[2]]),]
dim(testbaxter)


testbaxter = lognorm_baxter
testbaxter = testbaxter[,grep("g__", colnames(testbaxter))]
testbaxter = testbaxter[,colSums(testbaxter) != 0]
dim(testbaxter)

dim(testpower_baxter)


all.equal(testpower_baxter[length(testpower_baxter)], 
          test_powerfunction_list[[4]][length(test_powerfunction_list[[4]])])

testpower_LM_baxter = powerRarefaction_LM_function(maxiteration = 1, 
                             subsamplesize = 490,
                             lognormFile = testbaxter,
                             metadata = lognorm_baxter$Age)
t.test(as.numeric(powerResults_subsample100_lm[[4]]), as.numeric(testpower_LM_baxter))
wilcox.test(as.numeric(powerResults_subsample100_lm[[4]]), as.numeric(testpower_LM_baxter))

#exclude the counts at the end
t.test((as.numeric(powerResults_subsample100_lm[[4]][-nrow(powerResults_subsample100_lm[[4]])])),
as.numeric(testpower_LM_baxter[-nrow(testpower_LM_baxter)]))

wilcox.test((as.numeric(powerResults_subsample100_lm[[4]][-nrow(powerResults_subsample100_lm[[4]])])),
       as.numeric(testpower_LM_baxter[-nrow(testpower_LM_baxter)]))

boxplot(t.test((as.numeric(powerResults_subsample100_lm[[4]][-nrow(powerResults_subsample100_lm[[4]])])),
               as.numeric(testpower_LM_baxter[-nrow(testpower_LM_baxter)]))$p.value,
        
        wilcox.test((as.numeric(powerResults_subsample100_lm[[4]][-nrow(powerResults_subsample100_lm[[4]])])),
                    as.numeric(testpower_LM_baxter[-nrow(testpower_LM_baxter)]))$p.value)



# fullIteration_matrix_test = matrix(nrow = (ncol(testbaxter)+1), ncol = 10)
# for(i in 1:10) {
#   
#   
#   statsResults_FDR = stat_correlation_function(testbaxter, age = lognorm_baxter$Age, "kendall")
#   FDR_pvalues = statsResults_FDR[,3]
#   totalSignifTaxa=sum(as.numeric(FDR_pvalues) < 0.05)
#   FDR_and_sigifTaxaSum = c(FDR_pvalues, "totalSignifTaxa" =totalSignifTaxa)
#   fullIteration_matrix_test[,i] = FDR_and_sigifTaxaSum
#   
# }

powerResults_subsample100

subsample100_totalCounts_list = list()

output_matrix = matrix(nrow = 11, ncol = 3)
  
for(n in 1:11){
  
output_matrix[n,1] = mean(as.numeric(powerResults_subsample100[[n]][nrow(powerResults_subsample100[[n]]),]))
output_matrix[n,2] = sd(as.numeric(powerResults_subsample100[[n]][nrow(powerResults_subsample100[[n]]),]))

subsample = 100
output_matrix[n,3] = subsample

}

colnames(output_matrix) = c("mean", "sd", "samplesize")
rownames(output_matrix) = coh_l

## LM______________subsample 100, iterations: 50
output_matrix_LM = matrix(nrow = 11, ncol = 3)

for(n in 1:11){
  
  output_matrix_LM[n,1] = mean(as.numeric(powerResults_subsample100_lm[[n]][nrow(powerResults_subsample100_lm[[n]]),]))
  output_matrix_LM[n,2] = sd(as.numeric(powerResults_subsample100_lm[[n]][nrow(powerResults_subsample100_lm[[n]]),]))
  
  subsample = 100
  output_matrix_LM[n,3] = subsample
  
}

colnames(output_matrix_LM) = c("mean", "sd", "samplesize")
rownames(output_matrix_LM) = coh_l



##compare the outputs of subsample 100 btwn kendall corr and LM per cohort

for(i in 1:11){
jpeg(paste0("/Users/anhil/Desktop/powerboxplots/",coh_l[i], "_run50.jpeg"))

stattest_results=  t.test(as.numeric(powerResults_subsample100[[i]][-nrow(powerResults_subsample100[[i]]),50]),
         as.numeric(powerResults_subsample100_lm[[i]][-nrow(powerResults_subsample100_lm[[i]]),50]))
    
boxplot(as.numeric(powerResults_subsample100[[i]][-nrow(powerResults_subsample100[[i]]),50]),
as.numeric(powerResults_subsample100_lm[[i]][-nrow(powerResults_subsample100_lm[[i]]),50]),
main=paste("P-values from Power Analysis,Subsample=100,Iterations=50\n
P-value=", format(stattest_results$p.value, digits=3), "Cohort:", coh_l[i]),
names=c("Non-parametric(Kendall Corr.)", "Parametric(LM)")) 
dev.off() 
}

##testing the above...

boxplot(as.numeric(powerResults_subsample100[[1]][-nrow(powerResults_subsample100[[1]]),]),
        as.numeric(powerResults_subsample100_lm[[1]][-nrow(powerResults_subsample100_lm[[1]]),]),
        main="P-values from Rarefaction Analysis, Subsample=100, Iterations=50\n
P-value < 2.2e-16 (T-test) ",
        names=c("Non-parametric(Kendall Corr.)", "Parametric(LM)")) 

##___________________________________

t.test(as.numeric(powerResults_subsample100[[1]][-nrow(powerResults_subsample100[[1]]),]),
        as.numeric(powerResults_subsample100_lm[[1]][-nrow(powerResults_subsample100_lm[[1]]),]))

wilcox.test(as.numeric(powerResults_subsample100[[1]][-nrow(powerResults_subsample100[[1]]),]),
        as.numeric(powerResults_subsample100_lm[[1]][-nrow(powerResults_subsample100_lm[[1]]),]))


testpower_merge_df = merge(output_matrix, output_matrix_LM, by=0)
boxplot(testpower_merge_df$mean.x, testpower_merge_df$mean.y,
        main="Comparison of Average Number of Significant Taxa
        at Subsample Size 100 (50 iterations)",
        names=c("Non-parametric", "Parametric"))
names(testpower_merge_df) = c("Cohorts", "Mean(NP)", "SD(NP)", "SampleSize",
                              "Mean(P)", "SD(P)", "SampleSize")

write.table(testpower_merge_df, "PowerAnalysis_LMvsCorr_subsample100.txt", quote = F, row.names = F)
t.test(testpower_merge_df$`Mean(NP)`, testpower_merge_df$`Mean(P)`)



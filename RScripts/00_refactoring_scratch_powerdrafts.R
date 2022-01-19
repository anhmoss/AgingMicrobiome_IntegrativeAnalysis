powercheck_agp = lognorm_agp
powercheck_agp_taxa = powercheck_agp[, grep("g__", colnames(powercheck_agp))]
powercheck_agp_taxa = powercheck_agp_taxa[, colSums(powercheck_agp_taxa) != 0]
powercheck_agp_age = powercheck_agp$Age

dim(powercheck_agp)
dim(powercheck_agp_taxa)

testing_powerFunction_LM_agp = powerRarefaction_LM_function_2(maxiteration = 5,
                                                            subsamplesize = 700,
                                                            lognormFile = powercheck_agp_taxa,
                                                            metadata = powercheck_agp_age)

testing_powerFunction_LM_agp[nrow(testing_powerFunction_LM_agp),]

agp_testing = subsampling_corr_function_function( subsamplesize = 700,
                      lognormFile = powercheck_agp_taxa,
                      metadata = powercheck_agp_age)
agp_testing[length(agp_testing)]

agp_testing_2_corr = repeatSubsampling_corr_function(maxiteration = 50,
                  subsamplesize = 700,
                  lognormFile = powercheck_agp_taxa,
                  metadata = powercheck_agp_age, corrMethod = "kendall")

agp_testing_2_lm = repeatSubsampling_LM_function(maxiteration = 50,
                                                subsamplesize = 700,
                                                lognormFile = powercheck_agp_taxa,
                                                metadata = powercheck_agp_age)

agp_testing_3_lm = repeatSubsampling_LM_function(maxiteration = 50,
                                                 subsamplesize = 100,
                                                 lognormFile = powercheck_agp_taxa,
                                                 metadata = powercheck_agp_age)
agp_testing_4_lm = repeatSubsampling_LM_function(maxiteration = 50,
                                                 subsamplesize = 700,
                                                 lognormFile = powercheck_agp_taxa,
                                                 metadata = powercheck_agp_age)
agp_testing_3_corr = repeatSubsampling_corr_function(maxiteration = 50,
                                                     subsamplesize = 100,
                                                     lognormFile = powercheck_agp_taxa,
                                                     metadata = powercheck_agp_age, corrMethod = "kendall")
agp_testing_4_corr = repeatSubsampling_corr_function(maxiteration = 50,
                                                     subsamplesize = 100,
                                                     lognormFile = powercheck_agp_taxa,
                                                     metadata = powercheck_agp_age, corrMethod = "kendall")
agp_testing_5_corr = repeatSubsampling_corr_function(maxiteration = 50,
                                                     subsamplesize = 100,
                                                     lognormFile = powercheck_agp_taxa,
                                                     metadata = powercheck_agp_age, corrMethod = "kendall")
agp_testing_6_corr = repeatSubsampling_corr_function(maxiteration = 50,
                                                     subsamplesize = 100,
                                                     lognormFile = powercheck_agp_taxa,
                                                     metadata = powercheck_agp_age, corrMethod = "kendall")
agp_testing_7_corr = repeatSubsampling_corr_function(maxiteration = 50,
                                                     subsamplesize = 700,
                                                     lognormFile = powercheck_agp_taxa,
                                                     metadata = powercheck_agp_age, corrMethod = "kendall")


as.data.frame(powercheck_agp_taxa)

short_seqs_agp = c(10,20,30)
long_seqs_agp= c(100,200,700)


testpower_list_agp = list()

for(i in 1:length(long_seqs_agp)){
  
    power_results_LM = repeatSubsampling_LM_function(maxiteration = 5,
                                                              subsamplesize = long_seqs_agp[i],
                                                              lognormFile = powercheck_agp_taxa,
                                                             metadata = powercheck_agp_age)
    testpower_list_agp[[i]]=power_results_LM
  
  
}

names(testpower_list_agp) = long_seqs_agp
testlistagp_df = as.data.frame(testpower_list_agp)

plot(c(testpower_list_agp$`100`$SampleSize, testpower_list_agp$`200`$SampleSize, testpower_list_agp$`700`$SampleSize),
c(testpower_list_agp$`100`$SubSample_Mean, testpower_list_agp$`200`$SubSample_Mean, testpower_list_agp$`700`$SubSample_Mean))

all_samplesize = NULL
all_myMean = NULL

for(k in 1:length(testpower_list_agp)){
  all_samplesize[k] = testpower_list_agp[[k]][1]
  all_myMean[k] = testpower_list_agp[[k]][2]

}
  plot(all_samplesize, all_myMean)

##as df
  all_samplesize = vector()
  all_myMean = vector()
  
  
  for(k in 1:length(testpower_list_agp)){
    all_samplesize[k] = unlist(testpower_list_agp[[k]][1])
    all_myMean[k] = unlist(testpower_list_agp[[k]][2])
    
  }
  
  my_df = data.frame(all_samplesize, all_myMean)
  
  plot(my_df$all_samplesize, my_df$all_myMean, main="testdf")


# dim(testpower_list_agp[[1]])
# testpower_list_agp[[1]][nrow(testpower_list_agp[[1]]),]

  
  ########## ______________________________________________________________________________________________________________________
  
  powerRarefaction_LM_function(
    
    maxiteration = 1,
    subsamplesize = 100,
    lognormFile = lognorm_file,
    metadata = age_meta)
  
  mySubsampleIntervals = c(10,20,30,
                           40,48,
                           50,63,
                           70,84,
                           100,129,156,200, 228, 
                           250,300,350,400,450, 490, 
                           500, 600, 700, 803, 835, 
                           900, 1000, 1371)
  test_subsample = c(10,20,30)
  test_subsample[3]
  
  # 30, 48, 63, 84, 129, 156, 228, 490, 803, 835, 1371
  # maxSamplesize 
  
  testing_maxSampleSize = list()
  
  dim(testing_maxSampleSize)
  
  #lognorm and stats 
  for(i in 1:11){
    raw_file = raw_q2021_files[[i]]
    taxa_table = raw_file[,grep("g__",colnames(raw_file))]
    age_meta = raw_file$Age[rowSums(taxa_table)!=0]
    
    lognorm_file = lognorm_function(taxa_table)
    names_short=sapply(strsplit(colnames(lognorm_file),"g__"),"[[",2)
    colnames(lognorm_file)=names_short
    lognorm_file=lognorm_file[,colSums(lognorm_file)!=0]
    
    
    
    ##power analysis testing....
    
    test_subsample = c(10,20,30)
    coh_l
    
    testnames = paste(coh_l[1], "Subsamplesize:", test_subsample[1])
    
    for(k in test_subsample){
      
      power_result_k = powerRarefaction_LM_function( maxiteration = 2,
                                                     subsamplesize = k,
                                                     lognormFile = lognorm_file,
                                                     metadata = age_meta)
      
      names(power_result_k)=paste(coh_l[i], "Subsamplesize:",test_subsample[k])
      testing_maxSampleSize[[i]] = power_result_k
      
    }
    
  }
  
  
  ##------- 
  
  lognormfiles_q2021
  test_subsample = c(10,20,30)
  test_method = matrix(nrow = 11*length(test_subsample), ncol = 3)
  test_method_list = list()
  test_method_list_2 = list()
  
  for(i in 1:11){
    for(j in 1:length(test_subsample)) {
      
      taxa_input = lognormfiles_q2021[[i]][,-ncol(lognormfiles_q2021[[i]])]
      age_input = lognormfiles_q2021[[i]][,ncol(lognormfiles_q2021[[i]])]
      
      # power_corr_output = repeatSubsampling_corr_function(maxiteration = 50,
      #                                                     subsamplesize = 700,
      #                                                     lognormFile = taxa_input,
      #                                                     metadata = age_input,
      #                                                     corrMethod="kendall")
      # 
      # powerResults_subsample700_corr[[i]] = power_corr_output
      
      power_LM_output = repeatSubsampling_LM_function(maxiteration = 30,
                                                      subsamplesize = test_subsample[j],
                                                      lognormFile = taxa_input,
                                                      metadata = age_input)
      
      tested_subset_matrix = unlist(power_LM_output)
      
      test_method_list[i] = tested_subset_matrix
      
    }
  }
  
  
  
  
## sandbox for drafts and testing/ scratch

#import all files from one file input (filepaths, cohort names as columns)
raw_q2021_filespaths = read.table("/Users/anhil/Desktop/q2021_filepaths.txt", header = TRUE, sep = "\t")
raw_q2021_files = list()

for(i in 1:nrow(raw_q2021_filespaths)) {
  raw_q2021_files[[i]] = read.table(raw_q2021_filespaths[i,1], header = TRUE, sep = "\t")
}

names(raw_q2021_files) = raw_q2021_filespaths[,2]

stats_fdrs_list = list()
coh_l = vector()
# corr_output_list = list()

head(stats_fdrs_list[[11]])

#lognorm and stats 
for(i in 1:11){
  raw_file = raw_q2021_files[[i]]
  taxa_table = raw_file[,grep("g__",colnames(raw_file))]
  # meta_index= raw_file[1:3]
  age_meta = raw_file$Age[rowSums(taxa_table)!=0]
  
  # taxa_table= raw_file_genustaxaonly[-(c(1:3))]
  lognorm_file = lognorm_function(taxa_table)
  names_short=sapply(strsplit(colnames(lognorm_file),"g__"),"[[",2)
  colnames(lognorm_file)=names_short
  lognorm_file=lognorm_file[,colSums(lognorm_file)!=0]

  # print(c(which(rowSums(lognorm_file)==0), dim(lognorm_file),
  #         length(age_meta), length(raw_file$Age[rowSums(taxa_table)!=0]))) }

  stat_output = stat_simpleLM_function(counts=lognorm_file, age = age_meta)
  stats_fdrs_list[[i]] = stat_output
  # corr_output = stat_correlation_function(counts=lognorm_file, age = age_meta, corrMethod ="kendall")
  # corr_output_list[[i]] = corr_output
  coh_l[i]=raw_file[1,1]
}

##print out pval comparison as files

for(n in 1:10){
  for(m in (n+1):11){
    jpeg(paste0("/Users/anhil/Desktop/lm_test_3/",coh_l[[n]],"_vs_",coh_l[[m]],"_lm.jpeg"))
    p_compare(stats_fdrs_list[[n]], stats_fdrs_list[[m]],p_col1=2,p_col2=2,indicator1=4,indicator2=4,point_color="black",lab_cutoff=0.005,cor_method="kendall")
    dev.off()
  }
}


#testing and scratch-------------------------
test.index = c(1:3)
test.logfile = lognorm_function(raw_q2021_files[[11]][-(test.index)])
test.morgan = lognorm_ross
test.morgan.genus = test.morgan[,grep("g__", colnames(test.morgan))]
test.morgan.metaindex = test.morgan[(ncol(test.morgan)-1):ncol(test.morgan)]

raw_q2021_files[[1]][3]
testa= raw_q2021_files[[1]]
all.equal(testa[3], testa[,3])
typeof(as.data.frame(testa[-(test.index)]))
dim(testa)

all.equal(test.logfile, lognorm_agp[1:725])
testnames = test.logfile[,grep("g__",colnames(test.logfile))]
testshort = sapply(strsplit(colnames(testnames),"g__"),"[[",2)
colnames(testnames) = testshort
colnames(testnames)
head(colnames(testnames[-test.index]))
typeof(raw_q2021_files[[11]])
test_b = testnames[-test.index]
colnames(test_b[1:3])

head(lognorm_morgan[216:217])
length(kendallcorrResults_morgan$kendalltau %in% corr_output_list[[1]][,1])
##fix the agp file..
#read in files into one matrix
#loop to go through matrix per cohort
### get metadata_index....raw_q2021_files[[i]][1:3]
### lognorm counts
### clean names: testnames = sapply(strsplit(colnames(raw_q2021_files[[1]][1:3]), "g__", "[[",2))
### stat test of interest
### out of loop: plot pairs
### _pcomparison function

test_stats_a = t(head(stats_fdrs_list[[1]]))
test_stats_b = t(head(stats_fdrs_list[[2]]))
test_stats_merge = merge(test_stats_a, test_stats_b, by=0)

mytest_a = stats_fdrs_list[[1]]

### testing out linear  model functions
## obtaining lm results pval and estimate from corr test: OTU vs age 

stat_linearmodel_function = function(counts, age) {
  
  statResults = matrix(nrow=ncol(counts), ncol=5)
  simple_lm = NULL
  taxadf = data.frame(counts)
  #run lm and obtain stat and pval
  for(i in 1:ncol(taxadf))
  {
    simple_lm = try(lm(taxadf[,i] ~ age))
    
    if(class(simple_lm)=="try-error"){
      statResults[i,]  =NA
    }
    
    statResults[i,1] = summary(simple_lm)$coefficients[2,3]
    statResults[i,2] = summary(simple_lm)$coefficients[2,4]
    
    #fdr correction for pvals
    
    statResults[,3] = p.adjust(statResults[,2], method="BH")
    statResults[,4] = -log10(as.numeric(statResults[,2]))
    
    #account for pval direction based on stat  
    if(statResults[i,1] < 0){
      statResults[i,5] = as.numeric(statResults[i,4]) * (-1)
    } else {statResults[i,5] = as.numeric(statResults[i,4])}
    
  }
  
  colnames(statResults) = c("stats", "pval", "pvalFDR", "neglog10pvals", "pvalsEnriched")
  row.names(statResults) = colnames(counts)
  
  return(statResults)
  
} 


gloor_lm = stat_simpleLM_function(lognorm_gloor[1:663], lognorm_gloor$Age)
goodrich_lm = stat_simpleLM_function(lognorm_goodrich[1:439], lognorm_goodrich$Age)
baxter_lm = stat_simpleLM_function(lognorm_baxter[1:503], lognorm_baxter$Age)
agp_lm = stat_simpleLM_function(raw_q2021_files[[11]][-c(1:3)], age=raw_q2021_files[[11]][3])

all.equal(raw_q2021_files[[11]][2], rawCounts_agp_nodups[1])

head(raw_q2021_files[[11]][3])
rm(agp_lm)
summary(lm(lognorm_agp[,1] ~lognorm_agp$Age))

test_p_comp_function(gloor_lm, goodrich_lm)
test_p_comp_function(gloor_lm, baxter_lm)
test_p_comp_function(goodrich_lm, baxter_lm)

test_gloor_lm = stat_simpleLM_function(lognorm_gloor[1:663], lognorm_gloor$Age)


all.equal(stats_fdrs_list[[6]][,5], test_gloor_lm[,5])

head(stats_fdrs_list[[6]],n=3)
head(gloor_lm,n=3)

typeof(gloor_lm)
head(gloor_lm)
taxa_gloor = rownames(gloor_lm)
taxa_goodrich = rownames(goodrich_lm)

taxa_gloor_goodrich = list(taxa_goodrich, taxa_gloor)

which(taxa_goodrich %in% taxa_gloor)
common_taxa_gloor_goodrich = Reduce(intersect,taxa_gloor_goodrich)
common_genera_gloor_goodrich = common_taxa_gloor_goodrich[grep("g__", common_taxa_gloor_goodrich)]

table(duplicated(common_genera_gloor_goodrich))

gloor_lm_common = gloor_lm[rownames(gloor_lm) %in% common_genera_gloor_goodrich,]
goodrich_lm_common = goodrich_lm[rownames(goodrich_lm) %in% common_genera_gloor_goodrich,]

shortnames_genus = sapply(strsplit(rownames(gloor_lm_common),"g__"),"[[",2)

rownames(goodrich_lm_common) =shortnames_genus
rownames(gloor_lm_common)=shortnames_genus


dim(gloor_lm_common)
dim(goodrich_lm_common)
dim(gloor_lm)
dim(goodrich_lm)


# matchedpvals_gloorgoodrich = NULL
# matchedpvals_goodrichgloor = NULL
# 
# for(i in 1:nrow(gloor_lm)){
#    if(rownames(gloor_lm) %in% common_genera_gloor_goodrich){
#      if(gloor_lm[i,5] =="positive"){
#        matchedpvals_gloorgoodrich[i] = gloor_lm[i,4]
#      } else{matchedpvals_gloorgoodrich[i] = gloor_lm[i,4] * (-1)}
#    }
# }
#   
#   for(k in 1:nrow(goodrich_lm)){
#    if(intersect(rownames(goodrich_lm),rownames(gloor_lm))){
#      if(goodrich_lm[k,5] =="positive"){
#        matchedpvals_goodrichgloor[k] = as.numeric(goodrich_lm[i,4])
#      } else{matchedpvals_goodrichgloor[k] = as.numeric(goodrich_lm[i,4]) * (-1)}
#    }
#  }

par(mfrow=c(1,1))
jpeg(paste0("/Users/anhil/Desktop/lm_results/", "goodrichvgloor.jpg"))
plot(as.numeric(gloor_lm_common[,5]),as.numeric(goodrich_lm_common[,5]),
     main=paste("P-value=", format(cor.test(as.numeric(gloor_lm_common[,5]),as.numeric(goodrich_lm_common[,5]), method="kendall")$p.value, digits=3),
                "Cor=", format(cor.test(as.numeric(gloor_lm_common[,5]),as.numeric(goodrich_lm_common[,5]), method="kendall")$estimate, digits=3),
                sep=" ")
)
dev.off()



test_p_comp_function = function(table1, table2){
  
  all_taxaNames = list(rownames(table1), rownames(table2))
  matched_taxaNames = Reduce(intersect, all_taxaNames)
  
  matched_taxaNames_genus = matched_taxaNames[grep("g__", matched_taxaNames)]
  
  table1_common = table1[rownames(table1) %in% matched_taxaNames_genus,]
  table2_common = table2[rownames(table2) %in% matched_taxaNames_genus,]
  
  
  
  plot(as.numeric(table1_common[,5]),as.numeric(table2_common[,5]),
       main=paste("P-value=", format(cor.test(as.numeric(table1_common[,5]),as.numeric(table2_common[,5]), method="kendall")$p.value, digits=3),
                  "Cor=", format(cor.test(as.numeric(table1_common[,5]),as.numeric(table2_common[,5]), method="kendall")$estimate, digits=3),
                  sep=" "))
  
}


test_p_comp_function(stats_fdrs_list[[6]],stats_fdrs_list[[7]])


table_a = stats_fdrs_list[[6]]
table_b = stats_fdrs_list[[7]]

head(rownames(table_a),n=3)
head(rownames(table_b),n=3)

short_a = rownames(table_a[1:3,])
short_b = rownames(table_b[1:3,])

short_list = list(short_a, short_b)

short_test_match = Reduce(intersect,short_list)

short_test_a = as.data.frame(table_a[rownames(table_a) %in% short_test_match,])
short_test_b = as.data.frame(table_b[rownames(table_b) %in% short_test_match,])

plot(as.numeric(short_test_a[,2]), as.numeric(short_test_b[,2]))

#make another vector 

pvalEnriched_shorta = vector()
pvalEnriched_shortb = vector()

for(i in 1:nrow(short_test_a)){
  if(as.numeric(short_test_a[i,2]) < 0){
    pvalEnriched_shorta[i] = as.numeric(short_test_a[i,2])* -1
    
  } else {pvalEnriched_shorta[i] = as.numeric(short_test_a[i,2])}
}

for(j in 1:nrow(short_test_b)){
  if(as.numeric(short_test_b[j,2]) < 0){
    pvalEnriched_shortb[j] = as.numeric(short_test_b[j,2]) * -1
  } else{pvalEnriched_shortb[j] = as.numeric(short_test_b[j,2])}
}

short_test_a$pvalsEnriched = pvalEnriched_shorta
short_test_b$pvalsEnriched = pvalEnriched_shortb

plot(short_test_a$pvalsEnriched, short_test_b$pvalsEnriched)

p_compare(stats_fdrs_list[[6]],stats_fdrs_list[[7]],p_col1=2,p_col2=2,indicator1=4,indicator2=4,point_color="black",lab_cutoff=0.005,cor_method="kendall")


rownames(short_test_a)
##_____________
test_listnames = list(rownames(table_a), rownames(table_b))
match_testlist = Reduce(intersect,test_listnames)

secondtest = which(rownames(table_a) %in% rownames(table_b))
dim(table_a[rownames(table_a) %in% rownames(table_b),])

table(duplicated(match_testlist))

common_table_a = table_a[rownames(table_a) %in% match_testlist,]
common_table_b = table_b[rownames(table_b) %in% match_testlist,]

length(rownames(table_b))
dim(common_table_a)

dim(common_table_b)

###

# raw_morgan_qiime_genus <- read.csv("/Users/anhil/Desktop/qiime2_rawCounts_anh/morgan_genus_qiime_raw.csv", header=TRUE, sep=",")
# raw_zellerfrance16s_qiime_genus <- read.csv("/Users/anhil/Desktop/qiime2_rawCounts_anh/zellerfrance_genus_qiime_raw.csv", header=TRUE, sep=",")
# raw_zellergermany16s_qiime_genus <- read.csv("/Users/anhil/Desktop/qiime2_rawCounts_anh/zellergermany_genus_qiime_raw.csv", header=TRUE, sep=",")
# 
# rawCounts_filePath = read.csv("/Users/anhil/Desktop/qiime2_rawCounts_anh/morgan_genus_qiime_raw.csv", 
#                        "/Users/anhil/Desktop/qiime2_rawCounts_anh/zellerfrance_genus_qiime_raw.csv",
#                        "/Users/anhil/Desktop/qiime2_rawCounts_anh/zellergermany_genus_qiime_raw.csv",
#                        header=TRUE, sep=",")
# rawCounts_filePath_imported = lapply(rawCounts_filePath, read.csv)
# studynames_test = c("MorganRaw", "ZellerFranceRaw", "ZellerGermanyRaw")

# testingimport_file = read.table("/Users/anhil/Desktop/testingimport.txt", header = TRUE, sep = "\t")
# 
# #reading_in_multipleFiles = apply(testingimport_file, 1, read.csv)
# 
# reading_in_files = list()
# 
# for(i in 1:nrow(testingimport_file))
# {
#   
#   reading_in_files[[i]] = read.table(testingimport_file[i,1], header = TRUE, sep = "\t")
#   
# }
# 
# names(reading_in_files) = testingimport_file[,2]
# 
# all.equal(reading_in_files$Morgan, rawCounts_morgan)
# all.equal(reading_in_files$ZellerFrance, rawCounts_zellerfrance)
# all.equal(reading_in_files$ZellerGermany, rawCounts_zellergermany)

##scracth from power analysis for agp: 
##testing

myseq = seq(20,60, by=20)
for(number in myseq){
  print(number)
}

mytest_a = c(20,1,0.5)
mytest_b = c(40, 2, 0.3)

mytestlist = list(mytest_a, mytest_b)

mytestdf = cbind(mytest_a, mytest_b)
mytestdf_byrow = as.matrix(rbind(mytest_a, mytest_b))
plot(mytestdf_byrow[,1], mytestdf_byrow[,2])

repeat_samples_function = function(permute, subsamplesize) {
  repeat_subsamples = NULL
  for(i in 1:permute)
  {
    repeat_subsamples[i] = subsample_function_agp(subsamplesize)
  }
  
  permuted_subsamples_result = c(subsamplesize, mean(repeat_subsamples), sd(repeat_subsamples))
  return(permuted_subsamples_result)
}

# total_power_results = NULL
# 
# for(i in myseq){
#   
# total_power_results=repeat_samples_function(50, i)
# 
# return(total_power_results)  
# }
sapply(myseq, repeat_samples_function(permute = 50, subsamplesize = myseq))
##
mean(agp_repeated_x50_n60)
sd(agp_repeated_x50_n60)


total_test= matrix(nrow=length(myseq), ncol=3)

mydftest = data.frame(nrow=length(myseq), ncol=3)

for(i in 1:length(myseq)){
  
  mydftest[i] = repeat_samples_function(50, myseq[i])
  
}

total_test_df = as.data.frame(total_test)

plot(total_test[,1], total_test[,2])
plot(total_test_df$V1, total_test_df$V2)

### scractch draft code for power functions...#power analysis: subsample and repeat

powerRarefaction_function = function(maxiteration, subsamplesize,
                                     lognormFile, metadata, corrMethod) {
  
  file = lognormFile
  file$Age = metadata
  lognorm_subsample = slice_sample(file, n=subsamplesize)
  
  lognorm_subsample_taxa = lognorm_subsample[,-(ncol(lognorm_subsample))]
  lognorm_subsample_taxa= lognorm_subsample_taxa[, colSums(lognorm_subsample_taxa)!=0]
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


##lm

powerRarefaction_LM_function = function(maxiteration, subsamplesize,
                                        lognormFile, metadata) {
  
  fullIteration_matrix = matrix(nrow = ncol(lognorm_subsample_taxa)+1, ncol=maxiteration)
  
  for(i in 1:maxiteration) {
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
    fullIteration_matrix[,i] = FDR_and_sigifTaxaSum
    
  }
  
  return(fullIteration_matrix)
}


###_____________

powerRarefaction_LM_function_2 = function(maxiteration, subsamplesize, lognormFile, metadata) {
  
  fullIteration_matrix = matrix(nrow = ncol(lognorm_subsample_taxa)+1, ncol=maxiteration)
  file = lognormFile
  file$Age = metadata
  
  for(i in 1:maxiteration) {
    
    lognorm_subsample = slice_sample(file, n=subsamplesize)
    
    lognorm_subsample_taxa = lognorm_subsample[,-(ncol(lognorm_subsample))]
    lognorm_subsample_taxa= lognorm_subsample_taxa[, colSums(lognorm_subsample_taxa)!=0]
    lognorm_subsample_taxa= lognorm_subsample_taxa[, grep("g__", colnames(lognorm_subsample_taxa))]
    
    meta_subsample = lognorm_subsample$Age
    
    statsResults_FDR = stat_simpleLM_function(lognorm_subsample_taxa, meta_subsample)
    FDR_pvalues = statsResults_FDR[,3]
    totalSignifTaxa=sum(as.numeric(FDR_pvalues) < 0.05)
    FDR_and_sigifTaxaSum = c(FDR_pvalues, "totalSignifTaxa" =totalSignifTaxa)
    fullIteration_matrix[,i] = FDR_and_sigifTaxaSum
    
  }
  
  return(fullIteration_matrix)
}

###
##_________________________________________________testing outside the main loop 


taxa_mytest = lognormfiles_q2021[[5]][,-ncol(lognormfiles_q2021[[5]])]

age_mytest = lognormfiles_q2021[[5]][,ncol(lognormfiles_q2021[[5]])]

power_LM_mytest = repeatSubsampling_LM_function(maxiteration = 50,
                                                subsamplesize = 700,
                                                lognormFile = taxa_mytest,
                                                metadata = age_mytest)
dim(taxa_mytest)

lognormfiles_q2021
powerResults_subsample100_corr = list()
powerResults_subsample100_lm = list()

powerResults_subsample700_corr = list()
powerResults_subsample700_lm = list()

powerResults_subsample700_lm_testmethod1 = matrix(nrow=11,ncol=3)
testmethod_listingMatrixOutputs = list()

for(i in 1:11){
  taxa_input = lognormfiles_q2021[[i]][,-ncol(lognormfiles_q2021[[i]])]
  age_input = lognormfiles_q2021[[i]][,ncol(lognormfiles_q2021[[i]])]
  
  
  # power_corr_output = repeatSubsampling_corr_function(maxiteration = 50,
  #                                                     subsamplesize = 700,
  #                                                     lognormFile = taxa_input,
  #                                                     metadata = age_input,
  #                                                     corrMethod="kendall")
  # 
  # powerResults_subsample700_corr[[i]] = power_corr_output
  
  power_LM_output = repeatSubsampling_LM_function(maxiteration = 50,
                                                  subsamplesize = 100,
                                                  lognormFile = taxa_input,
                                                  metadata = age_input)
  
  powerResults_subsample700_lm_testmethod1[i,]=unlist(power_LM_output)
  
}

# testagp_700_lm =repeatSubsampling_LM_function(maxiteration = 50,
#                                               subsamplesize = 100,
#                                               lognormFile = lognormfiles_q2021[[11]][,-ncol(lognormfiles_q2021[[11]])],
#                                               metadata = lognormfiles_q2021[[11]][,ncol(lognormfiles_q2021[[11]])])
# 
# testagp_700_lm_a =repeatSubsampling_LM_function(maxiteration = 100,
#                                               subsamplesize = 200,
#                                               lognormFile = lognormfiles_q2021[[11]][,-ncol(lognormfiles_q2021[[11]])],
#                                               metadata = lognormfiles_q2021[[11]][,ncol(lognormfiles_q2021[[11]])])
# 
# 
# testagp_700_lm_2 =repeatSubsampling_LM_function(maxiteration = 50,
#                                               subsamplesize = 100,
#                                               lognormFile = lognorm_agp[1:725],
#                                               metadata = lognorm_agp$Age)
# testagp_700_lm_3 =repeatSubsampling_LM_function(maxiteration = 50,
#                                                 subsamplesize = 200,
#                                                 lognormFile = lognorm_agp[1:725],
#                                                 metadata = lognorm_agp$Age)
# testoutput_list = matrix(nrow=3, ncol=3)
# myshortlist = c(10,20,30)
# for(i in 1:length(myshortlist)){
#   testoutput_list[i,] = unlist(repeatSubsampling_LM_function(maxiteration = 50,
#                                 subsamplesize = myshortlist[i],
#                                 lognormFile = lognorm_agp[1:725],
#                                 metadata = lognorm_agp$Age))
# }

rownames(powerResults_subsample700_lm_testmethod1) = coh_l
colnames(powerResults_subsample700_lm_testmethod1) = c("Subsamplesize", "Mean", "SD")



powerResults_subsample100_lm[[5]]
names(powerResults_subsample700_corr) = coh_l
names(powerResults_subsample700_lm) = coh_l


tail(lognormfiles_q2021[[11]][,ncol(lognormfiles_q2021[[11]])],n=7)
tail(powercheck_agp_taxa[,1],n=7)

powerResults_subsample100_corr

unlist(powerResults_subsample700_corr[[11]])
powerResults_subsample700_corr$AGP$SubSample_Mean
powerResults_subsample700_lm

##parsing results from list and putting into dataframe
subsample700_mean_corr = NULL
subsample700_mean_lm = NULL

subsample700_sd_corr = NULL
subsample700_sd_lm = NULL

corr_700 = matrix(nrow = 11, ncol = 3)
lm_700 = matrix(nrow = 11, ncol = 3)

colnames(corr_700) = c("SubSampleSize", "Mean(NP)", "SD(NP)")
rownames(corr_700) = coh_l
colnames(lm_700) = c("SubSampleSize", "Mean(P)", "SD(P)")
rownames(lm_700) = coh_l

for(i in 1:11){
  # subsample700_mean_corr = unlist(powerResults_subsample700_corr[[i]][2])
  # subsample700_mean_lm = unlist(powerResults_subsample700_lm[[i]][2])
  # 
  # subsample700_sd_corr = unlist(powerResults_subsample700_corr[[i]][3])
  # subsample700_sd_lm = unlist(powerResults_subsample700_lm[[i]][3])
  # 
  
  corr_700[i,] = unlist(powerResults_subsample700_corr[[i]])
  lm_700[i,] = unlist(powerResults_subsample700_lm[[i]])
  
}

subsample700_merge = merge(corr_700, lm_700, by=0)



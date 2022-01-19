## with only age as variable 

#step 1: reads in files and parses through for genus 
#step 2: run stats test, power analysis

#import all files from one file input (filepaths, cohort names as columns)
raw_q2021_filespaths = read.table("/Users/anhil/Desktop/q2021_filepaths.txt", header = TRUE, sep = "\t")
raw_q2021_files = list()

for(i in 1:nrow(raw_q2021_filespaths)) {
  raw_q2021_files[[i]] = read.table(raw_q2021_filespaths[i,1], header = TRUE, sep = "\t")
}

names(raw_q2021_files) = raw_q2021_filespaths[,2]


###___________________________________________________________________________


stats_fdrs_list = list()
lm_results_list = list()

coh_l = vector()

lognormfiles_q2021 = list()

permanovaResult_cohorts_matrix = matrix(nrow=11, ncol=3)
colnames(permanovaResult_cohorts_matrix) = c("Cohort", "PERMANOVA_Pvalue", "PERMANOVA_R2")

richness_list = list()
shannondiv_list = list()
myPlotColors = c("black", "cyan", "darkgreen", "lightblue",
                 "limegreen", "orange", "purple",
                 "olivedrab", "red", "darkred", "plum")


#lognorm and stats at the genus level
for(i in 1:11){
  raw_file = raw_q2021_files[[i]]
  taxa_table = raw_file[,grep("g__",colnames(raw_file))]
  taxa_table = taxa_table[rowSums(taxa_table)!=0,]
  age_meta = raw_file$Age[rowSums(taxa_table)!=0]
  
  lognorm_file = lognorm_function(taxa_table)
  names_short=sapply(strsplit(colnames(lognorm_file),"g__"),"[[",2)
  colnames(lognorm_file)=names_short
  lognorm_file=lognorm_file[,colSums(lognorm_file)!=0]

  coh_l[i]=raw_file[1,1]

  #exporting lognorm files for indept testing...
  lognorm_file$Age = age_meta
  lognormfiles_q2021[[i]] = lognorm_file
  
  stat_output = stat_correlation_function(counts=lognorm_file, age = age_meta, corrMethod = "kendall")
  stats_fdrs_list[[i]] = stat_output
  
  lm_output = stat_simpleLM_function(counts=lognorm_file, age = age_meta)
  lm_results_list[[i]] = lm_output
  
  #rarefy raw counts and alpha div and richness analysis
  rarefied_1000 = rrarefy(taxa_table,1000)
  rarefied_1000 = rarefied_1000[,colSums(rarefied_1000) != 0]
  
  #for richness
  richness_result = specnumber(rarefied_1000)
  richness_result[which(rowSums(rarefied_1000)<1000)] =NA
  richness_list[[i]] = richness_result
  
  #for shannon diversity
  shannondiv_result = diversity(rarefied_1000,index = "shannon", MARGIN = 1, base = exp(1))
  shannondiv_result[which(rowSums(rarefied_1000)<1000)] =NA
  shannondiv_list[[i]] = shannondiv_result
  
  kcorr_agevrichness = cor.test(age_meta, richness_result, method="kendall")
  kcorr_agevshannondiv = cor.test(age_meta, shannondiv_result, method="kendall")
  
  #plot age v richness
  
  jpeg(paste0("/Users/anhil/Desktop/figures_jan10/", coh_l[[i]], "_Richness_Age_scatterplots.jpeg"))
  if(kcorr_agevrichness$p.value <0.05)
  { plot(age_meta, richness_result,
         xlab= "Age",
         ylab = "Richness",
         main=c(paste(coh_l[[i]]), 
                paste("P-value=", format(kcorr_agevrichness$p.value,digits=3),
                      "R=", format(kcorr_agevrichness$estimate,digits=3),sep=" ")),
         col.main="red",
         xlim = c(0,115),col=myPlotColors[i])} else {plot(age_meta, richness_result,
                                                          xlab= "Age",
                                                          ylab = "Richness",
                                                          main=c(paste(coh_l[[i]]), 
                                                                 paste("P-value=", format(kcorr_agevrichness$p.value,digits=3),
                                                                       "R=", format(kcorr_agevrichness$estimate,digits=3),sep=" ")),
                                                          xlim = c(0,115),col=myPlotColors[i])}
  dev.off()
  
  #plot age v shannon div
  jpeg(paste0("/Users/anhil/Desktop/figures_jan10/", coh_l[[i]], "_ShannonDiv_Age_scatterplots.jpeg"))
  
  if(kcorr_agevshannondiv$p.value <0.05)
  { plot(age_meta, shannondiv_result,
         xlab= "Age",
         ylab = "Shannon Diversity",
         main=c(paste(coh_l[[i]]), 
                paste("P-value=", format(kcorr_agevshannondiv$p.value,digits=3),
                      "R=", format(kcorr_agevshannondiv$estimate,digits=3),sep=" ")),
         col.main="red",
         xlim = c(0,115),col=myPlotColors[i])} else {plot(age_meta, shannondiv_result,
                                                          xlab= "Age",
                                                          ylab = "Shannon Diversity",
                                                          main=c(paste(coh_l[[i]]), 
                                                                 paste("P-value=", format(kcorr_agevshannondiv$p.value,digits=3),
                                                                       "R=", format(kcorr_agevshannondiv$estimate,digits=3),sep=" ")),
                                                          xlim = c(0,115),col=myPlotColors[i])}
  
  dev.off()
  

  ##permanova
  pcoa_indivCohort = capscale(lognorm_file~ 1, distance = "bray")
  kcorr_agevMDS1 = cor.test(age_meta, pcoa_indivCohort$CA$u[,1],method="kendall")
  kcorr_agevMDS2 = cor.test(age_meta, pcoa_indivCohort$CA$u[,2],method="kendall")
  
  #permanova pval and r2 for each cohort, into a matrix
  permanovaResults_indivCohort = adonis(lognorm_file ~ age_meta,permutations=999)
  permanovaPval_indivCohort = permanovaResults_indivCohort$aov.tab$`Pr(>F)`[1]
  permanovaR2_indivCohort = permanovaResults_indivCohort$aov.tab$R2[1]
  
  permanovaResult_cohorts_matrix[i,1] = coh_l[[i]]
  permanovaResult_cohorts_matrix[i,2] = permanovaPval_indivCohort
  permanovaResult_cohorts_matrix[i,3] = permanovaR2_indivCohort 
  
  #plot age v MDS1 and MDS2 scatterplots 
  jpeg(paste0("/Users/anhil/Desktop/figures_jan10/", coh_l[[i]], "_AgevMDS1_scatterplots.jpeg"))
  if(kcorr_agevMDS1$p.value < 0.05) 
  { plot(age_meta, pcoa_indivCohort$CA$u[,1],
         xlab = "Age",
         ylab = "MDS1",
         main=c(paste(coh_l[[i]]),
                paste("P-value=", format(kcorr_agevMDS1$p.value,digits=3),
                      "R=", format(kcorr_agevMDS1$estimate,digits=3),sep=" ")),
         col.main="red", xlim = c(0,115),col=myPlotColors[i])} else {plot(age_meta, pcoa_indivCohort$CA$u[,1],
                                                                          xlab = "Age",
                                                                          ylab = "MDS1",
                                                                          main=c(paste(coh_l[[i]]),
                                                                                 paste("P-value=", format(kcorr_agevMDS1$p.value,digits=3),
                                                                                       "R=", format(kcorr_agevMDS1$estimate,digits=3),sep=" ")),
                                                                          xlim = c(0,115),col=myPlotColors[i])}
  dev.off()
  
  jpeg(paste0("/Users/anhil/Desktop/figures_jan10/", coh_l[[i]], "_AgevMDS2_scatterplots.jpeg"))
  if(kcorr_agevMDS2$p.value < 0.05) 
  { plot(age_meta, pcoa_indivCohort$CA$u[,2],
         xlab = "Age",
         ylab = "MDS2",
         main=c(paste(coh_l[[i]]),
                paste("P-value=", format(kcorr_agevMDS2$p.value,digits=3),
                      "R=", format(kcorr_agevMDS2$estimate,digits=3),sep=" ")),
         col.main="red", xlim = c(0,115),col=myPlotColors[i])} else {plot(age_meta, pcoa_indivCohort$CA$u[,2],
                                                                          xlab = "Age",
                                                                          ylab = "MDS2",
                                                                          main=c(paste(coh_l[[i]]),
                                                                                 paste("P-value=", format(kcorr_agevMDS2$p.value,digits=3),
                                                                                       "R=", format(kcorr_agevMDS2$estimate,digits=3),sep=" ")),
                                                                          xlim = c(0,115),col=myPlotColors[i])}
  dev.off()
}


##_________________________________________________

##print out pval comparison as files

for(n in 1:10){
  for(m in (n+1):11){
    jpeg(paste0("/Users/anhil/Desktop/corr_test_check_4/",coh_l[[n]],"_vs_",coh_l[[m]],"_lm.jpeg"))
    p_compare(stats_fdrs_list[[n]], stats_fdrs_list[[m]],p_col1=2,p_col2=2,indicator1=4,indicator2=4,point_color="black",lab_cutoff=0.005,cor_method="kendall")
    dev.off()
  }
}

## histograms for kendall pvals

stats_fdrs_list
par(mfrow=c(3,4))
for(i in 1:11){
  totalSignifTaxa_sum = sum(as.numeric(stats_fdrs_list[[i]][,3]) < 0.05)
  hist_title = paste0(coh_l[[i]], "\n", "Significant Taxa: ", totalSignifTaxa_sum)
  hist_xaxis = paste0("Kendall P-values")
  hist(as.numeric(stats_fdrs_list[[i]][,2]),
       main=hist_title, 
       xlab = hist_xaxis)
}

## hist organized by row as subgroups: young, middle, full range
par(mfrow=c(3,5))
hist(as.numeric(stats_fdrs_list[[8]][,2]), breaks = 15,
     main = paste0(coh_l[[8]], "\n", " Significant Taxa: ", sum(as.numeric(stats_fdrs_list[[8]][,3]) < 0.05)),
     xlab = "Kendall P-values")
hist(as.numeric(stats_fdrs_list[[9]][,2]), breaks = 15,
     main = paste0(coh_l[[9]], "\n", " Significant Taxa: ", sum(as.numeric(stats_fdrs_list[[9]][,3]) < 0.05)),
     xlab = "Kendall P-values")
hist(as.numeric(stats_fdrs_list[[5]][,2]), breaks = 15,
     main = paste0(coh_l[[5]], "\n", " Significant Taxa: ", sum(as.numeric(stats_fdrs_list[[5]][,3]) < 0.05)),
     xlab = "Kendall P-values")
plot.new()
plot.new()

hist(as.numeric(stats_fdrs_list[[2]][,2]), breaks = 15,
     main = paste0(coh_l[[2]], "\n", " Significant Taxa: ", sum(as.numeric(stats_fdrs_list[[2]][,3]) < 0.05)),
     xlab = "Kendall P-values")
hist(as.numeric(stats_fdrs_list[[3]][,2]), breaks = 15,
     main = paste0(coh_l[[3]], "\n", " Significant Taxa: ", sum(as.numeric(stats_fdrs_list[[3]][,3]) < 0.05)),
     xlab = "Kendall P-values")
hist(as.numeric(stats_fdrs_list[[7]][,2]), breaks = 15,
     main = paste0(coh_l[[7]], "\n", " Significant Taxa: ", sum(as.numeric(stats_fdrs_list[[7]][,3]) < 0.05)),
     xlab = "Kendall P-values")
hist(as.numeric(stats_fdrs_list[[4]][,2]), breaks = 15,
     main = paste0(coh_l[[4]], "\n", " Significant Taxa: ", sum(as.numeric(stats_fdrs_list[[4]][,3]) < 0.05)),
     xlab = "Kendall P-values")
hist(as.numeric(stats_fdrs_list[[10]][,2]), breaks = 15,
     main = paste0(coh_l[[10]], "\n", " Significant Taxa: ", sum(as.numeric(stats_fdrs_list[[10]][,3]) < 0.05)),
     xlab = "Kendall P-values")


hist(as.numeric(stats_fdrs_list[[6]][,2]), breaks = 15,
     main = paste0(coh_l[[6]], "\n", " Significant Taxa: ", sum(as.numeric(stats_fdrs_list[[6]][,3]) < 0.05)),
     xlab = "Kendall P-values")
hist(as.numeric(stats_fdrs_list[[1]][,2]), breaks = 15,
     main = paste0(coh_l[[1]], "\n", " Significant Taxa: ", sum(as.numeric(stats_fdrs_list[[1]][,3]) < 0.05)),
     xlab = "Kendall P-values")
hist(as.numeric(stats_fdrs_list[[11]][,2]), breaks = 15,
     main = paste0(coh_l[[11]], "\n", " Significant Taxa: ", sum(as.numeric(stats_fdrs_list[[11]][,3]) < 0.05)),
     xlab = "Kendall P-values")
plot.new()
plot.new()

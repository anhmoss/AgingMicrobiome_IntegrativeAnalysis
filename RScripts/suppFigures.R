generateSuppFigures = function(myFilePath, resultDirPath) {
  
  #Step 1: reads in counts table, labeled by cohort name
  #input: one file input (filepaths and cohort names as columns)
  #output: list of OTU counts table for all 11 datasets
  
  raw_q2021_filespaths = read.table(myFilePath, header = TRUE, sep = "\t")
  raw_q2021_files = list()
  
  for(i in 1:nrow(raw_q2021_filespaths)) {
    raw_q2021_files[[i]] = read.table(raw_q2021_filespaths[i,1], header = TRUE, sep = "\t")
  }
  names(raw_q2021_files) = raw_q2021_filespaths[,2]
  
  #extracts only stool samples from Morgan cohort for subsequent analysis
  raw_q2021_files[[1]] = subset(raw_q2021_files[[1]], raw_q2021_files[[1]]$stool==1)
  
  pdf(paste0(resultDirPath, "AgeRangeDist_Boxplot.pdf"),onefile = T)
  #plot age range distribution boxplots
  boxplot(
  raw_q2021_files$NogBCN0$Age,
  raw_q2021_files$NogSTK$Age,
  raw_q2021_files$Escobar$Age,
  raw_q2021_files$ZellerFrance$Age,
  raw_q2021_files$ZellerGermany$Age,
  raw_q2021_files$Goodrich$Age,
  raw_q2021_files$Baxter$Age,
  raw_q2021_files$Ross$Age,
  raw_q2021_files$Gloor$Age,
  raw_q2021_files$Morgan$Age,
  raw_q2021_files$AGP$Age,
  cex.axis=0.8,
  names = c(
    "Nog_BCN0", "Nog_STK", "Escobar",
    "ZellerFrance", "ZellerGermany", "Goodrich", "Baxter","Ross", 
    "Gloor", "Morgan","AGP"), main= "Age Range",
  ylab= "Age (Years)")

allstudies_agerange_list <- list(
  "Nog_BCN0" = raw_q2021_files$NogBCN0$Age, 
  "Nog_STK" = raw_q2021_files$NogSTK$Age,
  "Escobar" = raw_q2021_files$Escobar$Age, 
  "ZellerFrance" = raw_q2021_files$ZellerFrance$Age,
  "ZellerGermany" = raw_q2021_files$ZellerGermany$Age, 
  "Goodrich" = raw_q2021_files$Goodrich$Age,
  "Baxter" = raw_q2021_files$Baxter$Age,
  "Ross" = raw_q2021_files$Ross$Age,
  "Gloor" = raw_q2021_files$Gloor$Age,
  "Morgan" = raw_q2021_files$Morgan$Age,
  "AGP"=raw_q2021_files$AGP$Age)
stripchart(allstudies_agerange_list, method = "jitter", pch = 20, vertical= TRUE, add = TRUE)
dev.off()
  
# shannon div comparison on raw vs rarefied counts ______________________________________________________________________________________________
raw_merge_allCounts = full_join(raw_q2021_files[[1]], raw_q2021_files[[2]])

for(i in 3:11){ 
  file = raw_q2021_files[[i]]
  raw_merge_allCounts = full_join(raw_merge_allCounts, file)
}

raw_merge_allCounts[is.na(raw_merge_allCounts)] = 0

raw_merge_allCounts_genusTaxa = raw_merge_allCounts[,grep("g__", colnames(raw_merge_allCounts))]
  
set.seed(316)
#all pooled samples, raw (not normalized), not rarefied
allGenera_raw = raw_merge_allCounts_genusTaxa
allGenera_raw_shannonDiv = diversity(allGenera_raw,index = "shannon", MARGIN = 1, base = exp(1))

#shannon div on rarefied counts
allGenera_raw_rarefied= rrarefy(allGenera_raw, 1000)
allGenera_raw_rarefied_shannonDiv = diversity(allGenera_raw_rarefied,index = "shannon", MARGIN = 1, base = exp(1))
allGenera_raw_rarefied_shannonDiv[which(rowSums(allGenera_raw_rarefied)<1000)] = NA

pdf(paste0(resultDirPath, "ShannonDiversity_RawVSRarefiedCounts.pdf"),onefile = T)
plot(allGenera_raw_shannonDiv,
     allGenera_raw_rarefied_shannonDiv,
     main = "Shannon Diversity on Raw vs Rarefied Counts",
     xlab = "Shannon Diversity (Raw Counts)",
     ylab = "Shannon Diversity (Rarefied Counts)")
 dev.off()
 
  
 #____________________________________________________________________________________________________________________________________________________

  #Step 2: lognormalize counts for each dataset, perform stats test 
  #input: counts table for each dataset (ie each element from the list)
  #output: power analysis (non-parametric and parametric), 
  
  coh_l = vector()
  myLognormFiles = list()
  myPlotColors = c("black", "cyan", "darkgreen", "lightblue",
                   "limegreen", "orange", "purple",
                   "olivedrab", "red", "darkred", "plum")

  for(i in 1:11){
    raw_file = raw_q2021_files[[i]]
    taxa_table = raw_file[,grep("g__",colnames(raw_file))]
    age_meta = raw_file$Age[which(rowSums(taxa_table)!=0)]
    age_fullset = raw_file$Age
    
    lognorm_file = lognorm_function(taxa_table)
    names_short=sapply(strsplit(colnames(lognorm_file),"g__"),"[[",2)
    colnames(lognorm_file)=names_short
    lognorm_file=lognorm_file[,colSums(lognorm_file)!=0]
    
    mylogfile = lognorm_file
    mylogfile$Age = age_meta
    myLognormFiles[[i]] = mylogfile 
    
    coh_l[i]=raw_file[1,1] 
    
    stat_output = stat_correlation_function(counts=lognorm_file, age = age_meta, corrMethod = "kendall")
    stats_fdrs_list[[i]] = stat_output   
    lm_output = stat_simpleLM_function(counts=lognorm_file, age = age_meta)
    lm_results_list[[i]] = lm_output

    #plot histograms for age distribution
      pdf(paste0(resultDirPath, coh_l[i], "_AgeDistribution_Histogram.pdf"),onefile = T,height=10,width=10)
      hist(age_fullset, breaks = 10, xlim = c(0,120),col = "mediumspringgreen",
         main = coh_l[i], xlab = "Age (Years)", cex.main=1.5, cex.lab=1.5) 
      dev.off()  
      
      }
     
     
# plot histograms for LM P-values
for(i in 1:11){
    pdf(paste0(resultDirPath,coh_l[i], "_KendallHist_LM.pdf"),onefile = T,height=10,width=10)
    totalSignifTaxa_sum = sum(as.numeric(lm_results_list[[i]][,3]) < 0.05)
    hist_title = paste0(coh_l[[i]], "\n", "Significant Taxa: ", totalSignifTaxa_sum)
    hist_xaxis = paste0("LM P-values")
    hist(as.numeric(lm_results_list[[i]][,2]),
         main=hist_title, 
         xlab = hist_xaxis)
    dev.off()
  }
  
  
    

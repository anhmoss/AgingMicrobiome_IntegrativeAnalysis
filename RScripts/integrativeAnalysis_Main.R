#Generates plots for integrative analysis
  #input: filepath to input file, filepath to result directory
  #output: Shannon Diversity v Age plots, Richness v Age plots, MDS1 v Age plots, MDS2 v Age plots, Kendall corr p-value histograms, and p-value v p-value plots, Pearson Corr Plot for Bifidobacterium

# The following code is for age as the only variable 

generatePlots = function(myFilePath, resultDirPath) {
  
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
  
  #Step 2: lognormalize counts for each dataset, perform stats test (taxa ~ age), permanova test, alpha diversity analysis at the genus level
  #input: counts table for each dataset (ie each element from the list)
  #output: stats calculations (non-parametric and parametric), permanova results per cohort, alpha diversity results and plots (richness and shannon diversity)
  
  stats_fdrs_list = list()
  lm_results_list = list()
  coh_l = vector()
  permanovaResult_cohorts_matrix = matrix(nrow=11, ncol=3)
  colnames(permanovaResult_cohorts_matrix) = c("Cohort", "PERMANOVA_Pvalue", "PERMANOVA_R2")
  richness_list = list()
  shannondiv_list = list()
  myLognormFiles = list()
  myPlotColors = c("black", "cyan", "darkgreen", "lightblue",
                   "limegreen", "orange", "purple",
                   "olivedrab", "red", "darkred", "plum")
  
  set.seed(123)
  
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
    
    #rarefy raw counts and alpha div and richness analysis
    rarefied_1000 = rrarefy(taxa_table,1000)
    rarefied_1000 = rarefied_1000[,colSums(rarefied_1000) != 0]
    
    richness_result = specnumber(rarefied_1000)
    richness_result[which(rowSums(rarefied_1000)<1000)] =NA
    richness_list[[i]] = richness_result
    
    
    shannondiv_result = diversity(rarefied_1000,index = "shannon", MARGIN = 1, base = exp(1))
    shannondiv_result[which(rowSums(rarefied_1000)<1000)] =NA
    shannondiv_list[[i]] = shannondiv_result
    
    kcorr_agevrichness = cor.test(age_fullset, richness_result, method="kendall")
    kcorr_agevshannondiv = cor.test(age_fullset, shannondiv_result, method="kendall") 
    
    #plot age v richness 
    pdf(paste0(resultDirPath, coh_l[[i]], "_Richness_Age_scatterplots.pdf"),onefile = T,height=10,width=10)
    if(kcorr_agevrichness$p.value <0.05)
    { plot(age_fullset, richness_result,
           xlab= "Age",
           ylab = "Richness", cex.main=1.4,
           main=c(paste(coh_l[[i]]), 
                  paste("P-value=", format(kcorr_agevrichness$p.value,digits=3),
                        "Kendall Tau=", format(kcorr_agevrichness$estimate,digits=3),sep=" ")),
           col.main="red",
           xlim = c(0,115),col=myPlotColors[i])} else {plot(age_fullset, richness_result,
                                                            xlab= "Age",
                                                            ylab = "Richness",cex.main=1.4,
                                                            main=c(paste(coh_l[[i]]), 
                                                                   paste("P-value=", format(kcorr_agevrichness$p.value,digits=3),
                                                                         "Kendall Tau=", format(kcorr_agevrichness$estimate,digits=3),sep=" ")),
                                                            xlim = c(0,115),col=myPlotColors[i])}
    dev.off()
    
    #plot age v shannon div
    pdf(paste0(resultDirPath, coh_l[[i]], "_ShannonDiv_Age_scatterplots.pdf"),onefile = T,height=10,width=10)
    
    if(kcorr_agevshannondiv$p.value <0.05)
    { plot(age_fullset, shannondiv_result,
           xlab= "Age",
           ylab = "Shannon Diversity",cex.main=1.4,
           main=c(paste(coh_l[[i]]), 
                  paste("P-value=", format(kcorr_agevshannondiv$p.value,digits=3),
                        "Kendall Tau=", format(kcorr_agevshannondiv$estimate,digits=3),sep=" ")),
           col.main="red",
           xlim = c(0,115),col=myPlotColors[i])} else {plot(age_fullset, shannondiv_result,
                                                            xlab= "Age",
                                                            ylab = "Shannon Diversity",cex.main=1.4,
                                                            main=c(paste(coh_l[[i]]), 
                                                                   paste("P-value=", format(kcorr_agevshannondiv$p.value,digits=3),
                                                                         "Kendall Tau=", format(kcorr_agevshannondiv$estimate,digits=3),sep=" ")),
                                                            xlim = c(0,115),col=myPlotColors[i])}
    
    dev.off()
    
    
    ##permanova for each cohort
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
    pdf(paste0(resultDirPath, coh_l[[i]], "_AgevMDS1_scatterplots.pdf"),onefile = T,height=10,width=10)
    if(kcorr_agevMDS1$p.value < 0.05) 
    { plot(age_meta, pcoa_indivCohort$CA$u[,1],
           xlab = "Age",
           ylab = "MDS1",cex.main=1.4,
           main=c(paste(coh_l[[i]]),
                  paste("P-value=", format(kcorr_agevMDS1$p.value,digits=3),
                        "Kendall Tau=", format(kcorr_agevMDS1$estimate,digits=3),sep=" ")),
           col.main="red", xlim = c(0,115),col=myPlotColors[i])} else {plot(age_meta, pcoa_indivCohort$CA$u[,1],
                                                                            xlab = "Age",
                                                                            ylab = "MDS1",cex.main=1.4,
                                                                            main=c(paste(coh_l[[i]]),
                                                                                   paste("P-value=", format(kcorr_agevMDS1$p.value,digits=3),
                                                                                         "Kendall Tau=", format(kcorr_agevMDS1$estimate,digits=3),sep=" ")),
                                                                            xlim = c(0,115),col=myPlotColors[i])}
    dev.off()
    
    pdf(paste0(resultDirPath, coh_l[[i]], "_AgevMDS2_scatterplots.pdf"),onefile = T,height=10,width=10)
    if(kcorr_agevMDS2$p.value < 0.05) 
    { plot(age_meta, pcoa_indivCohort$CA$u[,2],
           xlab = "Age",
           ylab = "MDS2",cex.main=1.4,
           main=c(paste(coh_l[[i]]),
                  paste("P-value=", format(kcorr_agevMDS2$p.value,digits=3),
                        "Kendall Tau=", format(kcorr_agevMDS2$estimate,digits=3),sep=" ")),
           col.main="red", xlim = c(0,115),col=myPlotColors[i])} else {plot(age_meta, pcoa_indivCohort$CA$u[,2],
                                                                            xlab = "Age",
                                                                            ylab = "MDS2",cex.main=1.4,
                                                                            main=c(paste(coh_l[[i]]),
                                                                                   paste("P-value=", format(kcorr_agevMDS2$p.value,digits=3),
                                                                                         "Kendall Tau=", format(kcorr_agevMDS2$estimate,digits=3),sep=" ")),
                                                                            xlim = c(0,115),col=myPlotColors[i])}
    dev.off()
  }
  
  write.table(permanovaResult_cohorts_matrix, paste0(resultDirPath,"PermanovaResults_perCohort.txt"), quote = F, row.names = F)
  
  #Step 3: Using the stat results to create p-value v p-value plots and Kendall p-value histograms
  ## p-value v p-value plots
  for(n in 1:10){
    for(m in (n+1):11){
      pdf(paste0(resultDirPath, coh_l[[n]],"_vs_",coh_l[[m]],".pdf"),onefile = T,height=10,width=10)
      p_compare_ver2(stats_fdrs_list[[n]], stats_fdrs_list[[m]],p_col1=2,p_col2=2,indicator1=4,indicator2=4,
                     point_color="black",lab_cutoff=0.005,cor_method="kendall",
                     data1name = coh_l[[n]], data2name = coh_l[[m]])
      dev.off()
    }
  }
  
  ## histograms for kendall pvals
  stats_fdrs_list
  for(i in 1:11){
    pdf(paste0(resultDirPath,coh_l[i], "_KendallHist.pdf"),onefile = T,height=10,width=10)
    totalSignifTaxa_sum = sum(as.numeric(stats_fdrs_list[[i]][,3]) < 0.05)
    hist_title = paste0(coh_l[[i]], "\n", "Significant Taxa: ", totalSignifTaxa_sum)
    hist_xaxis = paste0("Kendall P-values")
    hist(as.numeric(stats_fdrs_list[[i]][,2]),
         main=hist_title, 
         xlab = hist_xaxis)
    dev.off()
  }
  
  ## optional:: hist organized by row as subgroups: young, middle, full range
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
  
  
  #Step 4: checking common taxa across cohorts based on non-parametric statistical tests (kendall corr)
  allTaxaNames = list()
  
  for(i in 1:11){
    allTaxaNames[[i]] = rownames(stats_fdrs_list[[i]])
  }
  
  commonGenusNames_allCohorts = Reduce(intersect, allTaxaNames)
  commonGenusNames_allCohorts_genera = commonGenusNames_allCohorts[-(grep("UCG|uncultured", commonGenusNames_allCohorts))]
  
  # writes text file with all genera common across all cohorts
  write.table(commonGenusNames_allCohorts_genera, paste0(resultDirPath, "commonGenera.txt"), quote = F, row.names = F)
  
  #pulls stats results for only taxa that are common in all 11 cohorts (commonGenusNames_allCohorts_genera)
  commonGenus_narrow_list = list()
  
  for(i in 1:11){
    commonGenus_narrow_list[[i]] = stats_fdrs_list[[i]][rownames(stats_fdrs_list[[i]]) %in% commonGenusNames_allCohorts_genera,]
  }
  
  # corrects for R considering single row matrix as a vector (which drops rownames), so this maintains single-row output for Morgan as a matrix instead of a vector
  commonGenus_narrow_list[[1]][(commonGenus_narrow_list[[1]][,3] < 0.05),,drop=F]
  
  # writes text file w stats matrix (kendall tau, pvalue, fdr pvals, and enrichment) for taxa shared across all datasets
  write.table(commonGenus_narrow_list, paste0(resultDirPath, "commonGenera_statsResults.txt"), quote = F, row.names = F)
  
  #among the common taxa, pull the ones that have fdr adjusted pvals  < 0.05 for each cohort
  commonAndSignificant_list = list()
  
  for(i in 1:11){
    commonAndSignificant_list[[i]] = commonGenus_narrow_list[[i]][as.numeric(commonGenus_narrow_list[[i]][,3]) < 0.05,,drop=F]
    
  }
  names(commonAndSignificant_list) = coh_l
  
  cs_taxa_morgan = rownames(commonAndSignificant_list$Morgan)
  cs_taxa_baxter = rownames(commonAndSignificant_list$Baxter)
  cs_taxa_goodrich = rownames(commonAndSignificant_list$Goodrich)
  cs_taxa_gloor = rownames(commonAndSignificant_list$Gloor)
  cs_taxa_agp = rownames(commonAndSignificant_list$AGP)
  cs_taxa_nogbcn0 = rownames(commonAndSignificant_list$Nog_BCN0)
  
  cs_taxa_allList = list(cs_taxa_morgan, cs_taxa_baxter, cs_taxa_goodrich, cs_taxa_gloor,
                         cs_taxa_agp, cs_taxa_nogbcn0)
  names(cs_taxa_allList) = c("Morgan", "Baxter", "Goodrich", "Gloor", "AGP", "NogBCN0")
  
  intersect(cs_taxa_gloor[cs_taxa_gloor %in% cs_taxa_agp],cs_taxa_gloor[cs_taxa_gloor %in% cs_taxa_nogbcn0])
  #checking if taxa in baxter match w any other ones
  cs_taxa_baxter %in% cs_taxa_goodrich 
  cs_taxa_baxter %in% cs_taxa_gloor 
  cs_taxa_baxter %in% cs_taxa_agp 
  cs_taxa_baxter %in% cs_taxa_nogbcn0 
  
  #goodrich
  sum(cs_taxa_goodrich %in% cs_taxa_gloor) 
  cs_taxa_goodrich[cs_taxa_goodrich %in% cs_taxa_gloor] 
  cs_taxa_goodrich[cs_taxa_goodrich %in% cs_taxa_agp]  
  sum(cs_taxa_goodrich %in% cs_taxa_nogbcn0) 
  
  #gloor
  sum(cs_taxa_gloor %in% cs_taxa_agp) 
  sum(cs_taxa_gloor %in% cs_taxa_nogbcn0) 
  
  #agp
  sum(cs_taxa_agp %in% cs_taxa_nogbcn0) # 3 matches: "Bacteroides" "Howardella"  "Akkermansia"
  
  "Akkermansia" %in% cs_taxa_morgan 
  "Akkermansia" %in% cs_taxa_baxter
  "Akkermansia" %in% cs_taxa_goodrich #true
  "Akkermansia" %in% cs_taxa_gloor #true
  "Akkermansia" %in% cs_taxa_agp #true
  "Akkermansia" %in% cs_taxa_nogbcn0
  
  # Plot pearson correlation (Bifidobacterium ~ Age) for each cohort
  pcorr_bifido = list()
  zeroCounts_perCohort = NULL
  fullSampleSize_perCohort = NULL
  
  for(i in 1:11){
    zeroCounts_perCohort[i] = sum(myLognormFiles[[i]][,-ncol(myLognormFiles[[i]])]$Bifidobacterium==0) 
    fullSampleSize_perCohort[i] = nrow(myLognormFiles[[i]]) 
    pcorr_matrix = matrix(nrow=1, ncol=4)
    colnames(pcorr_matrix) = c("pval", "estimate", "ci_left", "ci_right")
    
    pcorr_bifido_output =  cor.test(myLognormFiles[[i]][,-ncol(myLognormFiles[[i]])]$Bifidobacterium, 
                                    myLognormFiles[[i]][,ncol(myLognormFiles[[i]])], method="pearson")
    pcorr_matrix[,1] = pcorr_bifido_output$p.value
    pcorr_matrix[,2] = pcorr_bifido_output$estimate
    pcorr_matrix[,3] = pcorr_bifido_output$conf.int[1]
    pcorr_matrix[,4] = pcorr_bifido_output$conf.int[2]
    
    pcorr_bifido[[i]] = pcorr_matrix
    
  }   
  names(pcorr_bifido) = coh_l
  pdf(paste0(resultDirPath, "PearsonCorrPlot_Bifidobacterium.pdf"),onefile = T,height=10,width=10) 
  boxplot(
    pcorr_bifido$Nog_BCN0[c(3,4)],
    pcorr_bifido$Nog_STK[c(3,4)],
    pcorr_bifido$Escobar[c(3,4)],
    pcorr_bifido$ZellerFrance[c(3,4)],
    pcorr_bifido$ZellerGermany[c(3,4)],
    pcorr_bifido$Goodrich[c(3,4)],
    pcorr_bifido$Baxter[c(3,4)],
    pcorr_bifido$Ross[c(3,4)],
    pcorr_bifido$Gloor[c(3,4)],
    pcorr_bifido$Morgan[c(3,4)],
    pcorr_bifido$AGP[c(3,4)],
    names = c("NogBCN0\n(n=156)", "NogSTK\n(n=84)", "Escobar\n(n=30)",
              "Zeller\nFrance\n(n=129)", "Zeller\nGermany\n(n=48)", "Goodrich\n(n=835)", "Baxter\n(n=490)", "Ross\n(n=63)",
              "Gloor\n(n=803)", "Morgan\n(n=134)", "AGP\n(n=1,368)"),
    horizontal = T,
    ylim=c(-0.6, 0.6), 
    main="Bifidobacterium",
    xlab="Pearson R Value
  Age ~ OTU Count", 
    las=2, font =2, cex.axis=0.6
    
  )
  text(-0.29,11, "*", cex = 2)
  text(-0.61,9,"*", cex = 2)
  text(-0.23,6, "*", cex = 2)
  abline(v=0, lty=2)
  
  dev.off()  
  
  ## Plot pearson correlation (Faecalibacterium ~ Age) for each cohort
  pcorr_faecal = list()
  zeroCounts_faecal = NULL
  fullSampleSize_perCohort_faecal = NULL
  
  for(i in 1:11){
    
    zeroCounts_faecal[i] = sum(myLognormFiles[[i]][,-ncol(myLognormFiles[[i]])]$Faecalibacterium==0) 
    fullSampleSize_perCohort_faecal[i] = nrow(myLognormFiles[[i]]) 
    pcorr_matrix = matrix(nrow=1, ncol=4)
    colnames(pcorr_matrix) = c("pval", "estimate", "ci_left", "ci_right")
    
    pcorr_faecal_output =  cor.test(myLognormFiles[[i]][,-ncol(myLognormFiles[[i]])]$Faecalibacterium, 
                                    myLognormFiles[[i]][,ncol(myLognormFiles[[i]])], method="pearson")
    pcorr_matrix[,1] = pcorr_faecal_output$p.value
    pcorr_matrix[,2] = pcorr_faecal_output$estimate
    pcorr_matrix[,3] = pcorr_faecal_output$conf.int[1]
    pcorr_matrix[,4] = pcorr_faecal_output$conf.int[2]
    
    pcorr_faecal[[i]] = pcorr_matrix
    
  }
  
  names(pcorr_faecal) = coh_l
  pdf(paste0(resultDirPath, "PearsonCorrPlot_Faecalibacterium.pdf"),onefile = T,height=10,width=10) 
  boxplot(
    pcorr_faecal$Nog_BCN0[c(3,4)],
    pcorr_faecal$Nog_STK[c(3,4)],
    pcorr_faecal$Escobar[c(3,4)],
    pcorr_faecal$ZellerFrance[c(3,4)],
    pcorr_faecal$ZellerGermany[c(3,4)],
    pcorr_faecal$Goodrich[c(3,4)],
    pcorr_faecal$Baxter[c(3,4)],
    pcorr_faecal$Ross[c(3,4)],
    pcorr_faecal$Gloor[c(3,4)],
    pcorr_faecal$Morgan[c(3,4)],
    pcorr_faecal$AGP[c(3,4)],
    names = c("NogBCN0\n(n=156)", "NogSTK\n(n=84)", "Escobar\n(n=30)",
              "Zeller\nFrance\n(n=129)", "Zeller\nGermany\n(n=48)", "Goodrich\n(n=835)", "Baxter\n(n=490)", "Ross\n(n=63)",
              "Gloor\n(n=803)", "Morgan\n(n=134)", "AGP\n(n=1,368)"),
    horizontal = T,
    ylim=c(-0.6, 0.6), 
    main="Faecalibacterium",
    xlab="Pearson R Value
  Age ~ OTU Count", 
    las=2, font =2, cex.axis=0.6
    
  )
  text(-0.18,11, "*", cex = 2)
  text(-0.24,9, "*", cex = 2)
  text(-0.25,7,"*", cex = 2)
  abline(v=0, lty=2)
  
  dev.off()
  
  # Comparing correlation trends across all cohorts
  
  ## Pulling enrichment (positive/negative correlation) for taxa common in all cohorts
  enrichment_list = list()
  as.data.frame(enrichment_list)
  
  
  for(i in 1:11){
    enrichment_commontaxa = matrix(nrow=length(commonGenusNames_allCohorts_genera), ncol=2) 
    common_stats_matrix= stats_fdrs_list[[i]][rownames(stats_fdrs_list[[i]]) %in% commonGenusNames_allCohorts_genera,]
    enrichment_commontaxa[,1] =rownames(common_stats_matrix)
    enrichment_commontaxa[,2] = common_stats_matrix[,4] 
    enrichment_list[[i]]= as.data.frame(enrichment_commontaxa)
  }
  
  all_enrichment_df = cbind(enrichment_list[[1]][,2], enrichment_list[[2]][,2], 
                            enrichment_list[[3]][,2], enrichment_list[[4]][,2],
                            enrichment_list[[5]][,2], enrichment_list[[6]][,2], 
                            enrichment_list[[7]][,2], enrichment_list[[8]][,2], 
                            enrichment_list[[9]][,2], enrichment_list[[10]][,2],
                            enrichment_list[[11]][,2])
  rownames(all_enrichment_df) = commonGenusNames_allCohorts_genera
  
  ## Counting how many taxa trend across all 11 cohorts
  total_en_sum = matrix(nrow=length(commonGenusNames_allCohorts_genera), ncol=2)
  for(i in 1:nrow(all_enrichment_df))
  {
    total_en_sum[i,2] = sum(all_enrichment_df[i,]== "negative") ==11
    total_en_sum[i,1] = sum(all_enrichment_df[i,]== "positive") ==11
  }
  rownames(total_en_sum) = commonGenusNames_allCohorts_genera
  sum(total_en_sum[,1]=="TRUE")
  sum(total_en_sum[,2]=="TRUE")
  
  ## Counting how many taxa trend across 10 out of 11 cohorts
  total_for10Cohorts = matrix(nrow=length(commonGenusNames_allCohorts_genera), ncol=2)
  for(i in 1:nrow(all_enrichment_df))
  {
    total_for10Cohorts[i,2] = sum(all_enrichment_df[i,]== "negative") ==10
    total_for10Cohorts[i,1] = sum(all_enrichment_df[i,]== "positive") ==10
  }
  rownames(total_for10Cohorts) = commonGenusNames_allCohorts_genera
  sum(total_for10Cohorts[,1]=="TRUE")
  sum(total_for10Cohorts[,2]=="TRUE")
  
}

## PCoA plots for pooled samples (all counts) and for each subgroup (Younger, Middle, Full Range)

#input: one file input (filepaths and cohort names as columns)
#output: pcoa plots for each subgroup

#Step 1: reads in counts table, labeled by cohort name
generatePcoaPlots = function (myFilePath, outputPath) {
  
  raw_q2021_filespaths = read.table(myFilePath, header = TRUE, sep = "\t")
  raw_q2021_files = list()
  
  for(i in 1:nrow(raw_q2021_filespaths)) {
    raw_q2021_files[[i]] = read.table(raw_q2021_filespaths[i,1], header = TRUE, sep = "\t")
  }
  
  names(raw_q2021_files) = raw_q2021_filespaths[,2]


#extracting only stool samples from Morgan dataset 
  raw_q2021_files[[1]] = subset(raw_q2021_files[[1]], raw_q2021_files[[1]]$stool==1)


#Step 2: generate PCoA plots 
  
## For pooled samples (all counts), first create df with merging two datasets, then looping the rest in
raw_merge_allCounts = full_join(raw_q2021_files[[1]], raw_q2021_files[[2]])

for(i in 3:11){ 
  file = raw_q2021_files[[i]]
  raw_merge_allCounts = full_join(raw_merge_allCounts, file)
}

raw_merge_allCounts[is.na(raw_merge_allCounts)] = 0

raw_merge_allCounts_genusTaxa = raw_merge_allCounts[,grep("g__", colnames(raw_merge_allCounts))]

lognorm_allCounts = lognorm_function(raw_merge_allCounts_genusTaxa)

zeroCountsRows = which(rowSums(raw_merge_allCounts_genusTaxa) ==0)
withoutZeroCountsRows = raw_merge_allCounts[-(zeroCountsRows),]

lognorm_allCounts$Study = withoutZeroCountsRows$Study
lognorm_allCounts$Age = withoutZeroCountsRows$Age
lognorm_allCounts$SampleID = withoutZeroCountsRows$SampleID

#permanova test
permanova_allcounts = adonis(lognorm_allCounts[,1:911] ~ factor(lognorm_allCounts$Study), permutations=999)
#could do...ncol(raw_merge_allCounts_genusTaxa)

#pcoa plot
lognorm_allCounts_q2021_pcoa = capscale(lognorm_allCounts[,1:911] ~ 1, distance = "bray")  

percentVariance_all <- eigenvals(lognorm_allCounts_q2021_pcoa)/sum(eigenvals(lognorm_allCounts_q2021_pcoa))

pdf(paste0(outputPath, "Pcoa_AllSamples.pdf"),onefile = T)
lognorm_allCounts_q2021_pcoa12 <- ordiplot(lognorm_allCounts_q2021_pcoa, choices = c(1,2), 
                                           main=c(paste("All Samples, Genus Level"),
                                                  paste("R2=", format(permanova_allcounts$aov.tab$R2[1],digits=3),
                                                        ", P-value=", permanova_allcounts$aov.tab$`Pr(>F)`[1])),
                                           type = "none", display="sites", cex.lab=1.5, xlab=paste("MDS1 (",format(percentVariance_all[1]*100,digits=3),"%)", sep=""), 
                                           ylab=paste("MDS 2 (",format(percentVariance_all[2]*100,digits=3),"%)", sep=""))

points(lognorm_allCounts_q2021_pcoa12, "sites",
       col=c(
         adjustcolor( "plum", alpha.f = 0.45), 
         adjustcolor( "lightblue", alpha.f = 0.45),
         adjustcolor( "limegreen", alpha.f = 0.45),
         adjustcolor( "orange", alpha.f = 0.45),
         adjustcolor( "purple", alpha.f = 0.45),
         adjustcolor( "black", alpha.f = 0.45),
         adjustcolor( "olivedrab", alpha.f = 0.45),
         adjustcolor( "red", alpha.f = 0.45),
         adjustcolor( "darkred", alpha.f = 0.45),
         adjustcolor( "cyan", alpha.f = 0.45),
         adjustcolor( "darkgreen", alpha.f = 0.45))[factor(lognorm_allCounts$Study)],pch=20, cex=2)

legend("topleft",
       levels(factor(lognorm_allCounts$Study)),
       col=c("plum", "lightblue", "limegreen", "orange", "purple",  
             "black", "olivedrab", "red", "darkred", "cyan","darkgreen"),
       pch=15,bty = "n")

dev.off() 

# boxplots for MDS1 and MDS2 with All Sample PCOA
generateMDSBoxplots_function = function(pcoaResult, mds, cohort){
  mds_output= pcoaResult$CA$u[,mds][lognorm_allCounts$Study==cohort]
  
  return(mds_output)
}

cohort_names = coh_l
mds1_per_cohort = sapply(cohort_names, generateMDSBoxplots_function,
                           pcoaResult= lognorm_allCounts_q2021_pcoa,
                           mds=1)
mds2_per_cohort = sapply(cohort_names, generateMDSBoxplots_function,
                         pcoaResult= lognorm_allCounts_q2021_pcoa,
                         mds=2)

myBoxplotColors = c("olivedrab", "red", "limegreen","cyan", "darkgreen",
                    "purple", "lightblue", "darkred", "orange", "darkgrey", "plum")
##boxplot : rearranged in order of asc sample size
pdf(paste0(outputPath, "AllSamples_MDS1_Boxplot.pdf"),onefile = T, width=14)
boxplot(
  mds1_per_cohort$Nog_BCN0,
  mds1_per_cohort$Nog_STK,
  mds1_per_cohort$Escobar,
  mds1_per_cohort$ZellerFrance,
  mds1_per_cohort$ZellerGermany,
  mds1_per_cohort$Goodrich,
  mds1_per_cohort$Baxter,
  mds1_per_cohort$Ross,
  mds1_per_cohort$Gloor,
  mds1_per_cohort$Morgan,
  mds1_per_cohort$AGP,
  cex.axis=0.8,
  names=c("NogBCN0", "NogSTK", "Escobar", "ZellerFrance", "ZellerGermany",
          "Goodrich", "Baxter","Ross","Gloor", "Morgan", "AGP"),
  main= "MDS1 per Cohort",
  xlab="Cohort",
  ylab="MDS1",
  border = myBoxplotColors)

mds1_allcohorts_list = list (
  "NogBCN0" = mds1_per_cohort$Nog_BCN0,
  "NogSTK" = mds1_per_cohort$Nog_STK,
  "Escobar" = mds1_per_cohort$Escobar,
  "Zeller France" = mds1_per_cohort$ZellerFrance,
  "Zeller Germany" = mds1_per_cohort$ZellerGermany,
  "Goodrich" = mds1_per_cohort$Goodrich,
  "Baxter" = mds1_per_cohort$Baxter,
  "Ross" = mds1_per_cohort$Ross,
  "Gloor" = mds1_per_cohort$Gloor,
  "Morgan" = mds1_per_cohort$Morgan,
  "AGP" = mds1_per_cohort$AGP)

stripchart(mds1_allcohorts_list, method = "jitter", 
           pch = 20, vertical= TRUE, add = TRUE, cex=1, col=myBoxplotColors)

dev.off()

pdf(paste0(outputPath, "AllSamples_MDS2_Boxplot.pdf"),onefile = T, width=14)
boxplot(
  mds2_per_cohort$Nog_BCN0,
  mds2_per_cohort$Nog_STK,
  mds2_per_cohort$Escobar,
  mds2_per_cohort$ZellerFrance,
  mds2_per_cohort$ZellerGermany,
  mds2_per_cohort$Goodrich,
  mds2_per_cohort$Baxter,
  mds2_per_cohort$Ross,
  mds2_per_cohort$Gloor,
  mds2_per_cohort$Morgan,
  mds2_per_cohort$AGP,
  cex.axis=0.8,
  names=c("NogBCN0", "NogSTK", "Escobar", "ZellerFrance", "ZellerGermany",
          "Goodrich", "Baxter","Ross","Gloor", "Morgan", "AGP"),
  main= "MDS2 per Cohort",
  xlab="Cohort",
  ylab="MDS2",
  border = myBoxplotColors)

mds2_allcohorts_list = list (
  "NogBCN0" = mds2_per_cohort$Nog_BCN0,
  "NogSTK" = mds2_per_cohort$Nog_STK,
  "Escobar" = mds2_per_cohort$Escobar,
  "Zeller France" = mds2_per_cohort$ZellerFrance,
  "Zeller Germany" = mds2_per_cohort$ZellerGermany,
  "Goodrich" = mds2_per_cohort$Goodrich,
  "Baxter" = mds2_per_cohort$Baxter,
  "Ross" = mds2_per_cohort$Ross,
  "Gloor" = mds2_per_cohort$Gloor,
  "Morgan" = mds2_per_cohort$Morgan,
  "AGP" = mds2_per_cohort$AGP)

stripchart(mds2_allcohorts_list, method = "jitter", 
           pch = 20, vertical= TRUE, add = TRUE, cex=1, col=myBoxplotColors)

dev.off()

  
## For subgroups(Younger, Middle, Full Range)

raw_merge_youngerRange = full_join(raw_q2021_files$NogBCN0, raw_q2021_files$NogSTK)
raw_merge_youngerRange_2 = full_join(raw_merge_youngerRange, raw_q2021_files$Escobar)

raw_merge_youngerRange_2[is.na(raw_merge_youngerRange_2)] = 0

raw_merge_youngerRange_2_genusTaxa = raw_merge_youngerRange_2[,grep("g__", colnames(raw_merge_youngerRange_2))]

lognorm_youngerRange = lognorm_function(raw_merge_youngerRange_2_genusTaxa)

lognorm_youngerRange$Study <- raw_merge_youngerRange_2$Study
lognorm_youngerRange$Age <- raw_merge_youngerRange_2$Age

#permanova test; could refer to ncol(raw_merge_youngerRange_2_genusTaxa) for adonis...
youngerRange_permanovaResults_cohort <- adonis(lognorm_youngerRange[,1:346] ~ factor(lognorm_youngerRange$Study), permutations = 999)

#pcoa plot
lognorm_youngerRange_pcoa <- capscale(lognorm_youngerRange[,1:346] ~ 1, distance = "bray")

percentVariance_youngerRange <- eigenvals(lognorm_youngerRange_pcoa)/sum(eigenvals(lognorm_youngerRange_pcoa))

pdf(paste0(outputPath, "Pcoa_YoungerRange.pdf"),onefile = T)
lognorm_youngerRange12 <- ordiplot(lognorm_youngerRange_pcoa, choices = c(1,2), 
                                   main=c(paste("Younger Range"),
                                          paste("R2=", format(youngerRange_permanovaResults_cohort$aov.tab$R2[1],digits=3),
                                                ", P-value=", youngerRange_permanovaResults_cohort$aov.tab$`Pr(>F)`[1])),
                                   type = "none", display="sites", cex.lab=1.5, xlab=paste("MDS1 (",format(percentVariance_youngerRange[1]*100,digits=3),"%)", sep=""), 
                                   ylab=paste("MDS 2 (",format(percentVariance_youngerRange[2]*100,digits=3),"%)", sep=""))
points(lognorm_youngerRange12, "sites",
       col=c(adjustcolor( "olivedrab", alpha.f = 0.35), 
             adjustcolor( "red", alpha.f = 0.35),
             adjustcolor( "limegreen", alpha.f = 0.35))[factor(lognorm_youngerRange$Study)],pch=20, cex=2)
legend("topleft",
       levels(factor(lognorm_youngerRange$Study)),
       col=c("olivedrab", "red", "limegreen"),
       pch=15,bty = "n")
dev.off()


## middle range ______________________________

raw_merge_middleRange <- full_join(raw_q2021_files$Baxter, raw_q2021_files$Ross)
raw_merge_middleRange_2 <- full_join(raw_merge_middleRange, raw_q2021_files$Goodrich)
raw_merge_middleRange_3 <- full_join(raw_merge_middleRange_2, raw_q2021_files$ZellerFrance)
raw_merge_middleRange_4 <- full_join(raw_merge_middleRange_3, raw_q2021_files$ZellerGermany)

raw_merge_middleRange_4[is.na(raw_merge_middleRange_4)] = 0

raw_middleRange_taxaGenus = raw_merge_middleRange_4[,grep("g__",colnames(raw_merge_middleRange_4))]
raw_middleRange_taxaGenus = raw_middleRange_taxaGenus[,colSums(raw_middleRange_taxaGenus) !=0]

lognorm_middleRange = lognorm_function(raw_middleRange_taxaGenus)

lognorm_middleRange$Study <- raw_merge_middleRange_4$Study
lognorm_middleRange$Age <- raw_merge_middleRange_4$Age

#permanova test 
middleRange_permanovaResults_cohort <- adonis(lognorm_middleRange[,1:602] ~ factor(lognorm_middleRange$Study), permutations = 999)

#pcoa plot
lognorm_middleRange_pcoa <- capscale(lognorm_middleRange[,1:602] ~ 1, distance = "bray")

percentVariance_middleRange <- eigenvals(lognorm_middleRange_pcoa)/sum(eigenvals(lognorm_middleRange_pcoa))

pdf(paste0(outputPath, "Pcoa_MiddleRange.pdf"),onefile = T)
lognorm_middleRange12 <- ordiplot(lognorm_middleRange_pcoa, choices = c(1,2), 
                                  main=c(paste("Middle Range"),
                                         paste("R2=", format(middleRange_permanovaResults_cohort$aov.tab$R2[1],digits=3),
                                               ", P-value=", middleRange_permanovaResults_cohort$aov.tab$`Pr(>F)`[1])),
                                  type = "none", display="sites", cex.lab=1.5, xlab=paste("MDS1 (",format(percentVariance_middleRange[1]*100,digits=3),"%)", sep=""), 
                                  ylab=paste("MDS 2 (",format(percentVariance_middleRange[2]*100,digits=3),"%)", sep=""))
points(lognorm_middleRange12, "sites",
       col=c(adjustcolor( "lightblue", alpha.f = 0.45), 
             adjustcolor( "darkred", alpha.f = 0.45),
             adjustcolor( "purple", alpha.f = 0.45),
             adjustcolor( "cyan", alpha.f = 0.45),
             adjustcolor( "darkgreen", alpha.f = 0.45))[factor(lognorm_middleRange$Study)],pch=20, cex=2)


legend("topleft",
       levels(factor(lognorm_middleRange$Study)),
       col=c("lightblue", "darkred", "purple", "cyan", "darkgreen"),
       pch=15,bty = "n")

dev.off()
### full range _______________________
raw_merge_fullRange_a <- full_join(raw_q2021_files$Gloor,raw_q2021_files$Morgan)
raw_merge_fullRange = full_join(raw_merge_fullRange_a,raw_q2021_files$AGP)

raw_merge_fullRange[is.na(raw_merge_fullRange)] = 0

raw_fullRange_taxaGenus = raw_merge_fullRange[,grep("g__", colnames(raw_merge_fullRange))]
raw_fullRange_taxaGenus = raw_fullRange_taxaGenus[,colSums(raw_fullRange_taxaGenus)!= 0]

lognorm_fullRange = lognorm_function(raw_fullRange_taxaGenus)

agp.rowsWithZeroCounts.index = which(rowSums(raw_fullRange_taxaGenus) ==0)
rawCounts_agp_zeroRowsRemoved = raw_merge_fullRange[-(agp.rowsWithZeroCounts.index),]

lognorm_fullRange$Study <- rawCounts_agp_zeroRowsRemoved$Study
lognorm_fullRange$Age <- rawCounts_agp_zeroRowsRemoved$Age

#permanova test
fullRange_permanovaResults <- adonis(lognorm_fullRange[,1:ncol(raw_fullRange_taxaGenus)] ~ factor(lognorm_fullRange$Study) + lognorm_fullRange$Age, permutations = 999)

#pcoa plot
lognorm_fullRange_pcoa <- capscale(lognorm_fullRange[,1:ncol(raw_fullRange_taxaGenus)] ~ 1, distance = "bray")

percentVariance_fullRange <- eigenvals(lognorm_fullRange_pcoa)/sum(eigenvals(lognorm_fullRange_pcoa))

pdf(paste0(outputPath, "Pcoa_FullRange.pdf"),onefile = T)
lognorm_fullRange12 <- ordiplot(lognorm_fullRange_pcoa, choices = c(1,2), 
                                main=c(paste("Full Range"),
                                       paste("R2=", format(fullRange_permanovaResults$aov.tab$R2[1],digits=3),
                                             ", P-value=", fullRange_permanovaResults$aov.tab$`Pr(>F)`[1])),
                                type = "none", display="sites", cex.lab=1.5, xlab=paste("MDS1 (",format(percentVariance_fullRange[1]*100,digits=3),"%)", sep=""), 
                                ylab=paste("MDS 2 (",format(percentVariance_fullRange[2]*100,digits=3),"%)", sep=""))
points(lognorm_fullRange12, "sites",
       col=c(adjustcolor( "plum", alpha.f = 0.35), 
             adjustcolor( "orange", alpha.f = 0.35),
             adjustcolor( "black", alpha.f = 0.35))[factor(lognorm_fullRange$Study)],pch=20, cex=2)
legend("topleft",
       levels(factor(lognorm_fullRange$Study)),
       col=c("plum", "orange", "black"),
       pch=15,bty = "n")

dev.off()
}


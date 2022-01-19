## pcoa for subgroups
#starting w raw counts and log10 normalizing per subgroup
#counts from qiime2 v2021.2

rawCounts_baxter = read.table("/Users/anhil/Desktop/qiime2_2021version_countstables/rawCounts_baxter.txt", header = TRUE)
rawCounts_escobar = read.table("/Users/anhil/Desktop/qiime2_2021version_countstables/rawCounts_escobar.txt", header = TRUE)
rawCounts_gloor = read.table("/Users/anhil/Desktop/qiime2_2021version_countstables/rawCounts_gloor.txt", header = TRUE)
rawCounts_goodrich = read.table("/Users/anhil/Desktop/qiime2_2021version_countstables/rawCounts_goodrich.txt", header = TRUE)
rawCounts_morgan = read.table("/Users/anhil/Desktop/qiime2_2021version_countstables/rawCounts_morgan.txt", header = TRUE)
rawCounts_ross = read.table("/Users/anhil/Desktop/qiime2_2021version_countstables/rawCounts_ross.txt", header = TRUE)
rawCounts_nogbcn0 = read.table("/Users/anhil/Desktop/qiime2_2021version_countstables/rawCounts_nogbcn0.txt", header = TRUE)
rawCounts_nogstk = read.table("/Users/anhil/Desktop/qiime2_2021version_countstables/rawCounts_nogstk.txt", header = TRUE)
rawCounts_zellerfrance = read.table("/Users/anhil/Desktop/qiime2_2021version_countstables/rawCounts_zellerfrance.txt", header = TRUE)
rawCounts_zellergermany = read.table("/Users/anhil/Desktop/qiime2_2021version_countstables/rawCounts_zellergermany.txt", header = TRUE)

raw_merge_youngerRange <- full_join(rawCounts_nogbcn0, rawCounts_nogstk)
raw_merge_youngerRange_2 <- full_join(raw_merge_youngerRange, rawCounts_escobar)

raw_merge_youngerRange_2[is.na(raw_merge_youngerRange_2)] = 0

raw_merge_youngerRange_2_genusTaxa = raw_merge_youngerRange_2[,grep("g__", colnames(raw_merge_youngerRange_2))]

lognorm_youngerRange = lognorm_function(raw_merge_youngerRange_2_genusTaxa)

lognorm_youngerRange$Study <- raw_merge_youngerRange_2$Study
lognorm_youngerRange$Age <- raw_merge_youngerRange_2$Age

lognorm_youngerRange_pcoa <- capscale(lognorm_youngerRange[,1:346] ~ 1, distance = "bray")

percentVariance_youngerRange <- eigenvals(lognorm_youngerRange_pcoa)/sum(eigenvals(lognorm_youngerRange_pcoa))

lognorm_youngerRange12 <- ordiplot(lognorm_youngerRange_pcoa, choices = c(1,2), 
                                   main="Genus Level, Younger Range,
                                  R2=0.171, P-value=0.001",
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

youngerRange_permanovaResults <- adonis(lognorm_youngerRange[,1:346] ~ factor(lognorm_youngerRange$Study) + lognorm_youngerRange$Age, permutations = 999)
youngerRange_permanovaResults_cohort <- adonis(lognorm_youngerRange[,1:346] ~ factor(lognorm_youngerRange$Study), permutations = 999)


## middle range ______________________________

raw_merge_middleRange <- full_join(rawCounts_baxter, rawCounts_ross)
raw_merge_middleRange_2 <- full_join(raw_merge_middleRange, rawCounts_goodrich)
raw_merge_middleRange_3 <- full_join(raw_merge_middleRange_2, rawCounts_zellerfrance)
raw_merge_middleRange_4 <- full_join(raw_merge_middleRange_3, rawCounts_zellergermany)

raw_merge_middleRange_4[is.na(raw_merge_middleRange_4)] = 0

raw_middleRange_taxaGenus = raw_merge_middleRange_4[,grep("g__",colnames(raw_merge_middleRange_4))]
raw_middleRange_taxaGenus = raw_middleRange_taxaGenus[,colSums(raw_middleRange_taxaGenus) !=0]

lognorm_middleRange = lognorm_function(raw_middleRange_taxaGenus)


lognorm_middleRange$Study <- raw_merge_middleRange_4$Study
lognorm_middleRange$Age <- raw_merge_middleRange_4$Age

lognorm_middleRange_pcoa <- capscale(lognorm_middleRange[,1:602] ~ 1, distance = "bray")

percentVariance_middleRange <- eigenvals(lognorm_middleRange_pcoa)/sum(eigenvals(lognorm_middleRange_pcoa))

lognorm_middleRange12 <- ordiplot(lognorm_middleRange_pcoa, choices = c(1,2), 
                                  main="Genus Level, Middle Range,
                                  R2=0.168, P-value=0.001",
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

middleRange_permanovaResults <- adonis(lognorm_middleRange[,1:602] ~ factor(lognorm_middleRange$Study) + lognorm_middleRange$Age, permutations = 999)
middleRange_permanovaResults_cohort <- adonis(lognorm_middleRange[,1:602] ~ factor(lognorm_middleRange$Study), permutations = 999)

### full range _______________________
current_agp = raw_q2021_files[[11]]
raw_merge_fullRange_a <- full_join(rawCounts_gloor,rawCounts_morgan)
raw_merge_fullRange = full_join(raw_merge_fullRange_a,current_agp)

raw_merge_fullRange[is.na(raw_merge_fullRange)] = 0

raw_fullRange_taxaGenus = raw_merge_fullRange[,grep("g__", colnames(raw_merge_fullRange))]
raw_fullRange_taxaGenus = raw_fullRange_taxaGenus[,colSums(raw_fullRange_taxaGenus)!= 0]

lognorm_fullRange = lognorm_function(raw_fullRange_taxaGenus)

agp.rowsWithZeroCounts.index = which(rowSums(raw_fullRange_taxaGenus) ==0)
rawCounts_agp_zeroRowsRemoved = raw_merge_fullRange[-(agp.rowsWithZeroCounts.index),]

###in the middle of testing if agp_ndups makes that big of a diff without the dup samples...
lognorm_fullRange$Study <- rawCounts_agp_zeroRowsRemoved$Study
lognorm_fullRange$Age <- rawCounts_agp_zeroRowsRemoved$Age

lognorm_fullRange_pcoa <- capscale(lognorm_fullRange[,1:712] ~ 1, distance = "bray")

percentVariance_fullRange <- eigenvals(lognorm_fullRange_pcoa)/sum(eigenvals(lognorm_fullRange_pcoa))

lognorm_fullRange12 <- ordiplot(lognorm_fullRange_pcoa, choices = c(1,2), 
                                main="Genus Level, Full Range,
                                  R2=0.178, P-value=0.001",
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

fullRange_permanovaResults <- adonis(lognorm_fullRange[,1:712] ~ factor(lognorm_fullRange$Study) + lognorm_fullRange$Age, permutations = 999)



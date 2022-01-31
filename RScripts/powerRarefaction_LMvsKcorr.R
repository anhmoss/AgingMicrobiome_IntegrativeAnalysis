generatePowerPlots_CompareModels = function(myFilePath,resultDirPath) {

raw_q2021_filespaths = read.table(myFilePath, header = TRUE, sep = "\t")
raw_q2021_files = list()

for(i in 1:nrow(raw_q2021_filespaths)) {
  raw_q2021_files[[i]] = read.table(raw_q2021_filespaths[i,1], header = TRUE, sep = "\t")
}

names(raw_q2021_files) = raw_q2021_filespaths[,2]

coh_l = vector()
lognormfiles_q2021 = list()

for(i in 1:11){
  raw_file = raw_q2021_files[[i]]
  taxa_table = raw_file[,grep("g__",colnames(raw_file))]
  age_meta = raw_file$Age[rowSums(taxa_table)!=0]
  
  lognorm_file = lognorm_function(taxa_table)
  names_short=sapply(strsplit(colnames(lognorm_file),"g__"),"[[",2)
  colnames(lognorm_file)=names_short
  lognorm_file=lognorm_file[,colSums(lognorm_file)!=0]

  lognorm_file$Age = age_meta
  lognormfiles_q2021[[i]] = lognorm_file
  
  coh_l[i]=raw_file[1,1]
}

#agp
lognorm_output_agp = lognormfiles_q2021[[11]]
lognorm_output_agp_taxa = lognormfiles_q2021[[11]][,-ncol(lognormfiles_q2021[[11]])]
lognorm_output_agp_age = lognormfiles_q2021[[11]][,ncol(lognormfiles_q2021[[11]])]

agp_intervals = c(20,40,60,80, 100, 120, 140, 160, 180, 
                  200, 300,400, 500, 600, 700, 800, 900, 1000,
                  1100, 1200, 1368)

lognorm_agp_powerMatrix_LM = matrix(nrow=length(agp_intervals), ncol = 3)
lognorm_agp_powerMatrix_KCorr = matrix(nrow=length(agp_intervals), ncol = 3)

for(i in 1:length(agp_intervals)){
  lognorm_agp_powerMatrix_LM[i,] = unlist(repeatSubsampling_LM_function(maxiteration = 50,
                                                                     subsamplesize = agp_intervals[i],
                                                                     lognormFile = lognorm_output_agp_taxa,
                                                                     metadata = lognorm_output_agp_age))

  lognorm_agp_powerMatrix_KCorr[i,] = unlist(repeatSubsampling_corr_function(maxiteration = 50,
                                                                             subsamplesize = agp_intervals[i],
                                                                             lognormFile = lognorm_output_agp_taxa,
                                                                             metadata = lognorm_output_agp_age,
                                                                             corrMethod="kendall"))                                                                   
}
colnames(lognorm_agp_powerMatrix_LM) = c("Subsamplesize", "Mean", "SD")
lognorm_agp_powerMatrix_LM = as.data.frame(lognorm_agp_powerMatrix_LM)

colnames(lognorm_agp_powerMatrix_KCorr) = c("Subsamplesize", "Mean", "SD")
lognorm_agp_powerMatrix_KCorr = as.data.frame(lognorm_agp_powerMatrix_KCorr)


#ross 
lognorm_output_ross = lognormfiles_q2021[[10]]
lognorm_output_ross_taxa = lognormfiles_q2021[[10]][,-ncol(lognormfiles_q2021[[10]])]
lognorm_output_ross_age = lognormfiles_q2021[[10]][,ncol(lognormfiles_q2021[[10]])]

ross_intervals = c(10,20,30,40,50,63)

lognorm_ross_powerMatrix_LM = matrix(nrow=length(ross_intervals), ncol = 3)
lognorm_ross_powerMatrix_KCorr = matrix(nrow=length(ross_intervals), ncol = 3)

for(i in 1:length(ross_intervals)){
  lognorm_ross_powerMatrix_LM[i,] = unlist(repeatSubsampling_LM_function(maxiteration = 50,
                                                                         subsamplesize = ross_intervals[i],
                                                                         lognormFile = lognorm_output_ross_taxa,
                                                                         metadata = lognorm_output_ross_age))
  lognorm_ross_powerMatrix_KCorr[i,] = unlist(repeatSubsampling_corr_function(maxiteration = 50,
                                                                              subsamplesize = ross_intervals[i],
                                                                              lognormFile = lognorm_output_ross_taxa,
                                                                              metadata = lognorm_output_ross_age,
                                                                              corrMethod="kendall"))
}

colnames(lognorm_ross_powerMatrix_LM) = c("Subsamplesize", "Mean", "SD")
lognorm_ross_powerMatrix_LM = as.data.frame(lognorm_ross_powerMatrix_LM)

colnames(lognorm_ross_powerMatrix_KCorr) = c("Subsamplesize", "Mean", "SD")
lognorm_ross_powerMatrix_KCorr = as.data.frame(lognorm_ross_powerMatrix_KCorr)

#nogstk
lognorm_output_nogstk = lognormfiles_q2021[[9]]
lognorm_output_nogstk_taxa = lognormfiles_q2021[[9]][,-ncol(lognormfiles_q2021[[9]])]
lognorm_output_nogstk_age = lognormfiles_q2021[[9]][,ncol(lognormfiles_q2021[[9]])]

nogstk_intervals = c(10,20,30,40,50,60,70,84)

lognorm_nogstk_powerMatrix_LM = matrix(nrow=length(nogstk_intervals), ncol = 3)
lognorm_nogstk_powerMatrix_KCorr = matrix(nrow=length(nogstk_intervals), ncol = 3)

for(i in 1:length(nogstk_intervals)){
  lognorm_nogstk_powerMatrix_LM[i,] = unlist(repeatSubsampling_LM_function(maxiteration = 50,
                                                                           subsamplesize = nogstk_intervals[i],
                                                                           lognormFile = lognorm_output_nogstk_taxa,
                                                                           metadata = lognorm_output_nogstk_age))
  lognorm_nogstk_powerMatrix_KCorr[i,] = unlist(repeatSubsampling_corr_function(maxiteration = 50,
                                                                                subsamplesize = nogstk_intervals[i],
                                                                                lognormFile = lognorm_output_nogstk_taxa,
                                                                                metadata = lognorm_output_nogstk_age,
                                                                                corrMethod="kendall"))
}

colnames(lognorm_nogstk_powerMatrix_LM) = c("Subsamplesize", "Mean", "SD")
lognorm_nogstk_powerMatrix_LM = as.data.frame(lognorm_nogstk_powerMatrix_LM)

colnames(lognorm_nogstk_powerMatrix_KCorr) = c("Subsamplesize", "Mean", "SD")
lognorm_nogstk_powerMatrix_KCorr = as.data.frame(lognorm_nogstk_powerMatrix_KCorr)


#nogbcn0
lognorm_output_nogbcn0 = lognormfiles_q2021[[8]]
lognorm_output_nogbcn0_taxa = lognormfiles_q2021[[8]][,-ncol(lognormfiles_q2021[[8]])]
lognorm_output_nogbcn0_age = lognormfiles_q2021[[8]][,ncol(lognormfiles_q2021[[8]])]

nogbcn0_intervals = c(20,40,60,80,100,120,140,156)

lognorm_nogbcn0_powerMatrix_LM = matrix(nrow=length(nogbcn0_intervals), ncol = 3)
lognorm_nogbcn0_powerMatrix_KCorr = matrix(nrow=length(nogbcn0_intervals), ncol = 3)

for(i in 1:length(nogbcn0_intervals)){
  lognorm_nogbcn0_powerMatrix_LM[i,] = unlist(repeatSubsampling_LM_function(maxiteration = 50,
                                                                            subsamplesize = nogbcn0_intervals[i],
                                                                            lognormFile = lognorm_output_nogbcn0_taxa,
                                                                            metadata = lognorm_output_nogbcn0_age))
  lognorm_nogbcn0_powerMatrix_KCorr[i,] = unlist(repeatSubsampling_corr_function(maxiteration = 50,
                                                                                 subsamplesize = nogbcn0_intervals[i],
                                                                                 lognormFile = lognorm_output_nogbcn0_taxa,
                                                                                 metadata = lognorm_output_nogbcn0_age,
                                                                                 corrMethod="kendall"))
}

colnames(lognorm_nogbcn0_powerMatrix_LM) = c("Subsamplesize", "Mean", "SD")
lognorm_nogbcn0_powerMatrix_LM = as.data.frame(lognorm_nogbcn0_powerMatrix_LM)

colnames(lognorm_nogbcn0_powerMatrix_KCorr) = c("Subsamplesize", "Mean", "SD")
lognorm_nogbcn0_powerMatrix_KCorr = as.data.frame(lognorm_nogbcn0_powerMatrix_KCorr)


#goodrich
lognorm_output_goodrich = lognormfiles_q2021[[7]]
lognorm_output_goodrich_taxa = lognormfiles_q2021[[7]][,-ncol(lognormfiles_q2021[[7]])]
lognorm_output_goodrich_age = lognormfiles_q2021[[7]][,ncol(lognormfiles_q2021[[7]])]

goodrich_intervals = c(20,40,60,80,100,120,140,160,180,
                       200,300,400,500,600,700,800,835)

lognorm_goodrich_powerMatrix_LM = matrix(nrow=length(goodrich_intervals), ncol = 3)
lognorm_goodrich_powerMatrix_KCorr = matrix(nrow=length(goodrich_intervals), ncol = 3)

for(i in 1:length(goodrich_intervals)){
  lognorm_goodrich_powerMatrix_LM[i,] = unlist(repeatSubsampling_LM_function(maxiteration = 50,
                                                                             subsamplesize = goodrich_intervals[i],
                                                                             lognormFile = lognorm_output_goodrich_taxa,
                                                                             metadata = lognorm_output_goodrich_age))

  lognorm_goodrich_powerMatrix_KCorr[i,] = unlist(repeatSubsampling_corr_function(maxiteration = 50,
                                                                                  subsamplesize = goodrich_intervals[i],
                                                                                  lognormFile = lognorm_output_goodrich_taxa,
                                                                                  metadata = lognorm_output_goodrich_age,
                                                                                  corrMethod="kendall"))
}

colnames(lognorm_goodrich_powerMatrix_LM) = c("Subsamplesize", "Mean", "SD")
lognorm_goodrich_powerMatrix_LM = as.data.frame(lognorm_goodrich_powerMatrix_LM)

colnames(lognorm_goodrich_powerMatrix_KCorr) = c("Subsamplesize", "Mean", "SD")
lognorm_goodrich_powerMatrix_KCorr = as.data.frame(lognorm_goodrich_powerMatrix_KCorr)



#gloor
lognorm_output_gloor = lognormfiles_q2021[[6]]
lognorm_output_gloor_taxa = lognormfiles_q2021[[6]][,-ncol(lognormfiles_q2021[[6]])]
lognorm_output_gloor_age = lognormfiles_q2021[[6]][,ncol(lognormfiles_q2021[[6]])]

gloor_intervals = c(20,40,60,80,100,120,140,160,180,
                    200,300,400,500,600,700,803)

lognorm_gloor_powerMatrix_LM = matrix(nrow=length(gloor_intervals), ncol = 3)
lognorm_gloor_powerMatrix_KCorr = matrix(nrow=length(gloor_intervals), ncol = 3)

for(i in 1:length(gloor_intervals)){
  lognorm_gloor_powerMatrix_LM[i,] = unlist(repeatSubsampling_LM_function(maxiteration = 50,
                                                                          subsamplesize = gloor_intervals[i],
                                                                          lognormFile = lognorm_output_gloor_taxa,
                                                                          metadata = lognorm_output_gloor_age))
  lognorm_gloor_powerMatrix_KCorr[i,] = unlist(repeatSubsampling_corr_function(maxiteration = 50,
                                                                               subsamplesize = gloor_intervals[i],
                                                                               lognormFile = lognorm_output_gloor_taxa,
                                                                               metadata = lognorm_output_gloor_age,
                                                                               corrMethod="kendall"))
}

colnames(lognorm_gloor_powerMatrix_LM) = c("Subsamplesize", "Mean", "SD")
lognorm_gloor_powerMatrix_LM = as.data.frame(lognorm_gloor_powerMatrix_LM)

colnames(lognorm_gloor_powerMatrix_KCorr) = c("Subsamplesize", "Mean", "SD")
lognorm_gloor_powerMatrix_KCorr = as.data.frame(lognorm_gloor_powerMatrix_KCorr)


#escobar
lognorm_output_escobar = lognormfiles_q2021[[5]]
lognorm_output_escobar_taxa = lognormfiles_q2021[[5]][,-ncol(lognormfiles_q2021[[5]])]
lognorm_output_escobar_age = lognormfiles_q2021[[5]][,ncol(lognormfiles_q2021[[5]])]

escobar_intervals = c(10,20,30)

lognorm_escobar_powerMatrix_LM = matrix(nrow=length(escobar_intervals), ncol = 3)
lognorm_escobar_powerMatrix_KCorr = matrix(nrow=length(escobar_intervals), ncol = 3)

for(i in 1:length(escobar_intervals)){
  lognorm_escobar_powerMatrix_LM[i,] = unlist(repeatSubsampling_LM_function(maxiteration = 50,
                                                                            subsamplesize = escobar_intervals[i],
                                                                            lognormFile = lognorm_output_escobar_taxa,
                                                                            metadata = lognorm_output_escobar_age))
  lognorm_escobar_powerMatrix_KCorr[i,] = unlist(repeatSubsampling_corr_function(maxiteration = 50,
                                                                                 subsamplesize = escobar_intervals[i],
                                                                                 lognormFile = lognorm_output_escobar_taxa,
                                                                                 metadata = lognorm_output_escobar_age,
                                                                                 corrMethod="kendall"))
}

colnames(lognorm_escobar_powerMatrix_LM) = c("Subsamplesize", "Mean", "SD")
lognorm_escobar_powerMatrix_LM = as.data.frame(lognorm_escobar_powerMatrix_LM)
colnames(lognorm_escobar_powerMatrix_KCorr) = c("Subsamplesize", "Mean", "SD")
lognorm_escobar_powerMatrix_KCorr = as.data.frame(lognorm_escobar_powerMatrix_KCorr)


#baxter
lognorm_output_baxter = lognormfiles_q2021[[4]]
lognorm_output_baxter_taxa = lognormfiles_q2021[[4]][,-ncol(lognormfiles_q2021[[4]])]
lognorm_output_baxter_age = lognormfiles_q2021[[4]][,ncol(lognormfiles_q2021[[4]])]

baxter_intervals = c(20,40,60,80,100,120,140,160,180,
                     200,220,240,260,280,300, 320, 340,
                     360,380,400, 420,440,460,490)

lognorm_baxter_powerMatrix_LM = matrix(nrow=length(baxter_intervals), ncol = 3)
lognorm_baxter_powerMatrix_KCorr = matrix(nrow=length(baxter_intervals), ncol = 3)

for(i in 1:length(baxter_intervals)){
  lognorm_baxter_powerMatrix_LM[i,] = unlist(repeatSubsampling_LM_function(maxiteration = 50,
                                                                           subsamplesize = baxter_intervals[i],
                                                                           lognormFile = lognorm_output_baxter_taxa,
                                                                           metadata = lognorm_output_baxter_age))
  lognorm_baxter_powerMatrix_KCorr[i,] = unlist(repeatSubsampling_corr_function(maxiteration = 50,
                                                                                subsamplesize = baxter_intervals[i],
                                                                                lognormFile = lognorm_output_baxter_taxa,
                                                                                metadata = lognorm_output_baxter_age,
                                                                                corrMethod="kendall"))
}

colnames(lognorm_baxter_powerMatrix_LM) = c("Subsamplesize", "Mean", "SD")
lognorm_baxter_powerMatrix_LM = as.data.frame(lognorm_baxter_powerMatrix_LM)
colnames(lognorm_baxter_powerMatrix_KCorr) = c("Subsamplesize", "Mean", "SD")
lognorm_baxter_powerMatrix_KCorr = as.data.frame(lognorm_baxter_powerMatrix_KCorr)


#zellergermany
lognorm_output_zellergermany = lognormfiles_q2021[[3]]
lognorm_output_zellergermany_taxa = lognormfiles_q2021[[3]][,-ncol(lognormfiles_q2021[[3]])]
lognorm_output_zellergermany_age = lognormfiles_q2021[[3]][,ncol(lognormfiles_q2021[[3]])]

zellergermany_intervals = c(10,20,30,40,48)

lognorm_zellergermany_powerMatrix_LM = matrix(nrow=length(zellergermany_intervals), ncol = 3)
lognorm_zellergermany_powerMatrix_KCorr = matrix(nrow=length(zellergermany_intervals), ncol = 3)

for(i in 1:length(zellergermany_intervals)){
  lognorm_zellergermany_powerMatrix_LM[i,] = unlist(repeatSubsampling_LM_function(maxiteration = 50,
                                                                                  subsamplesize = zellergermany_intervals[i],
                                                                                  lognormFile = lognorm_output_zellergermany_taxa,
                                                                                  metadata = lognorm_output_zellergermany_age))
lognorm_zellergermany_powerMatrix_KCorr[i,] = unlist(repeatSubsampling_corr_function(maxiteration = 50,
                                                                                       subsamplesize = zellergermany_intervals[i],
                                                                                       lognormFile = lognorm_output_zellergermany_taxa,
                                                                                       metadata = lognorm_output_zellergermany_age,
                                                                                       corrMethod="kendall"))
}

colnames(lognorm_zellergermany_powerMatrix_LM) = c("Subsamplesize", "Mean", "SD")
lognorm_zellergermany_powerMatrix_LM = as.data.frame(lognorm_zellergermany_powerMatrix_LM)
colnames(lognorm_zellergermany_powerMatrix_KCorr) = c("Subsamplesize", "Mean", "SD")
lognorm_zellergermany_powerMatrix_KCorr = as.data.frame(lognorm_zellergermany_powerMatrix_KCorr)


#zellerfrance
lognorm_output_zellerfrance = lognormfiles_q2021[[2]]
lognorm_output_zellerfrance_taxa = lognormfiles_q2021[[2]][,-ncol(lognormfiles_q2021[[2]])]
lognorm_output_zellerfrance_age = lognormfiles_q2021[[2]][,ncol(lognormfiles_q2021[[2]])]

zellerfrance_intervals = c(10,20,30,40,50,60,70,80,90,100,
                           120,129)

lognorm_zellerfrance_powerMatrix_LM = matrix(nrow=length(zellerfrance_intervals), ncol = 3)
lognorm_zellerfrance_powerMatrix_KCorr = matrix(nrow=length(zellerfrance_intervals), ncol = 3)

for(i in 1:length(zellerfrance_intervals)){
  lognorm_zellerfrance_powerMatrix_LM[i,] = unlist(repeatSubsampling_LM_function(maxiteration = 50,
                                                                                 subsamplesize = zellerfrance_intervals[i],
                                                                                 lognormFile = lognorm_output_zellerfrance_taxa,
                                                                                 metadata = lognorm_output_zellerfrance_age))
  lognorm_zellerfrance_powerMatrix_KCorr[i,] = unlist(repeatSubsampling_corr_function(maxiteration = 50,
                                                                                      subsamplesize = zellerfrance_intervals[i],
                                                                                      lognormFile = lognorm_output_zellerfrance_taxa,
                                                                                      metadata = lognorm_output_zellerfrance_age,
                                                                                      corrMethod="kendall"))
}

colnames(lognorm_zellerfrance_powerMatrix_LM) = c("Subsamplesize", "Mean", "SD")
lognorm_zellerfrance_powerMatrix_LM = as.data.frame(lognorm_zellerfrance_powerMatrix_LM)
colnames(lognorm_zellerfrance_powerMatrix_KCorr) = c("Subsamplesize", "Mean", "SD")
lognorm_zellerfrance_powerMatrix_KCorr = as.data.frame(lognorm_zellerfrance_powerMatrix_KCorr)


#morgan
lognorm_output_morgan = lognormfiles_q2021[[1]]
lognorm_output_morgan_taxa = lognormfiles_q2021[[1]][,-ncol(lognormfiles_q2021[[1]])]
lognorm_output_morgan_age = lognormfiles_q2021[[1]][,ncol(lognormfiles_q2021[[1]])]

morgan_intervals = c(10,20,30,40,50,60,70,80,90,100,
                     120,140,160,180,200,228)

lognorm_morgan_powerMatrix_LM = matrix(nrow=length(morgan_intervals), ncol = 3)
lognorm_morgan_powerMatrix_KCorr = matrix(nrow=length(morgan_intervals), ncol = 3)

for(i in 1:length(morgan_intervals)){
  lognorm_morgan_powerMatrix_LM[i,] = unlist(repeatSubsampling_LM_function(maxiteration = 50,
                                                                           subsamplesize = morgan_intervals[i],
                                                                           lognormFile = lognorm_output_morgan_taxa,
                                                                           metadata = lognorm_output_morgan_age))
  lognorm_morgan_powerMatrix_KCorr[i,] = unlist(repeatSubsampling_corr_function(maxiteration = 50,
                                                                                subsamplesize = morgan_intervals[i],
                                                                                lognormFile = lognorm_output_morgan_taxa,
                                                                                metadata = lognorm_output_morgan_age,
                                                                                corrMethod="kendall"))
}

colnames(lognorm_morgan_powerMatrix_LM) = c("Subsamplesize", "Mean", "SD")
lognorm_morgan_powerMatrix_LM = as.data.frame(lognorm_morgan_powerMatrix_LM)
colnames(lognorm_morgan_powerMatrix_KCorr) = c("Subsamplesize", "Mean", "SD")
lognorm_morgan_powerMatrix_KCorr = as.data.frame(lognorm_morgan_powerMatrix_KCorr)




# Combine LM and KCorr power results into one plot per cohort that had significant taxa_table

#agp
compare_lm_agp = lognorm_agp_powerMatrix_LM
compare_kcorr_agp = lognorm_agp_powerMatrix_KCorr

compare_lm_agp$test = rep("Parametric", nrow(compare_lm_agp))
compare_kcorr_agp$test = rep("Non-Parametric", nrow(compare_kcorr_agp))

compare_merge_agp = rbind(compare_lm_agp, compare_kcorr_agp)

compareStatModels_plot_agp = ggplot(compare_merge_agp, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = compare_merge_agp, aes(x= Subsamplesize, y = Mean, color=test)) +
  geom_errorbar(aes(ymin = Mean - SD, 
                    ymax = Mean + SD, color=test)) +
  ggtitle("AGP(Non-Parametric vs Parametric) n=1,368") +
  scale_x_continuous(trans = 'log10') + #scale_y_continuous(limits = c(-6, 260),breaks= seq(0, 260, by=20)) + 
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.title = element_text(size=11), axis.text=element_text(size=6)) + 
  scale_color_manual(values = c("red", "blue")) + 
  geom_vline(xintercept = 200, linetype=2) + 
  geom_vline(xintercept = 100, linetype=2) +
  geom_vline(xintercept = 150, linetype=2)
  
  jpeg(paste0(resultDirPath, "LMvsKcorr_AGP.jpeg")
  print(compareStatModels_plot_agp)
  dev.off()

## gloor
compare_lm_gloor = lognorm_gloor_powerMatrix_LM
compare_kcorr_gloor = lognorm_gloor_powerMatrix_KCorr

compare_lm_gloor$test = rep("Parametric", nrow(compare_lm_gloor))
compare_kcorr_gloor$test = rep("Non-Parametric", nrow(compare_kcorr_gloor))

compare_merge_gloor = rbind(compare_lm_gloor, compare_kcorr_gloor)

compareStatModels_plot_gloor = ggplot(compare_merge_gloor, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = compare_merge_gloor, aes(x= Subsamplesize, y = Mean, color=test)) +
  geom_errorbar(aes(ymin = Mean - SD, 
                    ymax = Mean + SD, color=test)) +
  ggtitle("Gloor(Non-Parametric vs Parametric) n=803") +
  scale_x_continuous(trans = 'log10') + #scale_y_continuous(limits = c(-6, 260),breaks= seq(0, 260, by=20)) + 
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.title = element_text(size=11), axis.text=element_text(size=6)) + 
  scale_color_manual(values = c("red", "blue")) + 
  geom_vline(xintercept = 200, linetype=2) + 
  geom_vline(xintercept = 100, linetype=2) +
  geom_vline(xintercept = 150, linetype=2)

  jpeg(paste0(resultDirPath, "LMvsKcorr_Gloor.jpeg")
  print(compareStatModels_plot_gloor)
  dev.off()

## goodrich

compare_lm_goodrich = lognorm_goodrich_powerMatrix_LM
compare_kcorr_goodrich = lognorm_goodrich_powerMatrix_KCorr

compare_lm_goodrich$test = rep("Parametric", nrow(compare_lm_goodrich))
compare_kcorr_goodrich$test = rep("Non-Parametric", nrow(compare_kcorr_goodrich))

compare_merge_goodrich = rbind(compare_lm_goodrich, compare_kcorr_goodrich)

compareStatModels_plot_goodrich = ggplot(compare_merge_goodrich, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = compare_merge_goodrich, aes(x= Subsamplesize, y = Mean, color=test)) +
  geom_errorbar(aes(ymin = Mean - SD, 
                    ymax = Mean + SD, color=test)) +
  ggtitle("Goodrich(Non-Parametric vs Parametric) n=835") +
  scale_x_continuous(trans = 'log10') + scale_y_continuous(breaks= seq(0, 260, by=20)) + 
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.title = element_text(size=11), axis.text=element_text(size=6)) + 
  scale_color_manual(values = c("red", "blue")) + 
  geom_vline(xintercept = 200, linetype=2) + 
  geom_vline(xintercept = 100, linetype=2) +
  geom_vline(xintercept = 150, linetype=2)

  jpeg(paste0(resultDirPath, "LMvsKcorr_Goodrich.jpeg")
  print(compareStatModels_plot_goodrich)
  dev.off() 

## baxter
compare_lm_baxter = lognorm_baxter_powerMatrix_LM
compare_kcorr_baxter = lognorm_baxter_powerMatrix_KCorr

compare_lm_baxter$test = rep("Parametric", nrow(compare_lm_baxter))
compare_kcorr_baxter$test = rep("Non-Parametric", nrow(compare_kcorr_baxter))

compare_merge_baxter = rbind(compare_lm_baxter, compare_kcorr_baxter)

compareStatModels_plot_baxter = ggplot(compare_merge_baxter, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = compare_merge_baxter, aes(x= Subsamplesize, y = Mean, color=test)) +
  geom_errorbar(aes(ymin = Mean - SD, 
                    ymax = Mean + SD, color=test)) +
  ggtitle("Baxter(Non-Parametric vs Parametric) n=490") +
  scale_x_continuous(trans = 'log10') + 
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.title = element_text(size=11), axis.text=element_text(size=6)) + 
  scale_color_manual(values = c("red", "blue")) + 
  geom_vline(xintercept = 200, linetype=2) + 
  geom_vline(xintercept = 100, linetype=2) +
  geom_vline(xintercept = 150, linetype=2)

  jpeg(paste0(resultDirPath, "LMvsKcorr_Baxter.jpeg")
  print(compareStatModels_plot_baxter)
  dev.off() 

## morgan

compare_lm_morgan = lognorm_morgan_powerMatrix_LM
compare_kcorr_morgan = lognorm_morgan_powerMatrix_KCorr

compare_lm_morgan$test = rep("Parametric", nrow(compare_lm_morgan))
compare_kcorr_morgan$test = rep("Non-Parametric", nrow(compare_kcorr_morgan))

compare_merge_morgan = rbind(compare_lm_morgan, compare_kcorr_morgan)

compareStatModels_plot_morgan = ggplot(compare_merge_morgan, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = compare_merge_morgan, aes(x= Subsamplesize, y = Mean, color=test)) +
  geom_errorbar(aes(ymin = Mean - SD, 
                    ymax = Mean + SD, color=test)) +
  ggtitle("Morgan(Non-Parametric vs Parametric) n=228") +
  scale_x_continuous(trans = 'log10') + #scale_y_continuous(limits = c(-6, 260),breaks= seq(0, 260, by=20)) + 
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.title = element_text(size=11), axis.text=element_text(size=6)) + 
  scale_color_manual(values = c("red", "blue")) + 
  geom_vline(xintercept = 200, linetype=2) + 
  geom_vline(xintercept = 100, linetype=2) +
  geom_vline(xintercept = 150, linetype=2)
       
  jpeg(paste0(resultDirPath, "LMvsKcorr_Morgan.jpeg")
  print(compareStatModels_plot_morgan)
  dev.off()      

## nogbcn0

compare_lm_nogbcn0 = lognorm_nogbcn0_powerMatrix_LM
compare_kcorr_nogbcn0 = lognorm_nogbcn0_powerMatrix_KCorr

compare_lm_nogbcn0$test = rep("Parametric", nrow(compare_lm_nogbcn0))
compare_kcorr_nogbcn0$test = rep("Non-Parametric", nrow(compare_kcorr_nogbcn0))

compare_merge_nogbcn0 = rbind(compare_lm_nogbcn0, compare_kcorr_nogbcn0)

compareStatModels_plot_nogbcn0 = ggplot(compare_merge_nogbcn0, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = compare_merge_nogbcn0, aes(x= Subsamplesize, y = Mean, color=test)) +
  geom_errorbar(aes(ymin = Mean - SD, 
                    ymax = Mean + SD, color=test)) +
  ggtitle("NogBCN0(Non-Parametric vs Parametric) n=156") +
  scale_x_continuous(trans = 'log10') + #scale_y_continuous(limits = c(-6, 260),breaks= seq(0, 260, by=20)) + 
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.title = element_text(size=11), axis.text=element_text(size=6)) + 
  scale_color_manual(values = c("red", "blue")) + 
  geom_vline(xintercept = 200, linetype=2) + 
  geom_vline(xintercept = 100, linetype=2) +
  geom_vline(xintercept = 150, linetype=2)
       
  jpeg(paste0(resultDirPath, "LMvsKcorr_NogBCN0.jpeg")
  print(compareStatModels_plot_nogbcn0)
  dev.off()         
}

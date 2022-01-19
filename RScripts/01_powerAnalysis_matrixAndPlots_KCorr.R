## using lognormfiles_q2021 output from main section....
## using KENDALL CORR for power analysis

lognorm_output_agp = lognormfiles_q2021[[11]]
lognorm_output_agp_taxa = lognormfiles_q2021[[11]][,-ncol(lognormfiles_q2021[[11]])]
lognorm_output_agp_age = lognormfiles_q2021[[11]][,ncol(lognormfiles_q2021[[11]])]

agp_intervals = c(20,40,60,80, 100, 120, 140, 160, 180, 
                  200, 300,400, 500, 600, 700, 800, 900, 1000,
                  1100, 1200, 1368)

lognorm_agp_powerMatrix_KCorr = matrix(nrow=length(agp_intervals), ncol = 3)

for(i in 1:length(agp_intervals)){
  lognorm_agp_powerMatrix_KCorr[i,] = unlist(repeatSubsampling_corr_function(maxiteration = 50,
                                                                             subsamplesize = agp_intervals[i],
                                                                             lognormFile = lognorm_output_agp_taxa,
                                                                             metadata = lognorm_output_agp_age,
                                                                             corrMethod="kendall"))
}

colnames(lognorm_agp_powerMatrix_KCorr) = c("Subsamplesize", "Mean", "SD")

lognorm_agp_powerMatrix_KCorr = as.data.frame(lognorm_agp_powerMatrix_KCorr)

powerPlot_agp_KCorr = ggplot(lognorm_agp_powerMatrix_KCorr, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = lognorm_agp_powerMatrix_KCorr, aes(x= Subsamplesize, y = Mean)) +
  geom_errorbar(aes(ymin = Mean - SD, 
                    ymax = Mean + SD)) +
  ggtitle("AGP (n=1,368)") +
  scale_x_continuous(trans = 'log10') + scale_y_continuous(limits=c(-6,260),breaks= seq(0, 260, by=20)) + 
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.title = element_text(size=11), axis.text=element_text(size=6))


####_______________________________________________________ROSS____________________________

lognorm_output_ross = lognormfiles_q2021[[10]]
lognorm_output_ross_taxa = lognormfiles_q2021[[10]][,-ncol(lognormfiles_q2021[[10]])]
lognorm_output_ross_age = lognormfiles_q2021[[10]][,ncol(lognormfiles_q2021[[10]])]

ross_intervals = c(10,20,30,40,50,63)

lognorm_ross_powerMatrix_KCorr = matrix(nrow=length(ross_intervals), ncol = 3)

for(i in 1:length(ross_intervals)){
  lognorm_ross_powerMatrix_KCorr[i,] = unlist(repeatSubsampling_corr_function(maxiteration = 50,
                                                                              subsamplesize = ross_intervals[i],
                                                                              lognormFile = lognorm_output_ross_taxa,
                                                                              metadata = lognorm_output_ross_age,
                                                                              corrMethod="kendall"))
}

colnames(lognorm_ross_powerMatrix_KCorr) = c("Subsamplesize", "Mean", "SD")

lognorm_ross_powerMatrix_KCorr = as.data.frame(lognorm_ross_powerMatrix_KCorr)


powerPlot_ross_KCorr = ggplot(lognorm_ross_powerMatrix_KCorr, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = lognorm_ross_powerMatrix_KCorr, aes(x= Subsamplesize, y = Mean)) +
  geom_errorbar(aes(ymin = Mean - SD, 
                    ymax = Mean + SD)) +
  ggtitle("Ross (n=63)") +
  scale_x_continuous(trans = 'log10') + scale_y_continuous(limits= c(-6,35), breaks= seq(0, 260, by=5)) + 
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.title = element_text(size=11))



##_______________________________________________________NOG_STK__________________________________________
lognorm_output_nogstk = lognormfiles_q2021[[9]]
lognorm_output_nogstk_taxa = lognormfiles_q2021[[9]][,-ncol(lognormfiles_q2021[[9]])]
lognorm_output_nogstk_age = lognormfiles_q2021[[9]][,ncol(lognormfiles_q2021[[9]])]

nogstk_intervals = c(10,20,30,40,50,60,70,84)

lognorm_nogstk_powerMatrix_KCorr = matrix(nrow=length(nogstk_intervals), ncol = 3)

for(i in 1:length(nogstk_intervals)){
  lognorm_nogstk_powerMatrix_KCorr[i,] = unlist(repeatSubsampling_corr_function(maxiteration = 50,
                                                                                subsamplesize = nogstk_intervals[i],
                                                                                lognormFile = lognorm_output_nogstk_taxa,
                                                                                metadata = lognorm_output_nogstk_age,
                                                                                corrMethod="kendall"))
}

colnames(lognorm_nogstk_powerMatrix_KCorr) = c("Subsamplesize", "Mean", "SD")
lognorm_nogstk_powerMatrix_KCorr = as.data.frame(lognorm_nogstk_powerMatrix_KCorr)

powerPlot_nogstk_KCorr = ggplot(lognorm_nogstk_powerMatrix_KCorr, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = lognorm_nogstk_powerMatrix_KCorr, aes(x= Subsamplesize, y = Mean)) +
  geom_errorbar(aes(ymin = Mean - SD, 
                    ymax = Mean + SD)) +
  ggtitle("NogSTK (n=84)") +
  scale_x_continuous(trans = 'log10') + scale_y_continuous(limits= c(-6,35), breaks= seq(0, 260, by=5)) + 
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.title = element_text(size=11))


##_______________________________________________________NOG_BCN0__________________________________________
lognorm_output_nogbcn0 = lognormfiles_q2021[[8]]
lognorm_output_nogbcn0_taxa = lognormfiles_q2021[[8]][,-ncol(lognormfiles_q2021[[8]])]
lognorm_output_nogbcn0_age = lognormfiles_q2021[[8]][,ncol(lognormfiles_q2021[[8]])]

nogbcn0_intervals = c(20,40,60,80,100,120,140,156)

lognorm_nogbcn0_powerMatrix_KCorr = matrix(nrow=length(nogbcn0_intervals), ncol = 3)

for(i in 1:length(nogbcn0_intervals)){
  lognorm_nogbcn0_powerMatrix_KCorr[i,] = unlist(repeatSubsampling_corr_function(maxiteration = 50,
                                                                                 subsamplesize = nogbcn0_intervals[i],
                                                                                 lognormFile = lognorm_output_nogbcn0_taxa,
                                                                                 metadata = lognorm_output_nogbcn0_age,
                                                                                 corrMethod="kendall"))
}

colnames(lognorm_nogbcn0_powerMatrix_KCorr) = c("Subsamplesize", "Mean", "SD")
lognorm_nogbcn0_powerMatrix_KCorr = as.data.frame(lognorm_nogbcn0_powerMatrix_KCorr)

powerPlot_nogbcn0_KCorr = ggplot(lognorm_nogbcn0_powerMatrix_KCorr, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = lognorm_nogbcn0_powerMatrix_KCorr, aes(x= Subsamplesize, y = Mean)) +
  geom_errorbar(aes(ymin = Mean - SD, 
                    ymax = Mean + SD)) +
  ggtitle("NogBCN0 (n=156)") +
  scale_x_continuous(trans = 'log10') + scale_y_continuous(limits= c(-6,35.5), breaks= seq(0, 260, by=5)) + 
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.title = element_text(size=11))


##_______________________________________________________GOODRICH__________________________________________
lognorm_output_goodrich = lognormfiles_q2021[[7]]
lognorm_output_goodrich_taxa = lognormfiles_q2021[[7]][,-ncol(lognormfiles_q2021[[7]])]
lognorm_output_goodrich_age = lognormfiles_q2021[[7]][,ncol(lognormfiles_q2021[[7]])]

goodrich_intervals = c(20,40,60,80,100,120,140,160,180,
                       200,300,400,500,600,700,800,835)

lognorm_goodrich_powerMatrix_KCorr = matrix(nrow=length(goodrich_intervals), ncol = 3)

for(i in 1:length(goodrich_intervals)){
  lognorm_goodrich_powerMatrix_KCorr[i,] = unlist(repeatSubsampling_corr_function(maxiteration = 50,
                                                                                  subsamplesize = goodrich_intervals[i],
                                                                                  lognormFile = lognorm_output_goodrich_taxa,
                                                                                  metadata = lognorm_output_goodrich_age,
                                                                                  corrMethod="kendall"))
}

colnames(lognorm_goodrich_powerMatrix_KCorr) = c("Subsamplesize", "Mean", "SD")
lognorm_goodrich_powerMatrix_KCorr = as.data.frame(lognorm_goodrich_powerMatrix_KCorr)

powerPlot_goodrich_KCorr = ggplot(lognorm_goodrich_powerMatrix_KCorr, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = lognorm_goodrich_powerMatrix_KCorr, aes(x= Subsamplesize, y = Mean)) +
  geom_errorbar(aes(ymin = Mean - SD, 
                    ymax = Mean + SD)) +
  ggtitle("Goodrich (n=835)") +
  scale_x_continuous(trans = 'log10') + scale_y_continuous(limits=c(-6,35),breaks= seq(0, 260, by=5)) + 
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.title = element_text(size=11), axis.text=element_text(size=6))


##_______________________________________________________GLOOR__________________________________________

lognorm_output_gloor = lognormfiles_q2021[[6]]
lognorm_output_gloor_taxa = lognormfiles_q2021[[6]][,-ncol(lognormfiles_q2021[[6]])]
lognorm_output_gloor_age = lognormfiles_q2021[[6]][,ncol(lognormfiles_q2021[[6]])]

gloor_intervals = c(20,40,60,80,100,120,140,160,180,
                    200,300,400,500,600,700,803)

lognorm_gloor_powerMatrix_KCorr = matrix(nrow=length(gloor_intervals), ncol = 3)

for(i in 1:length(gloor_intervals)){
  lognorm_gloor_powerMatrix_KCorr[i,] = unlist(repeatSubsampling_corr_function(maxiteration = 50,
                                                                               subsamplesize = gloor_intervals[i],
                                                                               lognormFile = lognorm_output_gloor_taxa,
                                                                               metadata = lognorm_output_gloor_age,
                                                                               corrMethod="kendall"))
}

colnames(lognorm_gloor_powerMatrix_KCorr) = c("Subsamplesize", "Mean", "SD")
lognorm_gloor_powerMatrix_KCorr = as.data.frame(lognorm_gloor_powerMatrix_KCorr)

powerPlot_gloor_KCorr = ggplot(lognorm_gloor_powerMatrix_KCorr, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = lognorm_gloor_powerMatrix_KCorr, aes(x= Subsamplesize, y = Mean)) +
  geom_errorbar(aes(ymin = Mean - SD, 
                    ymax = Mean + SD)) +
  ggtitle("Gloor (n=803)") +
  scale_x_continuous(trans = 'log10') + scale_y_continuous(breaks= seq(0, 260, by=20)) + 
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.title = element_text(size=11), axis.text=element_text(size=6))

##_______________________________________________________ESCOBAR_________________________________________

lognorm_output_escobar = lognormfiles_q2021[[5]]
lognorm_output_escobar_taxa = lognormfiles_q2021[[5]][,-ncol(lognormfiles_q2021[[5]])]
lognorm_output_escobar_age = lognormfiles_q2021[[5]][,ncol(lognormfiles_q2021[[5]])]

escobar_intervals = c(10,20,30)

lognorm_escobar_powerMatrix_KCorr = matrix(nrow=length(escobar_intervals), ncol = 3)

for(i in 1:length(escobar_intervals)){
  lognorm_escobar_powerMatrix_KCorr[i,] = unlist(repeatSubsampling_corr_function(maxiteration = 50,
                                                                                 subsamplesize = escobar_intervals[i],
                                                                                 lognormFile = lognorm_output_escobar_taxa,
                                                                                 metadata = lognorm_output_escobar_age,
                                                                                 corrMethod="kendall"))
}

colnames(lognorm_escobar_powerMatrix_KCorr) = c("Subsamplesize", "Mean", "SD")
lognorm_escobar_powerMatrix_KCorr = as.data.frame(lognorm_escobar_powerMatrix_KCorr)

powerPlot_escobar_KCorr = ggplot(lognorm_escobar_powerMatrix_KCorr, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = lognorm_escobar_powerMatrix_KCorr, aes(x= Subsamplesize, y = Mean)) +
  geom_errorbar(aes(ymin = Mean - SD, 
                    ymax = Mean + SD)) +
  ggtitle("Escobar (n=30)") +
  scale_x_continuous(trans = 'log10') + scale_y_continuous(limits= c(-6,35), breaks= seq(0, 260, by=5)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        plot.title = element_text(size=11)) 

##_______________________________________________________BAXTER_________________________________________
lognorm_output_baxter = lognormfiles_q2021[[4]]
lognorm_output_baxter_taxa = lognormfiles_q2021[[4]][,-ncol(lognormfiles_q2021[[4]])]
lognorm_output_baxter_age = lognormfiles_q2021[[4]][,ncol(lognormfiles_q2021[[4]])]

baxter_intervals = c(20,40,60,80,100,120,140,160,180,
                     200,220,240,260,280,300, 320, 340,
                     360,380,400, 420,440,460,490)

lognorm_baxter_powerMatrix_KCorr = matrix(nrow=length(baxter_intervals), ncol = 3)

for(i in 1:length(baxter_intervals)){
  lognorm_baxter_powerMatrix_KCorr[i,] = unlist(repeatSubsampling_corr_function(maxiteration = 50,
                                                                                subsamplesize = baxter_intervals[i],
                                                                                lognormFile = lognorm_output_baxter_taxa,
                                                                                metadata = lognorm_output_baxter_age,
                                                                                corrMethod="kendall"))
}

colnames(lognorm_baxter_powerMatrix_KCorr) = c("Subsamplesize", "Mean", "SD")
lognorm_baxter_powerMatrix_KCorr = as.data.frame(lognorm_baxter_powerMatrix_KCorr)

powerPlot_baxter_KCorr = ggplot(lognorm_baxter_powerMatrix_KCorr, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = lognorm_baxter_powerMatrix_KCorr, aes(x= Subsamplesize, y = Mean)) +
  geom_errorbar(aes(ymin = Mean - SD, 
                    ymax = Mean + SD)) +
  ggtitle("Baxter (n=490)") +
  scale_x_continuous(trans = 'log10') + scale_y_continuous(limits= c(-6,35), breaks= seq(0, 260, by=5)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        plot.title = element_text(size=11)) 

   ### baxter ---disease subset -----
disease.baxter=subset(baxter_status_lognorm, baxter_status_lognorm$status=="case")
disease.baxter.age = disease.baxter$Age
disease.baxter.genusTaxa = disease.baxter[,grep("g__", colnames(disease.baxter))]

disease_baxter_intervals = c(20,40,60,80,100,120,140,160,180,
                     200,220,240,260,280,300, 318)

disease_baxter_powerMatrix_KCorr = matrix(nrow=length(disease_baxter_intervals), ncol = 3)

for(i in 1:length(disease_baxter_intervals)){
  disease_baxter_powerMatrix_KCorr[i,] = unlist(repeatSubsampling_corr_function(maxiteration = 50,
                                                                                subsamplesize = disease_baxter_intervals[i],
                                                                                lognormFile = disease.baxter.genusTaxa,
                                                                                metadata = disease.baxter.age,
                                                                                corrMethod="kendall"))
}

colnames(disease_baxter_powerMatrix_KCorr) = c("Subsamplesize", "Mean", "SD")
disease_baxter_powerMatrix_KCorr = as.data.frame(disease_baxter_powerMatrix_KCorr)

powerPlot_baxter_KCorr_disease = ggplot(disease_baxter_powerMatrix_KCorr, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = disease_baxter_powerMatrix_KCorr, aes(x= Subsamplesize, y = Mean)) +
  geom_errorbar(aes(ymin = Mean - SD, 
                    ymax = Mean + SD)) +
  ggtitle("Baxter, Disease (n=318)") +
  scale_x_continuous(trans = 'log10') + scale_y_continuous(limits= c(-6,35), breaks= seq(0, 260, by=5)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        plot.title = element_text(size=11)) 

### baxter ---healthy subset -----

healthy.baxter=subset(baxter_status_lognorm, baxter_status_lognorm$status=="control")
healthy.baxter.age = healthy.baxter$Age
healthy.baxter.genusTaxa = healthy.baxter[,grep("g__", colnames(healthy.baxter))]
healthy_baxter_intervals = c(20,40,60,80,100,120,140,160,172)
healthy_baxter_powerMatrix_KCorr = matrix(nrow=length(healthy_baxter_intervals), ncol = 3)

for(i in 1:length(healthy_baxter_intervals)){
  healthy_baxter_powerMatrix_KCorr[i,] = unlist(repeatSubsampling_corr_function(maxiteration = 50,
                                                                                subsamplesize = healthy_baxter_intervals[i],
                                                                                lognormFile = healthy.baxter.genusTaxa,
                                                                                metadata = healthy.baxter.age,
                                                                                corrMethod="kendall"))
}

colnames(healthy_baxter_powerMatrix_KCorr) = c("Subsamplesize", "Mean", "SD")
healthy_baxter_powerMatrix_KCorr = as.data.frame(healthy_baxter_powerMatrix_KCorr)

powerPlot_baxter_KCorr_healthy = ggplot(healthy_baxter_powerMatrix_KCorr, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = healthy_baxter_powerMatrix_KCorr, aes(x= Subsamplesize, y = Mean)) +
  geom_errorbar(aes(ymin = Mean - SD, 
                    ymax = Mean + SD)) +
  ggtitle("Baxter, Healthy (n=172)") +
  scale_x_continuous(trans = 'log10') + #scale_y_continuous(limits= c(-6,35), breaks= seq(0, 260, by=5)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        plot.title = element_text(size=11)) 





##_______________________________________________________ZELLERGERMANY_________________________________________
lognorm_output_zellergermany = lognormfiles_q2021[[3]]
lognorm_output_zellergermany_taxa = lognormfiles_q2021[[3]][,-ncol(lognormfiles_q2021[[3]])]
lognorm_output_zellergermany_age = lognormfiles_q2021[[3]][,ncol(lognormfiles_q2021[[3]])]

zellergermany_intervals = c(10,20,30,40,48)

lognorm_zellergermany_powerMatrix_KCorr = matrix(nrow=length(zellergermany_intervals), ncol = 3)

for(i in 1:length(zellergermany_intervals)){
  lognorm_zellergermany_powerMatrix_KCorr[i,] = unlist(repeatSubsampling_corr_function(maxiteration = 50,
                                                                                       subsamplesize = zellergermany_intervals[i],
                                                                                       lognormFile = lognorm_output_zellergermany_taxa,
                                                                                       metadata = lognorm_output_zellergermany_age,
                                                                                       corrMethod="kendall"))
}

colnames(lognorm_zellergermany_powerMatrix_KCorr) = c("Subsamplesize", "Mean", "SD")
lognorm_zellergermany_powerMatrix_KCorr = as.data.frame(lognorm_zellergermany_powerMatrix_KCorr)

powerPlot_zellergermany_KCorr = ggplot(lognorm_zellergermany_powerMatrix_KCorr, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = lognorm_zellergermany_powerMatrix_KCorr, aes(x= Subsamplesize, y = Mean)) +
  geom_errorbar(aes(ymin = Mean - SD, 
                    ymax = Mean + SD)) +
  ggtitle("ZellerGermany (n=48)") +
  scale_x_continuous(trans = 'log10') + scale_y_continuous(limits= c(-6,35), breaks= seq(0, 260, by=5)) + 
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.title = element_text(size=11))

##_______________________________________________________ZELLERFRANCE_____________
lognorm_output_zellerfrance = lognormfiles_q2021[[2]]
lognorm_output_zellerfrance_taxa = lognormfiles_q2021[[2]][,-ncol(lognormfiles_q2021[[2]])]
lognorm_output_zellerfrance_age = lognormfiles_q2021[[2]][,ncol(lognormfiles_q2021[[2]])]

zellerfrance_intervals = c(10,20,30,40,50,60,70,80,90,100,
                           120,129)

lognorm_zellerfrance_powerMatrix_KCorr = matrix(nrow=length(zellerfrance_intervals), ncol = 3)

for(i in 1:length(zellerfrance_intervals)){
  lognorm_zellerfrance_powerMatrix_KCorr[i,] = unlist(repeatSubsampling_corr_function(maxiteration = 50,
                                                                                      subsamplesize = zellerfrance_intervals[i],
                                                                                      lognormFile = lognorm_output_zellerfrance_taxa,
                                                                                      metadata = lognorm_output_zellerfrance_age,
                                                                                      corrMethod="kendall"))
}

colnames(lognorm_zellerfrance_powerMatrix_KCorr) = c("Subsamplesize", "Mean", "SD")
lognorm_zellerfrance_powerMatrix_KCorr = as.data.frame(lognorm_zellerfrance_powerMatrix_KCorr)

powerPlot_zellerfrance_KCorr = ggplot(lognorm_zellerfrance_powerMatrix_KCorr, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = lognorm_zellerfrance_powerMatrix_KCorr, aes(x= Subsamplesize, y = Mean)) +
  geom_errorbar(aes(ymin = Mean - SD, 
                    ymax = Mean + SD)) +
  ggtitle("ZellerFrance (n=129)") +
  scale_x_continuous(trans = 'log10') + scale_y_continuous(limits= c(-6,35), breaks= seq(0, 260, by=5)) + 
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.title = element_text(size=11))

##____________________________________________________________MORGAN________________________
lognorm_output_morgan = lognormfiles_q2021[[1]]
lognorm_output_morgan_taxa = lognormfiles_q2021[[1]][,-ncol(lognormfiles_q2021[[1]])]
lognorm_output_morgan_age = lognormfiles_q2021[[1]][,ncol(lognormfiles_q2021[[1]])]

morgan_intervals = c(10,20,30,40,50,60,70,80,90,100,
                     120,140,160,180,200,228)

lognorm_morgan_powerMatrix_KCorr = matrix(nrow=length(morgan_intervals), ncol = 3)

for(i in 1:length(morgan_intervals)){
  lognorm_morgan_powerMatrix_KCorr[i,] = unlist(repeatSubsampling_corr_function(maxiteration = 50,
                                                                                subsamplesize = morgan_intervals[i],
                                                                                lognormFile = lognorm_output_morgan_taxa,
                                                                                metadata = lognorm_output_morgan_age,
                                                                                corrMethod="kendall"))
}

colnames(lognorm_morgan_powerMatrix_KCorr) = c("Subsamplesize", "Mean", "SD")
lognorm_morgan_powerMatrix_KCorr = as.data.frame(lognorm_morgan_powerMatrix_KCorr)

powerPlot_morgan_KCorr = ggplot(lognorm_morgan_powerMatrix_KCorr, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = lognorm_morgan_powerMatrix_KCorr, aes(x= Subsamplesize, y = Mean)) +
  geom_errorbar(aes(ymin = Mean - SD, 
                    ymax = Mean + SD)) +
  ggtitle("Morgan (n=228)") +
  scale_x_continuous(trans = 'log10') + scale_y_continuous(limits= c(-6,35), breaks= seq(0, 260, by=5)) + 
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.title = element_text(size=11))


### combine all plots

fullsets_comparison_KCorr = ggarrange(powerPlot_nogbcn0_KCorr,powerPlot_nogstk_KCorr, powerPlot_escobar_KCorr,
                                      plot.new(), plot.new(), 
                                      powerPlot_zellerfrance_KCorr, powerPlot_zellergermany_KCorr,powerPlot_goodrich_KCorr, powerPlot_baxter_KCorr,powerPlot_ross_KCorr, 
                                      powerPlot_gloor_KCorr,powerPlot_morgan_KCorr, powerPlot_agp_KCorr, plot.new(), plot.new(), 
                                      labels = c(LETTERS[1:3], "", "", LETTERS[4:11]),
                                      nrow = 3, ncol = 5) 

annotate_figure(fullsets_comparison_KCorr, bottom="Sample size", 
                left="Number of Significant Taxa")

fullsets_comparison_KCorr_list = list(lognorm_morgan_powerMatrix_KCorr, lognorm_zellerfrance_powerMatrix_KCorr, lognorm_zellergermany_powerMatrix_KCorr,
                                      lognorm_baxter_powerMatrix_KCorr, lognorm_escobar_powerMatrix_KCorr, lognorm_gloor_powerMatrix_KCorr, 
                                      lognorm_goodrich_powerMatrix_KCorr, lognorm_nogbcn0_powerMatrix_KCorr, lognorm_nogstk_powerMatrix_KCorr,
                                      lognorm_ross_powerMatrix_KCorr, lognorm_agp_powerMatrix_KCorr)

fullsets_comparison_KCorr_list[[11]]


for(i in 1:11){
  outfileName = paste0("/Users/anhil/Desktop/", coh_l[[i]], "_powerResults_kcorr.txt")
  write.table(fullsets_comparison_KCorr_list[[i]], outfileName, quote = F, row.names = F)
}













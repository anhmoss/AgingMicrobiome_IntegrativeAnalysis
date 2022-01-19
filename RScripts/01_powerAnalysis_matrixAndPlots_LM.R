## using lognormfiles_q2021 output from main section....

lognorm_output_agp = lognormfiles_q2021[[11]]
lognorm_output_agp_taxa = lognormfiles_q2021[[11]][,-ncol(lognormfiles_q2021[[11]])]
lognorm_output_agp_age = lognormfiles_q2021[[11]][,ncol(lognormfiles_q2021[[11]])]

agp_intervals = c(20,40,60,80, 100, 120, 140, 160, 180, 
                  200, 300,400, 500, 600, 700, 800, 900, 1000,
                  1100, 1200, 1368)

lognorm_agp_powerMatrix_LM = matrix(nrow=length(agp_intervals), ncol = 3)

for(i in 1:length(agp_intervals)){
  lognorm_agp_powerMatrix_LM[i,] = unlist(repeatSubsampling_LM_function(maxiteration = 50,
                                                                     subsamplesize = agp_intervals[i],
                                                                     lognormFile = lognorm_output_agp_taxa,
                                                                     metadata = lognorm_output_agp_age))
}

colnames(lognorm_agp_powerMatrix_LM) = c("Subsamplesize", "Mean", "SD")

lognorm_agp_powerMatrix_LM = as.data.frame(lognorm_agp_powerMatrix_LM)

powerPlot_agp_LM = ggplot(lognorm_agp_powerMatrix_LM, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = lognorm_agp_powerMatrix_LM, aes(x= Subsamplesize, y = Mean)) +
  geom_errorbar(aes(ymin = Mean - SD, 
                    ymax = Mean + SD)) +
  ggtitle("AGP (LM) n=1,368") +
  scale_x_continuous(trans = 'log10') + scale_y_continuous(limits = c(-6, 260),breaks= seq(0, 260, by=20)) + 
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.title = element_text(size=11), axis.text=element_text(size=6))


####_______________________________________________________ROSS____________________________

lognorm_output_ross = lognormfiles_q2021[[10]]
lognorm_output_ross_taxa = lognormfiles_q2021[[10]][,-ncol(lognormfiles_q2021[[10]])]
lognorm_output_ross_age = lognormfiles_q2021[[10]][,ncol(lognormfiles_q2021[[10]])]

ross_intervals = c(10,20,30,40,50,63)

lognorm_ross_powerMatrix_LM = matrix(nrow=length(ross_intervals), ncol = 3)

for(i in 1:length(ross_intervals)){
  lognorm_ross_powerMatrix_LM[i,] = unlist(repeatSubsampling_LM_function(maxiteration = 50,
                                                                         subsamplesize = ross_intervals[i],
                                                                         lognormFile = lognorm_output_ross_taxa,
                                                                         metadata = lognorm_output_ross_age))
}

colnames(lognorm_ross_powerMatrix_LM) = c("Subsamplesize", "Mean", "SD")

lognorm_ross_powerMatrix_LM = as.data.frame(lognorm_ross_powerMatrix_LM)


powerPlot_ross_LM = ggplot(lognorm_ross_powerMatrix_LM, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = lognorm_ross_powerMatrix_LM, aes(x= Subsamplesize, y = Mean)) +
  geom_errorbar(aes(ymin = Mean - SD, 
                    ymax = Mean + SD)) +
  ggtitle("Ross (LM) n=63") +
  scale_x_continuous(trans = 'log10') + scale_y_continuous(limits= c(-6,35), breaks= seq(0, 260, by=5)) + 
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.title = element_text(size=11))



##_______________________________________________________NOG_STK__________________________________________
lognorm_output_nogstk = lognormfiles_q2021[[9]]
lognorm_output_nogstk_taxa = lognormfiles_q2021[[9]][,-ncol(lognormfiles_q2021[[9]])]
lognorm_output_nogstk_age = lognormfiles_q2021[[9]][,ncol(lognormfiles_q2021[[9]])]

nogstk_intervals = c(10,20,30,40,50,60,70,84)

lognorm_nogstk_powerMatrix_LM = matrix(nrow=length(nogstk_intervals), ncol = 3)

for(i in 1:length(nogstk_intervals)){
  lognorm_nogstk_powerMatrix_LM[i,] = unlist(repeatSubsampling_LM_function(maxiteration = 50,
                                                                           subsamplesize = nogstk_intervals[i],
                                                                           lognormFile = lognorm_output_nogstk_taxa,
                                                                           metadata = lognorm_output_nogstk_age))
}

colnames(lognorm_nogstk_powerMatrix_LM) = c("Subsamplesize", "Mean", "SD")
lognorm_nogstk_powerMatrix_LM = as.data.frame(lognorm_nogstk_powerMatrix_LM)

powerPlot_nogstk_LM = ggplot(lognorm_nogstk_powerMatrix_LM, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = lognorm_nogstk_powerMatrix_LM, aes(x= Subsamplesize, y = Mean)) +
  geom_errorbar(aes(ymin = Mean - SD, 
                    ymax = Mean + SD)) +
  ggtitle("NogSTK (LM) n=84") +
  scale_x_continuous(trans = 'log10') + scale_y_continuous(limits= c(-6,35), breaks= seq(0, 260, by=5)) + 
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.title = element_text(size=11))


##_______________________________________________________NOG_BCN0__________________________________________
lognorm_output_nogbcn0 = lognormfiles_q2021[[8]]
lognorm_output_nogbcn0_taxa = lognormfiles_q2021[[8]][,-ncol(lognormfiles_q2021[[8]])]
lognorm_output_nogbcn0_age = lognormfiles_q2021[[8]][,ncol(lognormfiles_q2021[[8]])]

nogbcn0_intervals = c(20,40,60,80,100,120,140,156)

lognorm_nogbcn0_powerMatrix_LM = matrix(nrow=length(nogbcn0_intervals), ncol = 3)

for(i in 1:length(nogbcn0_intervals)){
  lognorm_nogbcn0_powerMatrix_LM[i,] = unlist(repeatSubsampling_LM_function(maxiteration = 50,
                                                                            subsamplesize = nogbcn0_intervals[i],
                                                                            lognormFile = lognorm_output_nogbcn0_taxa,
                                                                            metadata = lognorm_output_nogbcn0_age))
}

colnames(lognorm_nogbcn0_powerMatrix_LM) = c("Subsamplesize", "Mean", "SD")
lognorm_nogbcn0_powerMatrix_LM = as.data.frame(lognorm_nogbcn0_powerMatrix_LM)

powerPlot_nogbcn0_LM = ggplot(lognorm_nogbcn0_powerMatrix_LM, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = lognorm_nogbcn0_powerMatrix_LM, aes(x= Subsamplesize, y = Mean)) +
  geom_errorbar(aes(ymin = Mean - SD, 
                    ymax = Mean + SD)) +
  ggtitle("NogBCN0 (LM) n=156") +
  scale_x_continuous(trans = 'log10') + scale_y_continuous(limits= c(-6,35), breaks= seq(0, 260, by=5)) + 
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.title = element_text(size=11))


##_______________________________________________________GOODRICH__________________________________________
lognorm_output_goodrich = lognormfiles_q2021[[7]]
lognorm_output_goodrich_taxa = lognormfiles_q2021[[7]][,-ncol(lognormfiles_q2021[[7]])]
lognorm_output_goodrich_age = lognormfiles_q2021[[7]][,ncol(lognormfiles_q2021[[7]])]

goodrich_intervals = c(20,40,60,80,100,120,140,160,180,
                       200,300,400,500,600,700,800,835)

lognorm_goodrich_powerMatrix_LM = matrix(nrow=length(goodrich_intervals), ncol = 3)

for(i in 1:length(goodrich_intervals)){
  lognorm_goodrich_powerMatrix_LM[i,] = unlist(repeatSubsampling_LM_function(maxiteration = 50,
                                                                             subsamplesize = goodrich_intervals[i],
                                                                             lognormFile = lognorm_output_goodrich_taxa,
                                                                             metadata = lognorm_output_goodrich_age))
}

colnames(lognorm_goodrich_powerMatrix_LM) = c("Subsamplesize", "Mean", "SD")
lognorm_goodrich_powerMatrix_LM = as.data.frame(lognorm_goodrich_powerMatrix_LM)

powerPlot_goodrich_LM = ggplot(lognorm_goodrich_powerMatrix_LM, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = lognorm_goodrich_powerMatrix_LM, aes(x= Subsamplesize, y = Mean)) +
  geom_errorbar(aes(ymin = Mean - SD, 
                    ymax = Mean + SD)) +
  ggtitle("Goodrich (LM) n=835") +
  scale_x_continuous(trans = 'log10') + scale_y_continuous(limits= c(-6,35),breaks= seq(0, 260, by=5)) + 
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.title = element_text(size=11), axis.text=element_text(size=6))


##_______________________________________________________GLOOR__________________________________________

lognorm_output_gloor = lognormfiles_q2021[[6]]
lognorm_output_gloor_taxa = lognormfiles_q2021[[6]][,-ncol(lognormfiles_q2021[[6]])]
lognorm_output_gloor_age = lognormfiles_q2021[[6]][,ncol(lognormfiles_q2021[[6]])]

gloor_intervals = c(20,40,60,80,100,120,140,160,180,
                    200,300,400,500,600,700,803)

lognorm_gloor_powerMatrix_LM = matrix(nrow=length(gloor_intervals), ncol = 3)

for(i in 1:length(gloor_intervals)){
  lognorm_gloor_powerMatrix_LM[i,] = unlist(repeatSubsampling_LM_function(maxiteration = 50,
                                                                          subsamplesize = gloor_intervals[i],
                                                                          lognormFile = lognorm_output_gloor_taxa,
                                                                          metadata = lognorm_output_gloor_age))
}

colnames(lognorm_gloor_powerMatrix_LM) = c("Subsamplesize", "Mean", "SD")
lognorm_gloor_powerMatrix_LM = as.data.frame(lognorm_gloor_powerMatrix_LM)

powerPlot_gloor_LM = ggplot(lognorm_gloor_powerMatrix_LM, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = lognorm_gloor_powerMatrix_LM, aes(x= Subsamplesize, y = Mean)) +
  geom_errorbar(aes(ymin = Mean - SD, 
                    ymax = Mean + SD)) +
  ggtitle("Gloor (LM) n=803") +
  scale_x_continuous(trans = 'log10') + scale_y_continuous(breaks= seq(0, 260, by=20)) + 
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.title = element_text(size=11), axis.text=element_text(size=6))

##_______________________________________________________ESCOBAR_________________________________________

lognorm_output_escobar = lognormfiles_q2021[[5]]
lognorm_output_escobar_taxa = lognormfiles_q2021[[5]][,-ncol(lognormfiles_q2021[[5]])]
lognorm_output_escobar_age = lognormfiles_q2021[[5]][,ncol(lognormfiles_q2021[[5]])]

escobar_intervals = c(10,20,30)

lognorm_escobar_powerMatrix_LM = matrix(nrow=length(escobar_intervals), ncol = 3)

for(i in 1:length(escobar_intervals)){
  lognorm_escobar_powerMatrix_LM[i,] = unlist(repeatSubsampling_LM_function(maxiteration = 50,
                                                                            subsamplesize = escobar_intervals[i],
                                                                            lognormFile = lognorm_output_escobar_taxa,
                                                                            metadata = lognorm_output_escobar_age))
}

colnames(lognorm_escobar_powerMatrix_LM) = c("Subsamplesize", "Mean", "SD")
lognorm_escobar_powerMatrix_LM = as.data.frame(lognorm_escobar_powerMatrix_LM)

powerPlot_escobar_LM = ggplot(lognorm_escobar_powerMatrix_LM, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = lognorm_escobar_powerMatrix_LM, aes(x= Subsamplesize, y = Mean)) +
  geom_errorbar(aes(ymin = Mean - SD, 
                    ymax = Mean + SD)) +
  ggtitle("Escobar (LM) n=30") +
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

lognorm_baxter_powerMatrix_LM = matrix(nrow=length(baxter_intervals), ncol = 3)

for(i in 1:length(baxter_intervals)){
  lognorm_baxter_powerMatrix_LM[i,] = unlist(repeatSubsampling_LM_function(maxiteration = 50,
                                                                           subsamplesize = baxter_intervals[i],
                                                                           lognormFile = lognorm_output_baxter_taxa,
                                                                           metadata = lognorm_output_baxter_age))
}

colnames(lognorm_baxter_powerMatrix_LM) = c("Subsamplesize", "Mean", "SD")
lognorm_baxter_powerMatrix_LM = as.data.frame(lognorm_baxter_powerMatrix_LM)

powerPlot_baxter_LM = ggplot(lognorm_baxter_powerMatrix_LM, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = lognorm_baxter_powerMatrix_LM, aes(x= Subsamplesize, y = Mean)) +
  geom_errorbar(aes(ymin = Mean - SD, 
                    ymax = Mean + SD)) +
  ggtitle("Baxter (LM) n=490") +
  scale_x_continuous(trans = 'log10') + scale_y_continuous(limits= c(-6,35), breaks= seq(0, 260, by=5)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        plot.title = element_text(size=11)) 

##_______________________________________________________ZELLERGERMANY_________________________________________
lognorm_output_zellergermany = lognormfiles_q2021[[3]]
lognorm_output_zellergermany_taxa = lognormfiles_q2021[[3]][,-ncol(lognormfiles_q2021[[3]])]
lognorm_output_zellergermany_age = lognormfiles_q2021[[3]][,ncol(lognormfiles_q2021[[3]])]

zellergermany_intervals = c(10,20,30,40,48)

lognorm_zellergermany_powerMatrix_LM = matrix(nrow=length(zellergermany_intervals), ncol = 3)

for(i in 1:length(zellergermany_intervals)){
  lognorm_zellergermany_powerMatrix_LM[i,] = unlist(repeatSubsampling_LM_function(maxiteration = 50,
                                                                                  subsamplesize = zellergermany_intervals[i],
                                                                                  lognormFile = lognorm_output_zellergermany_taxa,
                                                                                  metadata = lognorm_output_zellergermany_age))
}

colnames(lognorm_zellergermany_powerMatrix_LM) = c("Subsamplesize", "Mean", "SD")
lognorm_zellergermany_powerMatrix_LM = as.data.frame(lognorm_zellergermany_powerMatrix_LM)

powerPlot_zellergermany_LM = ggplot(lognorm_zellergermany_powerMatrix_LM, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = lognorm_zellergermany_powerMatrix_LM, aes(x= Subsamplesize, y = Mean)) +
  geom_errorbar(aes(ymin = Mean - SD, 
                    ymax = Mean + SD)) +
  ggtitle("ZellerGermany (LM) n=48") +
  scale_x_continuous(trans = 'log10') + scale_y_continuous(limits= c(-6,35), breaks= seq(0, 260, by=5)) + 
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.title = element_text(size=11))

##_______________________________________________________ZELLERFRANCE_____________
lognorm_output_zellerfrance = lognormfiles_q2021[[2]]
lognorm_output_zellerfrance_taxa = lognormfiles_q2021[[2]][,-ncol(lognormfiles_q2021[[2]])]
lognorm_output_zellerfrance_age = lognormfiles_q2021[[2]][,ncol(lognormfiles_q2021[[2]])]

zellerfrance_intervals = c(10,20,30,40,50,60,70,80,90,100,
                           120,129)

lognorm_zellerfrance_powerMatrix_LM = matrix(nrow=length(zellerfrance_intervals), ncol = 3)

for(i in 1:length(zellerfrance_intervals)){
  lognorm_zellerfrance_powerMatrix_LM[i,] = unlist(repeatSubsampling_LM_function(maxiteration = 50,
                                                                                 subsamplesize = zellerfrance_intervals[i],
                                                                                 lognormFile = lognorm_output_zellerfrance_taxa,
                                                                                 metadata = lognorm_output_zellerfrance_age))
}

colnames(lognorm_zellerfrance_powerMatrix_LM) = c("Subsamplesize", "Mean", "SD")
lognorm_zellerfrance_powerMatrix_LM = as.data.frame(lognorm_zellerfrance_powerMatrix_LM)

powerPlot_zellerfrance_LM = ggplot(lognorm_zellerfrance_powerMatrix_LM, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = lognorm_zellerfrance_powerMatrix_LM, aes(x= Subsamplesize, y = Mean)) +
  geom_errorbar(aes(ymin = Mean - SD, 
                    ymax = Mean + SD)) +
  ggtitle("ZellerFrance (LM) n=129") +
  scale_x_continuous(trans = 'log10') + scale_y_continuous(limits= c(-6,35), breaks= seq(0, 260, by=5)) + 
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.title = element_text(size=11))

##____________________________________________________________MORGAN________________________
lognorm_output_morgan = lognormfiles_q2021[[1]]
lognorm_output_morgan_taxa = lognormfiles_q2021[[1]][,-ncol(lognormfiles_q2021[[1]])]
lognorm_output_morgan_age = lognormfiles_q2021[[1]][,ncol(lognormfiles_q2021[[1]])]

morgan_intervals = c(10,20,30,40,50,60,70,80,90,100,
                     120,140,160,180,200,228)

lognorm_morgan_powerMatrix_LM = matrix(nrow=length(morgan_intervals), ncol = 3)

for(i in 1:length(morgan_intervals)){
  lognorm_morgan_powerMatrix_LM[i,] = unlist(repeatSubsampling_LM_function(maxiteration = 50,
                                                                           subsamplesize = morgan_intervals[i],
                                                                           lognormFile = lognorm_output_morgan_taxa,
                                                                           metadata = lognorm_output_morgan_age))
}

colnames(lognorm_morgan_powerMatrix_LM) = c("Subsamplesize", "Mean", "SD")
lognorm_morgan_powerMatrix_LM = as.data.frame(lognorm_morgan_powerMatrix_LM)

powerPlot_morgan_LM = ggplot(lognorm_morgan_powerMatrix_LM, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = lognorm_morgan_powerMatrix_LM, aes(x= Subsamplesize, y = Mean)) +
  geom_errorbar(aes(ymin = Mean - SD, 
                    ymax = Mean + SD)) +
  ggtitle("Morgan (LM) n=228") +
  scale_x_continuous(trans = 'log10') + scale_y_continuous(limits= c(-6,35), breaks= seq(0, 260, by=5)) + 
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.title = element_text(size=11))


### combine all plots

fullsets_comparison_LM = ggarrange(powerPlot_nogbcn0_LM,powerPlot_nogstk_LM, powerPlot_escobar_LM,
   plot.new(), plot.new(), 
  powerPlot_zellerfrance_LM, powerPlot_zellergermany_LM,powerPlot_goodrich_LM, powerPlot_baxter_LM,powerPlot_ross_LM, 
  powerPlot_gloor_LM,powerPlot_morgan_LM, powerPlot_agp_LM, plot.new(), plot.new(), 
  labels = c(LETTERS[1:3], "", "", LETTERS[4:11]),
  nrow = 3, ncol = 5) 

annotate_figure(fullsets_comparison_LM, bottom="Sample size", 
                left="Number of Significant Taxa")


fullsets_comparison_LM_list = list(lognorm_morgan_powerMatrix_LM, lognorm_zellerfrance_powerMatrix_LM, lognorm_zellergermany_powerMatrix_LM,
                                   lognorm_baxter_powerMatrix_LM, lognorm_escobar_powerMatrix_LM, lognorm_gloor_powerMatrix_LM, 
                                   lognorm_goodrich_powerMatrix_LM, lognorm_nogbcn0_powerMatrix_LM, lognorm_nogstk_powerMatrix_LM,
                                   lognorm_ross_powerMatrix_LM, lognorm_agp_powerMatrix_LM)

for(i in 1:11){
  outfileName = paste0("/Users/anhil/Desktop/powerResults/", coh_l[[i]], "_powerResults_LM.txt")
  write.table(fullsets_comparison_LM_list[[i]], outfileName, quote = F, row.names = F)
}

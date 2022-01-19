## LM vs KCORR ggplots merged 


## agp 

test_lm = lognorm_agp_powerMatrix_LM
test_kcorr = lognorm_agp_powerMatrix_KCorr

test_lm$test = rep("Parametric", nrow(test_lm))
test_kcorr$test = rep("Non-Parametric", nrow(test_kcorr))

testpower_merge_agp = rbind(test_lm, test_kcorr)

compareStatModels_plot_agp = ggplot(testpower_merge_agp, aes(x=Subsamplesize, y=Mean)) +
  geom_point(data = testpower_merge_agp, aes(x= Subsamplesize, y = Mean, color=test)) +
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


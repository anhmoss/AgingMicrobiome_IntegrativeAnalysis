# Installing and loading packages for this study

packagesToInstall = c("ggplot2", "ggpubr", "vegan", "dplyr", "ggrepel")
install.packages(packagesToInstall) 
lapply(packagesToInstall, library, character.only=TRUE)

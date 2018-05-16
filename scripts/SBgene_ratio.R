#install and load relavant R packages
install.packages("faraway")
install.packages("ggplot2")
library(faraway)
library(ggplot2)

setwd("~/Tvedora_dev_RNAseq/")
datapath <- '~/Tvedora_dev_RNAseq/input/'
results <- '~/Tvedora_dev_RNAseq/output/'

#load the dataset
#compare the ratio of male biased gene among fold changes across stages
  sb_ratio <- read.table("~/Tvedora_dev_RNAseq/input/sb_ratio.txt", header = TRUE)
  str(sb_ratio)
  
  ggplot(sb_ratio, aes(x = foldc, y=ratio_mbias, group=stage), cex=2) + 
    geom_point(size = 5, aes(shape=factor(stage)))+ 
    scale_shape_manual(name="Gosner stage", values = c(18,15,17,16),labels=c("G27","G31","G43","G46")) +
    labs (x = "Fold change", y= "Proportion of male-biased genes", size =2) +
    theme(axis.text.x = element_text(size = 12,color = "black"))  + geom_line() +
    theme(axis.text.y = element_text(size = 12,color = "black")) +
    theme(text = element_text(size=14)) +
  theme(legend.text=element_text(size=11, color="black")) +
    labs(x = "Fold change", y = "Proportion male-biased genes") 
  

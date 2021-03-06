#install and load packages
install.packages("ggplot2")
install.packages("reshape2")

library(ggplot2)
library(reshape2)

#set up working directory
setwd <-'~/Tvedora_dev_RNAseq/'
datapath <- '~/Tvedora_dev_RNAseq/input/'
results <- '~/Tvedora_dev_RNAseq/output/'

#load datasets
sb <- read.table('~/Tvedora_dev_RNAseq/input/DE_gene_log.txt', header = TRUE)
str(sb)

sb <- data.frame(with(sb,sb[order(Gosner_stage,Fold_change,Sex_bias),]))
head(sb)

blues <- RColorBrewer::brewer.pal(4, "Blues")
reds <- RColorBrewer::brewer.pal(4, "Reds")

ggplot(data = sb, aes(x=Sex_bias, y= Gene_number, fill=interaction(factor(Fold_change),Sex_bias), label=Gene_number)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = c(reds,blues), name="Fold change",labels=c("FDR ≤0.05","Log2 ≥1","Log2 ≥2","Log2 ≥3","FDR ≤0.05","Log2 ≥1","Log2 ≥2","Log2 ≥3")) +
  facet_grid(~Gosner_stage) +
  geom_text(size = 4, position = position_stack(vjust = 0.5)) +
  scale_x_discrete(labels=c("XX", "XY'"),name="Sex bias") +
  scale_y_continuous(name = "Number of genes", limits = c(0,10500)) + 
  theme(axis.title.x = element_text(size=16,colour = "black"),axis.title.y = element_text(size=16,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12))




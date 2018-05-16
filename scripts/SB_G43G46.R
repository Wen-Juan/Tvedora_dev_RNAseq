#load dataset from output of edgeR_main.R.
violineData_de_G43 <- read.table("~/Tvedora_dev_RNAseq/input/violinData_de_G43.txt", header = TRUE)

ggplot(violineData_de_G43, aes(x=chr, y=logFC.XY43.XX43, fill=sexb)) + 
  scale_fill_manual(values = c("firebrick2","dodgerblue2"), name="Sex bias",labels=c("female","male")) +
  geom_rect(xmin = 1.5, xmax = 10.5, ymin=-Inf, ymax=Inf,alpha=0.5, fill="grey75") +
  geom_boxplot() +
  geom_vline(aes(xintercept=1.5), linetype="blank") +
  geom_vline(aes(xintercept=10.5), linetype="blank") +
  labs(x='Chromosome', y='Absolute Log2 values')

#G46
violineData_sort_g46 <- read.table("~/Tvedora_dev_RNAseq/input/violinData_de_G46.txt", header = TRUE)

ggplot(violineData_sort_g46, aes(x=chr, y=logFC.XY46.XX46, fill=sexb)) + 
  scale_fill_manual(values = c("firebrick2","dodgerblue2"), name="Sex bias",labels=c("female","male")) +
  geom_rect(xmin = 1.5, xmax = 10.5, ymin=-Inf, ymax=Inf,alpha=0.5, fill="grey75") +
  geom_boxplot() +
  geom_vline(aes(xintercept=1.5), linetype="blank") +
  geom_vline(aes(xintercept=10.5), linetype="blank") +
  labs(x='Chromosome', y='Absolute Log2 values')


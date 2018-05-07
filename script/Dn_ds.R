#install and load relevant R libraries and packages.
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("easyGgplot2", "kassambara")
library(easyGgplot2)

#set up the working directory
setwd("/Users/Wen-Juan/my_postdoc/useful_scripts/Tvedora_dev_RNAseq/input/")
#results directoy
setwd("/Users/Wen-Juan/my_postdoc/useful_scripts/Tvedora_dev_RNAseq/output/")

###OLD
dn_ds_all<-read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Tvedora_dev_RNAseq/input/dn_ds_all_stages.txt", header = T)
str(dn_ds_all)

###NEW
dn_ds_all_combine<-read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Tvedora_dev_RNAseq/input/dnds_allcombined_sb_un_g23g31g31g43g46.txt", header = T)
str(dn_ds_all_combine)

###OLD
share_dnds_4346 <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Tvedora_dev_RNAseq/input/share_sb_unbias_4346_absfi.txt", header = T)
str(share_dnds_4346)

###NEW
share_dnds_4346 <- read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Tvedora_dev_RNAseq/input/share_sb_unbias_4346_absfi_new.txt", header = T)
str(share_dnds_4346)

dn_ds <-read.table("/Users/Wen-Juan/my_postdoc/useful_scripts/Tvedora_dev_RNAseq/input/dn_ds.txt", header = T)
str(dn_ds)

dn_ds_auto <-read.table("/Users/Wen-Juan/git/rana_transcriptome/input/dnds/tv_all_dnds_new.txt", header = T)
str(dn_ds_auto)

#dn/ds plot between sex chromosome and autosomes combined.
ggplot(dn_ds, aes(x=reorder(type,dNdS), y=dNdS, fill=type)) + 
  scale_fill_manual(values = c("grey","firebrick2")) +
  theme(legend.position="none") +
  geom_boxplot() +
  labs(x='Chromosome', y='dn/ds') +
  scale_x_discrete(labels = c("Sex chromosome", "Autosome")) +
  theme(axis.text=element_text(size=12, color="black"),text = element_text(size=15,color="black"))

#boxplot for dn/ds of sex-biased genes along three developmental stages.
dn_ds_all <- dn_ds_all_combine
dn_ds_all <- subset (dn_ds_all,dn_ds_all$stage!='G31')
dn_ds_all$ID <- as.factor(dn_ds_all$ID)
head(dn_ds_all)

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Tvedora_dev_RNAseq/output/Figure_5_dNdS.pdf", width=8, height=8)
ggplot(dn_ds_all, aes(x=sb, y=dnds,fill=(sb))) + 
  geom_boxplot(notch = FALSE) +
  scale_fill_manual(values = c("firebrick2","dodgerblue2","grey")) +
  theme(legend.position="none") +
  facet_grid(~stage) +
  scale_x_discrete(labels=c("XX", "XY", "Unbias"),name="Sex bias") +
  scale_y_continuous(name = "dN/dS", limits = c(0,1.25)) + 
  theme(axis.title.x = element_text(size=14,colour = "black"),axis.title.y = element_text(size=14,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=9))
dev.off()



#combined all sex biased in three stages.
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Tvedora_dev_RNAseq/output/Figure_combined_dNdS.pdf", width=8, height=8)
ggplot(dn_ds_all_combine, aes(x=sb, y=dnds,fill=sb)) + 
  geom_boxplot(notch = FALSE) +
  ylim(0,1) +
  scale_fill_manual(values = c("firebrick2","dodgerblue2","grey")) +
  scale_x_discrete(labels=c("XX", "XY'", "Unbias"),name="Sex bias") +
  theme(legend.position="none")
dev.off()

wilcox.test(dn_ds_all$dnds[dn_ds_all$sb=='F-bias'], dn_ds_all$dnds[dn_ds_all$sb=='M-bias']) #W = 857470, p-value = 0.4227
wilcox.test(dn_ds_all$dnds[dn_ds_all$sb=='F-bias'], dn_ds_all$dnds[dn_ds_all$sb=='unbias']) #W = 2538500, p-value = 0.0935
wilcox.test(dn_ds_all$dnds[dn_ds_all$sb=='unbias'], dn_ds_all$dnds[dn_ds_all$sb=='M-bias']) #W = 1656800, p-value = 0.5779

#two stages for dN
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Tvedora_dev_RNAseq/output/Figure_S7_dN.pdf", width=8, height=8)
dn_ds_all_combine_sub <- subset(dn_ds_all_combine,dn_ds_all_combine$stage!='G31')
ggplot(dn_ds_all_combine_sub, aes(x=sb, y=dn,fill=(sb))) + 
  geom_boxplot(notch = FALSE) +
  scale_fill_manual(values = c("firebrick2","dodgerblue2","grey")) +
  theme(legend.position="none") +
  facet_grid(~stage) +
  scale_x_discrete(labels=c("XX", "XY", "Unbias"),name="Sex bias") +
  scale_y_continuous(name = "dN/dS", limits = c(0,2.1)) + 
  theme(axis.title.x = element_text(size=14,colour = "black"),axis.title.y = element_text(size=14,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=9))
dev.off()

#two stages for dS
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Tvedora_dev_RNAseq/output/Figure_S7_dS.pdf", width=8, height=8)
ggplot(dn_ds_all_combine_sub, aes(x=sb, y=ds,fill=(sb))) + 
  geom_boxplot(notch = FALSE) +
  scale_fill_manual(values = c("firebrick2","dodgerblue2","grey")) +
  theme(legend.position="none") +
  facet_grid(~stage) +
  scale_x_discrete(labels=c("XX", "XY", "Unbias"),name="Sex bias") +
  scale_y_continuous(name = "dN/dS", limits = c(0,2.3)) + 
  theme(axis.title.x = element_text(size=14,colour = "black"),axis.title.y = element_text(size=14,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=9))
dev.off()

#new
dn_ds_all_combine_sub2 <- subset (dn_ds_all_combine, dn_ds_all_combine$stage=='G46')
str(dn_ds_all_combine_sub2)
wilcox.test(dn_ds_all_combine_sub2$dn[dn_ds_all_combine_sub2$sb=='unbias'], dn_ds_all_combine_sub2$dn[dn_ds_all_combine_sub2$sb=='M-bias']) #W = 852650, p-value = 0.01714
wilcox.test(dn_ds_all_combine_sub2$ds[dn_ds_all_combine_sub2$sb=='unbias'], dn_ds_all_combine_sub2$ds[dn_ds_all_combine_sub2$sb=='M-bias']) #W = 900740, p-value = 7.28e-07

#dn/ds of shared sex-biased genes and unbiased between G43 and G46 stages. 
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Tvedora_dev_RNAseq/output/Figure6_sharedsb_unb.pdf", width=8, height=8)
ggplot(share_dnds_4346, aes(x=sb, y=dnds,fill=sb)) + 
  geom_boxplot(notch = FALSE) + 
  labs(x="Sex bias", y="dN/dS", aex.lab=2, aex.axis=1.5, aex.main=2) +
  ylim(0,0.9) +
  scale_fill_manual(values = c("firebrick2","dodgerblue2","grey")) +
  scale_x_discrete(labels=c("XX", "XY'", "Unbias"),name="Sex bias") +
  theme(legend.position="none")
dev.off()

#old  
wilcox.test(share_dnds_4346$dnds[share_dnds_4346$sb=='unbias'], share_dnds_4346$dnds[share_dnds_4346$sb=='female']) #W = 18554, p-value = 2.644e-05
wilcox.test(share_dnds_4346$dnds[share_dnds_4346$sb=='unbias'], share_dnds_4346$dnds[share_dnds_4346$sb=='male']) #W = 2659, p-value = 0.001536
wilcox.test(share_dnds_4346$dnds[share_dnds_4346$sb=='female'], share_dnds_4346$dnds[share_dnds_4346$sb=='male']) #W = 180, p-value = 0.2896

#new
wilcox.test(share_dnds_4346$dnds[share_dnds_4346$sb=='unbias'], share_dnds_4346$dnds[share_dnds_4346$sb=='female']) #W = 23424, p-value = 5.185e-05
wilcox.test(share_dnds_4346$dnds[share_dnds_4346$sb=='unbias'], share_dnds_4346$dnds[share_dnds_4346$sb=='male']) #W = 3403.5, p-value = 0.002034
wilcox.test(share_dnds_4346$dnds[share_dnds_4346$sb=='female'], share_dnds_4346$dnds[share_dnds_4346$sb=='male']) #W = 180, p-value = 0.2896

#scatter plot with outlier
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Tvedora_dev_RNAseq/output/Figure7a_g43_with.pdf", width=8, height=8)
ggplot2.scatterplot(data=share_dnds_4346, xName='abs43_XY_XX',yName='dnds', 
                    size=2,
                    addRegLine=TRUE, addConfidenceInterval=TRUE)  +
  labs(x="Log2 (XY/XX) ratio of gene expression at G43", y="dn/ds", color ="Sex bias")
dev.off()

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Tvedora_dev_RNAseq/output/Figure7a_g46_with.pdf", width=8, height=8)
ggplot2.scatterplot(data=share_dnds_4346, xName='abs46_XY_XX',yName='dnds', 
                    size=2,
                    addRegLine=TRUE, addConfidenceInterval=TRUE)  +
  labs(x="Log2 (XY/XX) ratio of gene expression at G46", y="dn/ds", color ="Sex bias")
dev.off()

###if removing the outliner
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Tvedora_dev_RNAseq/output/Figure7a_g43_rm.pdf", width=8, height=8)
ggplot2.scatterplot(data=share_dnds_4346rm, xName='abs43_XY_XX',yName='dnds', 
                    size=2,
                    addRegLine=TRUE, addConfidenceInterval=TRUE)  +
  labs(x="Log2 (XY/XX) ratio of gene expression at G43", y="dn/ds", color ="Sex bias")
dev.off()

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Tvedora_dev_RNAseq/output/Figure7a_g46_rm.pdf", width=8, height=8)
ggplot2.scatterplot(data=share_dnds_4346rm, xName='abs46_XY_XX',yName='dnds', 
                    size=2,
                    addRegLine=TRUE, addConfidenceInterval=TRUE)  +
  labs(x="Log2 (XY/XX) ratio of gene expression at G46", y="dn/ds", color ="Sex bias")
dev.off()

#subset only sex-bised genes
share_dnds_4346_sub <- subset(share_dnds_4346rm,share_dnds_4346rm$sb!='unbias')
pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Tvedora_dev_RNAseq/output/Figure7a_g43_sbonly.pdf", width=8, height=8)
ggplot2.scatterplot(data=share_dnds_4346_sub, xName='abs43_XY_XX',yName='dnds', 
                    size=2,
                    addRegLine=TRUE, addConfidenceInterval=TRUE)  +
  labs(x="Log2 (XY/XX) ratio of gene expression at G43", y="dn/ds", color ="Sex bias")
dev.off()

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Tvedora_dev_RNAseq/output/Figure7a_g46_sbonly.pdf", width=8, height=8)
ggplot2.scatterplot(data=share_dnds_4346_sub, xName='abs46_XY_XX',yName='dnds', 
                    size=2,
                    addRegLine=TRUE, addConfidenceInterval=TRUE)  +
  labs(x="Log2 (XY/XX) ratio of gene expression at G46", y="dn/ds", color ="Sex bias")
dev.off()

shapiro.test(sqrt(share_dnds_4346rm$abs43_XY_XX)) #non normal distribution

#with the outlier
summary(glm(sqrt(dnds)~abs43_XY_XX, family = binomial, data=share_dnds_4346)) #not significant P=0.6
summary(glm(sqrt(dnds)~abs43_XY_XX*sb, family = binomial, data=share_dnds_4346)) #not significant
summary(glm(sqrt(dnds)~abs46_XY_XX, family = binomial, data=share_dnds_4346)) #not significant
summary(glm(sqrt(dnds)~abs46_XY_XX*sb, family = binomial, data=share_dnds_4346))#not significant

#without outlier
summary(glm(sqrt(dnds)~abs43_XY_XX, family = binomial, data=share_dnds_4346rm)) #not significant P=0.9
summary(glm(sqrt(dnds)~abs43_XY_XX*sb, family = binomial, data=share_dnds_4346rm)) #not significant
summary(glm(sqrt(dnds)~abs46_XY_XX, family = binomial, data=share_dnds_4346rm)) #not significant
summary(glm(sqrt(dnds)~abs46_XY_XX*sb, family = binomial, data=share_dnds_4346rm))#not significant

#with only sb genes
summary(glm(sqrt(dnds)~abs43_XY_XX, family = binomial, data=share_dnds_4346_sub)) #not significant 
summary(glm(sqrt(dnds)~abs43_XY_XX*sb, family = binomial, data=share_dnds_4346_sub)) #not significant
summary(glm(sqrt(dnds)~abs46_XY_XX, family = binomial, data=share_dnds_4346_sub)) #not significant
summary(glm(sqrt(dnds)~abs46_XY_XX*sb, family = binomial, data=share_dnds_4346_sub))#not significant

#summary(glm(dnds~abs43_XY_XX*sb, family = binomial, data=share_dnds_4346))

pdf("/Users/Wen-Juan/my_postdoc/useful_scripts/Tvedora_dev_RNAseq/output/Figure7b_g46.pdf", width=8, height=8)
ggplot2.scatterplot(data=share_dnds_4346, xName='abs46_XY_XX',yName='dnds', 
                    size=2,
                    addRegLine=TRUE, addConfidenceInterval=TRUE)  +
  labs(x="Log2 (XY/XX) ratio of gene expression at G46", y="dn/ds", color ="Sex bias")
dev.off()

#dn/ds for all genes on sex chromsome vs autosomes combined.
ggplot(dn_ds, aes(x=type, y=dNdS,fill=type)) + 
  geom_boxplot(notch = FALSE) + 
  labs(x="Sex bias", y="dN/dS", aex.lab=2, aex.axis=1.5, aex.main=2) +
  scale_fill_manual(values = c("grey","firebrick2")) +
  scale_x_discrete(labels=c("Autosome", "Sex chromosome"),name="Chromosome") +
  theme(legend.position="none")

#G43, unbias vs. F-bias
dnds_g43 <- subset(dn_ds_all_combine, dn_ds_all_combine$stage=='G43')

#NEW g43
wilcox.test(dnds_g43$dnds[dnds_g43$ID=='G43_unbias'], dnds_g43$dnds[dnds_g43$ID=='43_F-bias'], exact=FALSE) #W = 25598, p-value = 1.253e-05
wilcox.test(dnds_g43$dnds[dnds_g43$ID=='G43_unbias'], dnds_g43$dnds[dnds_g43$ID=='43_M-bias'], exact=FALSE) #W = 4710.5, p-value = 0.000235
wilcox.test(dnds_g43$dnds[dnds_g43$ID=='43_M-bias'], dnds_g43$dnds[dnds_g43$ID=='43_F-bias'], exact=FALSE) #W = 434, p-value = 0.2222

#new G46
dnds_g46 <- subset(dn_ds_all_combine, dn_ds_all_combine$stage=='G46')
wilcox.test(dnds_g46$dnds[dnds_g46$ID=='G46_unbias'], dnds_g46$dnds[dnds_g46$ID=='46_F-bias'], exact=FALSE) #W = 1320400, p-value = 0.04221
wilcox.test(dnds_g46$dnds[dnds_g46$ID=='G46_unbias'], dnds_g46$dnds[dnds_g46$ID=='46_M-bias'], exact=FALSE) #W = 823710, p-value = 0.4029
wilcox.test(dnds_g46$dnds[dnds_g46$ID=='46_F-bias'], dnds_g46$dnds[dnds_g46$ID=='46_M-bias'], exact=FALSE) #W = 816360, p-value = 0.3217

#new combined
wilcox.test(dn_ds_all_combine$dnds[dn_ds_all_combine$sb=='M-bias'], dn_ds_all_combine$dnds[dn_ds_all_combine$sb=='F-bias'], exact=FALSE) #W = 890990, p-value = 0.4001
wilcox.test(dn_ds_all_combine$dnds[dn_ds_all_combine$sb=='unbias'], dn_ds_all_combine$dnds[dn_ds_all_combine$sb=='F-bias'], exact=FALSE) #W = 2692000, p-value = 0.0935
wilcox.test(dn_ds_all_combine$dnds[dn_ds_all_combine$sb=='unbias'], dn_ds_all_combine$dnds[dn_ds_all_combine$sb=='M-bias'], exact=FALSE) #W = 1656900, p-value = 0.6084

########old
#Wilcoxon rank sum test with continuity correction

#data:  dn_ds_all$dnds[dn_ds_all$ID == "43_unbias"] and dn_ds_all$dnds[dn_ds_all$ID == "43_F-bias"]
#W = 65316, p-value = 2.357e-06
#alternative hypothesis: true location shift is not equal to 0

#########
#G43, unbias vs. M-bias
#wilcox.test(dn_ds_all$dnds[dn_ds_all$ID=='43_unbias'], dn_ds_all$dnds[dn_ds_all$ID=='43_M-bias'])

#########
#Wilcoxon rank sum test with continuity correction
#data:  dn_ds_all$dnds[dn_ds_all$ID == "43_unbias"] and dn_ds_all$dnds[dn_ds_all$ID == "43_M-bias"]
#W = 11844, p-value = 0.0001248
#alternative hypothesis: true location shift is not equal to 0
#########

#G43, F-bias vs. M-bias
wilcox.test(dn_ds_all$dnds[dn_ds_all$ID=='43_F-bias'], dn_ds_all$dnds[dn_ds_all$ID=='43_M-bias'])

########
##Wilcoxon rank sum test with continuity correction
##data:  dn_ds_all$dnds[dn_ds_all$ID == "43_F-bias"] and dn_ds_all$dnds[dn_ds_all$ID == "43_M-bias"]
#W = 280, p-value = 0.2222
#alternative hypothesis: true location shift is not equal to 0
##############

#G46, F-bias vs. M-bias
wilcox.test(dn_ds_all$dnds[dn_ds_all$ID=='46_F-bias'], dn_ds_all$dnds[dn_ds_all$ID=='46_M-bias'])

########
#Wilcoxon rank sum test with continuity correction
#data:  dn_ds_all$dnds[dn_ds_all$ID == "46_F-bias"] and dn_ds_all$dnds[dn_ds_all$ID == "46_M-bias"]
#W = 816360, p-value = 0.3217
#alternative hypothesis: true location shift is not equal to 0
########

#G46, unbias vs. F-bias
wilcox.test(dn_ds_all$dnds[dn_ds_all$ID=='46_unbias'], dn_ds_all$dnds[dn_ds_all$ID=='46_F-bias'])

#########
##Wilcoxon rank sum test with continuity correction
#data:  dn_ds_all$dnds[dn_ds_all$ID == "46_unbias"] and dn_ds_all$dnds[dn_ds_all$ID == "46_F-bias"]
#W = 1327900, p-value = 0.03341
#alternative hypothesis: true location shift is not equal to 0
#########

#G46, unbias vs. M-bias
wilcox.test(dn_ds_all$dnds[dn_ds_all$ID=='46_unbias'], dn_ds_all$dnds[dn_ds_all$ID=='46_M-bias'])

#########
#Wilcoxon rank sum test with continuity correction
#data:  dn_ds_all$dnds[dn_ds_all$ID == "46_unbias"] and dn_ds_all$dnds[dn_ds_all$ID == "46_M-bias"]
#W = 828440, p-value = 0.3571
#ralternative hypothesis: true location shift is not equal to 0
#########

wilcox.test(dn_ds$dNdS[dn_ds$type=='auto'], dn_ds$dNdS[dn_ds$type=='sexchr']) #W = 4191400, p-value = 0.2871


#### chr02 vs auto
pdf("/Users/Wen-Juan/git/rana_transcriptome/output/dnds_tv_chr01chr02_auto.pdf", width=8, height=8)
ggplot(dn_ds_auto, aes(x=chr, y=dNdS,fill=(chr))) + 
  geom_boxplot(notch = TRUE) +
  scale_fill_manual(values = c("grey","firebrick2","firebrick2")) +
  theme(legend.position="none") +
  scale_y_continuous(name = "dN/dS", limits = c(0,1.25)) + 
  theme(axis.title.x = element_text(size=14,colour = "black"),axis.title.y = element_text(size=14,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=9))
dev.off()
wilcox.test(dn_ds_auto$dNdS[dn_ds_auto$chr=='Chr02'], dn_ds_auto$dNdS[dn_ds_auto$chr=='auto']) #W = 1029800, p-value = 0.2513
wilcox.test(dn_ds_auto$dNdS[dn_ds_auto$chr=='Chr01'], dn_ds_auto$dNdS[dn_ds_auto$chr=='auto']) #W = 1257600, p-value = 0.5552
wilcox.test(dn_ds_auto$dNdS[dn_ds_auto$chr=='Chr01'], dn_ds_auto$dNdS[dn_ds_auto$chr=='Chr02']) #W = 252580, p-value = 0.1535

#install and load relevant R libraries and packages.
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("easyGgplot2", "kassambara")
library(easyGgplot2)

#set up the working directory
setwd("~/Tvedora_dev_RNAseq/input/")
#results directoy
setwd("~/Tvedora_dev_RNAseq/output/")

dn_ds_all_combine<-read.table("~/Tvedora_dev_RNAseq/input/dnds_allcombined_sb_un_g23g31g31g43g46.txt", header = T)
str(dn_ds_all_combine)

share_dnds_4346 <- read.table("~/Tvedora_dev_RNAseq/input/share_sb_unbias_4346_absfi_new.txt", header = T)
str(share_dnds_4346)

dn_ds <-read.table("~/Tvedora_dev_RNAseq/input/dn_ds.txt", header = T)
str(dn_ds)

dn_ds_auto <-read.table("~/Tvedora_dev_RNAseq/input/dnds/tv_all_dnds_new.txt", header = T)
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

ggplot(dn_ds_all, aes(x=sb, y=dnds,fill=(sb))) + 
  geom_boxplot(notch = FALSE) +
  scale_fill_manual(values = c("firebrick2","dodgerblue2","grey")) +
  theme(legend.position="none") +
  facet_grid(~stage) +
  scale_x_discrete(labels=c("XX", "XY", "Unbias"),name="Sex bias") +
  scale_y_continuous(name = "dN/dS", limits = c(0,1.25)) + 
  theme(axis.title.x = element_text(size=14,colour = "black"),axis.title.y = element_text(size=14,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=9))


#combined all sex biased in three stages.
ggplot(dn_ds_all_combine, aes(x=sb, y=dnds,fill=sb)) + 
  geom_boxplot(notch = FALSE) +
  ylim(0,1) +
  scale_fill_manual(values = c("firebrick2","dodgerblue2","grey")) +
  scale_x_discrete(labels=c("XX", "XY'", "Unbias"),name="Sex bias") +
  theme(legend.position="none")

wilcox.test(dn_ds_all$dnds[dn_ds_all$sb=='F-bias'], dn_ds_all$dnds[dn_ds_all$sb=='M-bias']) #W = 857470, p-value = 0.4227
wilcox.test(dn_ds_all$dnds[dn_ds_all$sb=='F-bias'], dn_ds_all$dnds[dn_ds_all$sb=='unbias']) #W = 2538500, p-value = 0.0935
wilcox.test(dn_ds_all$dnds[dn_ds_all$sb=='unbias'], dn_ds_all$dnds[dn_ds_all$sb=='M-bias']) #W = 1656800, p-value = 0.5779

##dnds for female-biased genes resulting from the comparison between XY0 females and XY0 males at stage G43.
XY_sbun_dnds <- read.table("~/Tvedora_dev_RNAseq/input/XY_unsb_dnds.txt", header = T)
str(XY_sbun_dnds)

pdf("~/Tvedora_dev_RNAseq/output/XY_sbun_dnds.pdf", width=8, height=8)
ggplot(XY_sbun_dnds, aes(x=bias, y=dnds, fill=bias)) + 
  scale_fill_manual(values = c("firebrick2","grey")) +
  theme(legend.position="none") +
  ylim(0,0.9) +
  geom_boxplot() +
  labs(x='Sex bias', y='dN/dS') +
  theme(axis.text=element_text(size=12, color="black"),text = element_text(size=15,color="black"))
dev.off()

wilcox.test(XY_sbun_dnds$dnds[XY_sbun_dnds$bias=='unbias'], XY_sbun_dnds$dnds[XY_sbun_dnds$bias=='female']) #W = 7805, p-value = 0.003315

##dnds for sex-biased genes resulting from the comparison between XX male and XX females at G46.
XX_sbun_dnds <- read.table("~/Tvedora_dev_RNAseq/input/g46_sbun_dnds.txt", header = T)
str(XX_sbun_dnds)

pdf("~/Tvedora_dev_RNAseq/output/XXg46_sbun_dnds.pdf",width=8, height=8)
ggplot(XX_sbun_dnds, aes(x=bias, y=dnds, fill=bias)) + 
  scale_fill_manual(values = c("firebrick2","dodgerblue2","grey")) +
  theme(legend.position="none") +
  ylim(0,0.9) +
  geom_boxplot() +
  labs(x='Sex bias', y='dN/dS') +
  theme(axis.text=element_text(size=12, color="black"),text = element_text(size=15,color="black"))
dev.off()

wilcox.test(XX_sbun_dnds$dnds[XX_sbun_dnds$bias=='unbias'], XX_sbun_dnds$dnds[XX_sbun_dnds$bias=='female']) #W = 1166900, p-value = 0.3559
wilcox.test(XX_sbun_dnds$dnds[XX_sbun_dnds$bias=='unbias'], XX_sbun_dnds$dnds[XX_sbun_dnds$bias=='male']) #W = 182880, p-value = 0.1689
wilcox.test(XX_sbun_dnds$dnds[XX_sbun_dnds$bias=='male'], XX_sbun_dnds$dnds[XX_sbun_dnds$bias=='female']) #W = 193950, p-value = 0.06187

#two stages for dN
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

#two stages for dS
ggplot(dn_ds_all_combine_sub, aes(x=sb, y=ds,fill=(sb))) + 
  geom_boxplot(notch = FALSE) +
  scale_fill_manual(values = c("firebrick2","dodgerblue2","grey")) +
  theme(legend.position="none") +
  facet_grid(~stage) +
  scale_x_discrete(labels=c("XX", "XY", "Unbias"),name="Sex bias") +
  scale_y_continuous(name = "dN/dS", limits = c(0,2.3)) + 
  theme(axis.title.x = element_text(size=14,colour = "black"),axis.title.y = element_text(size=14,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=9))

#new
dn_ds_all_combine_sub2 <- subset (dn_ds_all_combine, dn_ds_all_combine$stage=='G46')
str(dn_ds_all_combine_sub2)
wilcox.test(dn_ds_all_combine_sub2$dn[dn_ds_all_combine_sub2$sb=='unbias'], dn_ds_all_combine_sub2$dn[dn_ds_all_combine_sub2$sb=='M-bias']) #W = 852650, p-value = 0.01714
wilcox.test(dn_ds_all_combine_sub2$ds[dn_ds_all_combine_sub2$sb=='unbias'], dn_ds_all_combine_sub2$ds[dn_ds_all_combine_sub2$sb=='M-bias']) #W = 900740, p-value = 7.28e-07

#dn/ds of shared sex-biased genes and unbiased between G43 and G46 stages. 
ggplot(share_dnds_4346, aes(x=sb, y=dnds,fill=sb)) + 
  geom_boxplot(notch = FALSE) + 
  labs(x="Sex bias", y="dN/dS", aex.lab=2, aex.axis=1.5, aex.main=2) +
  ylim(0,0.9) +
  scale_fill_manual(values = c("firebrick2","dodgerblue2","grey")) +
  scale_x_discrete(labels=c("XX", "XY'", "Unbias"),name="Sex bias") +
  theme(legend.position="none")

#new
wilcox.test(share_dnds_4346$dnds[share_dnds_4346$sb=='unbias'], share_dnds_4346$dnds[share_dnds_4346$sb=='female']) #W = 23424, p-value = 5.185e-05
wilcox.test(share_dnds_4346$dnds[share_dnds_4346$sb=='unbias'], share_dnds_4346$dnds[share_dnds_4346$sb=='male']) #W = 3403.5, p-value = 0.002034
wilcox.test(share_dnds_4346$dnds[share_dnds_4346$sb=='female'], share_dnds_4346$dnds[share_dnds_4346$sb=='male']) #W = 180, p-value = 0.2896


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

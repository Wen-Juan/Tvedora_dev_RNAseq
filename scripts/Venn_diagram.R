### venns in R
install.packages("VennDiagram")

library("VennDiagram")
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

#for the shared female-biased genes across five stages.
setwd <-'~/Tvedora_dev_RNAseq/'
datapath <- '~/Tvedora_dev_RNAseq/input/sex_bias/'
results <- '~/Tvedora_dev_RNAseq/input/sex_bias/'

fivestages_fbias <- read.table("~/Tvedora_dev_RNAseq/input/sexbias_female_inclusr.txt", header = TRUE)
str(fivestages_fbias)

g43XY_fbias <- read.table("~/Tvedora_dev_RNAseq/input/XY_fbias.txt", header = TRUE)

fbis_g23 <- fivestages_fbias$trans[fivestages_fbias$stage=='G23']
fbis_g27 <- fivestages_fbias$trans[fivestages_fbias$stage=='G27']
fbis_g31 <- fivestages_fbias$trans[fivestages_fbias$stage=='G31']
fbis_g43 <- fivestages_fbias$trans[fivestages_fbias$stage=='G43']
fbis_g46 <- fivestages_fbias$trans[fivestages_fbias$stage=='G46']
fbis_g43xy <- fivestages_fbias$trans[fivestages_fbias$stage=='G43XY']
fbis_g46xx <- fivestages_fbias$trans[fivestages_fbias$stage=='G46XX']

#for shared female-biased genes among five stages
venn.plot <- venn.diagram(list(G23 = as.character(fbis_g23), G27 = as.character(fbis_g27), G31 = as.character(fbis_g31), G43=as.character(fbis_g43), G46=as.character(fbis_g46)), filename =NULL,
                          fill=c("orange","red","black","green","blue"),
                          ext.line.lwd = 3,
                          cex = 2,
                          cat.cex = 2.5,
                          rotation.degree = 60)

# to draw to the screen:
grid.arrange(gTree(children=venn.plot),ncol = 1 )

# to output to pdf
venn_fbias_out <- arrangeGrob(gTree(children=venn.plot),ncol = 1 )
ggsave(file="shared_fembias.pdf", venn_fbias_out, path = "~/Tvedora_dev_RNAseq/output/")

#overlap of G43 female-biased genes resulting from comparison between XX females, or XY0 females and XY0 males.
venn.plot1 <- venn.diagram(list(XX_f=as.character(fbis_g43), XY0_f=as.character(fbis_g43xy)), filename =NULL,
                          fill=c("orange","red"),
                          ext.line.lwd = 2,
                          cex = 1,
                          cat.cex = 1,rotation.degree = 60)

grid.arrange(gTree(children=venn.plot1),ncol = 1 )
venn_fbias_out <- arrangeGrob(gTree(children=venn.plot1),ncol = 1 )
ggsave(file="g43_XYXXfemale_fembias.pdf", venn_fbias_out, path = "~/Tvedora_dev_RNAseq/")

##overlap of G46 female-biased genes resulting from comparison between XX females and either XY0 males or XX male.
venn.plot2 <- venn.diagram(list(XY0_m=as.character(fbis_g46), XX_m=as.character(fbis_g46xx)), filename =NULL,
                           fill=c("orange","red"),
                           ext.line.lwd = 2,
                           cex = 1,
                           cat.cex = 1,rotation.degree = 60)

grid.arrange(gTree(children=venn.plot2),ncol = 1 )
venn_fbias_out <- arrangeGrob(gTree(children=venn.plot2),ncol = 1 )
ggsave(file="g46_XXmale_fembias.pdf", venn_fbias_out, path = "~/Tvedora_dev_RNAseq/output/")


#test whether the overlap of female-biased genes across five stages is greater than by chance.
f_bias <- list(as.character(fbis_g23), as.character(fbis_g27), as.character(fbis_g31),as.character(fbis_g43),as.character(fbis_g46))
list(f_bias)
str(f_bias)

total1 <- 9680

####### TO DO for all interections

res=supertest(f_bias, n=total1)
plot(res, sort.by="size")
plot(res, Layout="landscape", degree=2:4, sort.by="size")
summary(res)

write.csv(summary(res)$Table, file="~/Tvedora_dev_RNAseq/output/f_bias_testbychance.csv", row.names=FALSE)


#for shared male-biased genes
mbias <- read.table("~/Tvedora_dev_RNAseq/input/share_mbias_incluesr.txt", header = TRUE)
str(mbias)


mbis_g27 <- mbias$trans[mbias$stage=='G27']
mbis_g31 <- mbias$trans[mbias$stage=='G31']
mbis_g43 <- mbias$trans[mbias$stage=='G43']
mbis_g46 <- mbias$trans[mbias$stage=='G46']
mbis_g46xx <- mbias$trans[mbias$stage=='G46XX']


venn.plot.2 <- venn.diagram(list(G27 = as.character(mbis_g27), G31 = as.character(mbis_g31), G43=as.character(mbis_g43), G46=as.character(mbis_g46)), filename =NULL,
                          fill=c("green","blue"),
                          ext.line.lwd = 3,
                          cex = 2,
                          cat.cex = 2.5)

## to draw to the screen:
grid.arrange(gTree(children=venn.plot.2),ncol = 1 )

# to output to pdf
venn_mbias_out <- arrangeGrob(gTree(children=venn.plot.2),ncol = 1 )
ggsave(file="shared_mbias.pdf", venn_mbias_out, path = "~/Tvedora_dev_RNAseq/output/")

#overlap of male-biased genes resulting from the comparison between XX females and either XY0 males, or XX male.
venn.plot.3 <- venn.diagram(list(G46=as.character(mbis_g46), G46XX=as.character(mbis_g46xx)), filename =NULL,
                            fill=c("green","blue"),
                            ext.line.lwd = 2,
                            cex = 1,
                            cat.cex = 1,rotation.degree = 60)

grid.arrange(gTree(children=venn.plot.3),ncol = 1 )
venn_mbias_out <- arrangeGrob(gTree(children=venn.plot.3),ncol = 1 )
ggsave(file="g46_xxmale_mbias.pdf", venn_mbias_out, path = "~/Tvedora_dev_RNAseq/output/")


#test whether the overlap of male-biased genes across five stages is greater than by chance.
m_bias <- list(as.character(mbis_g27), as.character(mbis_g31),as.character(mbis_g43),as.character(mbis_g46))
list(m_bias)
str(m_bias)

total2 <- 6217

####### TO DO for all interections

res=supertest(m_bias, n=total2)
plot(res, Layout="landscape", degree=2:4, sort.by="size")
summary(res)

write.csv(summary(res)$Table, file="~/Tvedora_dev_RNAseq/output/mbias_testbychance.csv", row.names=FALSE)


#for shared unbiased genes
fivestages_unbias <- read.table("~/Tvedora_dev_RNAseq/input/share_unbias.txt", header = TRUE)
str(fivestages_unbias)

unbias_g23 <- fivestages_unbias$unbiased[fivestages_unbias$stage=='G23']
unbias_g27 <- fivestages_unbias$unbiased[fivestages_unbias$stage=='G27']
unbias_g31 <- fivestages_unbias$unbiased[fivestages_unbias$stage=='G31']
unbias_g43 <- fivestages_unbias$unbiased[fivestages_unbias$stage=='G43']
unbias_g46 <- fivestages_unbias$unbiased[fivestages_unbias$stage=='G46']

venn.plot <- venn.diagram(list(G23 = as.character(unbias_g23), G27 = as.character(unbias_g27), G31 = as.character(unbias_g31), G43=as.character(unbias_g43), G46=as.character(unbias_g46)), filename =NULL,
                          cat.col=c("orange","red","black","green","blue"),
                          fill=c("orange","red","black","green","blue"),
                          ext.line.lwd = 3,
                          cex = 1.5,
                          cat.cex = 2.5,
                          rotation.degree = 60)

grid.arrange(gTree(children=venn.plot),ncol = 1 )

venn_unbias_out <- arrangeGrob(gTree(children=venn.plot),ncol = 1 )
ggsave(file="shared_unbias.pdf", venn_unbias_out, path = "~/Tvedora_dev_RNAseq/output/")


#for shared sex-biased and unbiased genes between G43 and G46
g4346_fbias <- read.table("~/Tvedora_dev_RNAseq/input/G4346_unbias_sb.txt", header = TRUE)
str(g4346_fbias)

sb_g43 <- g4346_fbias$trans[g4346_fbias$stage=='G43_sb']
unbis_g43 <- g4346_fbias$trans[g4346_fbias$stage=='G43_unbias']
sb_g46 <- g4346_fbias$trans[g4346_fbias$stage=='G46_sb']
unbis_g46 <- g4346_fbias$trans[g4346_fbias$stage=='G46_unbias']

venn.plot.1 <- venn.diagram(list(G43_SB = as.character(sb_g43),G43_UN = as.character(unbis_g43),G46_SB=as.character(sb_g46),G46_UN=as.character(unbis_g46)), filename =NULL,
                          fill=c("red","green","orange","blue"),
                          ext.line.lwd = 3,
                          cex = 1.5,
                          cat.cex = 2)

# to draw to the screen:
grid.arrange(gTree(children=venn.plot.1),ncol = 1 )
# to output to pdf
venn_shareg4346_out <- arrangeGrob(gTree(children=venn.plot.1),ncol = 1 )
ggsave(file="G4346_sharedunbias.pdf", venn_shareg4346_out, path = "~/Tvedora_dev_RNAseq/output/")

#test whether the overlap of sex-biased genes between G43 and G46 is greater than by chance.
share_G43G46 <- list(as.character(sb_g43), as.character(unbis_g43),as.character(sb_g46),as.character(unbis_g46))
list(share_G43G46)
str(share_G43G46)

total3 <- 42410

####### TO DO for all interections

res=supertest(share_G43G46, n=total3)
plot(res, Layout="landscape", degree=2:4, sort.by="size")
summary(res)

write.csv(summary(res)$Table, file="~/Tvedora_dev_RNAseq/output/g43g46_sbunbias_testbychance.csv", row.names=FALSE)

#correction for multiple test.
Pvals = c(0.000013, 0.0002, 0.04,0.4) ### vector of pvals
p.adjust(Pvals, method = "BH") 





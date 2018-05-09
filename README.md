# Tvedora_dev_RNAseq 
In this project, we investigate the evolution and dynamics of sex-biased genes and gene expression along the five developmental stages of the common frog Rana temporaria, in a population with proto-Y chromosomes.


## Instruction of this repository

input folder is for relevant input files.
output folder is for some output files generated by according scripts.
scripts folder is for relevant scripts used in this projects.


### Step1. Quality check and trim of raw reads of all RNAseq libraries.

~/Tvedora_dev_RNAseq/scripts/Trim_reads_qualitycheck.sh


### Step2. De novo transcriptome assembly and a series of filtering steps. 

~/Tvedora_dev_RNAseq/scripts/De_novo_assembly.sh
~/Tvedora_dev_RNAseq/scripts/Transcriptome_filtering_pipeline.sh


### Step3. Transcriptome annotation and Gene Ontology annotation.

~/Tvedora_dev_RNAseq/scripts/Tv_annotation.sh
~/Tvedora_dev_RNAseq/scripts/Remove_short_transcripts.pl
~/Tvedora_dev_RNAseq/scripts/Blast_RBH.py


### Step4. Quantify the abundance of transcripts in each samples with Kallito.

1) index the transcriptome
~/Tvedora_dev_RNAseq/scripts/Kallisto_index.sh

2) quantify the abundance of transcripts 
~/Tvedora_dev_RNAseq/scripts/Kallisto_abundance.sh

3) construct the matrix of read counts 
~/Tvedora_dev_RNAseq/scripts/Kallisto_matrix.sh


### Step5. Differential gene expression analysis with EdgeR (R package).

~/Tvedora_dev_RNAseq/scripts/EdgeR_main.r

#Instruction for use of EdgR_main.r is as the following:
#Installation
#Install the following libraries are available

install.packages("gplots", dependencies=TRUE)
install.packages("ggplot2", dependencies=TRUE)
install.packages("dynamicTreeCut", dependencies=TRUE)
source("https://bioconductor.org/biocLite.R")
biocLite('impute')
biocLite('topGO')
biocLite('limma')
biocLite('edgeR')

#Main analysis
cd ~/Tvedora_dev_RNAseq/input

Rscript scripts/edgeR_main.r file_name 0.05

#Main analysis results

The run produces the following files:

Chisq*.pdf - chisq test for DE, Up and Down gene subsets testing for difference in the number of genes in the sex chromosome compared to the autosomes. The expected number of genes shown is based on the number of genes on the chromosome categories being compared.

chr_location*.pdf - Outlier plots of log difference in expression for the contrast, plotted against the chromosomal location of homologues in Xenopus. FDR 5% outliers are shown in red.

de_analysis - summary table of results.

Dispersion - Dispersion plot, provides an idea on the power of the analysis.

DOWN_*.txt - gene names, logFC and FDR of the downregulated genes (significantly upregulated in the second group in the contrast name). Both 5% and 10% FDR versions are created.

FC-CPMplot.pdf - graphical representation of the DE genes in each contrast plotted agains their normalised expression.

gene*.txt - gene names of DE genes in contrast. Useful for venn diagrams to compare contrasts.

gname*.txt - alternative gene names for DE genes in contrast. Might give an idea of functions, and useful to search your favourite candidate gene.

Heatmap*.pdf - Heatmap of normalised counts in each contrast. Useful to check samples cluster as expected, and to identify interesting coexpressed gene subsets. 

logFC*.pdf - plot showing the proportion of genes in each chromosome that is DE. 

logFC_per_chr*.pdf - plot showing the logFC per chromosome, for each contrast. Stars above the chromosome name indicate the significance level of Man-Witney tests on the average logFC on that chromosome, and the remaining chromosomes.

MDS.pdf - MDS plot to check libraries make sense and how they cluster. 

Number_de.txt - Summary number of DE genes in each contrast for 5% FDR.

p_hist.pdf - Histogram of number of DE genes in each contrast. Not very useful.

pairwise_raw_count*.pdf - Pairwise correlation of gene expression amonst libraries of the same category. They should correlate well.

unbiased*.txt - List of non DE genes in each contrast.

UP_*.txt - gene names, logFC and FDR of the upregulated genes (significantly upregulated in the first group in the contrast name). Both 5% and 10% FDR versions are created.

violin_per_chr_* (DE).pdf - Violin plots of contrast expression difference, by chromosome.


### Step6. GO enrichment analysis with TopGO. 

scripts are here:
~/Tvedora_dev_RNAseq/scripts/GO_analysis_05.sh
~/Tvedora_dev_RNAseq/scripts/GO_BP.R
~/Tvedora_dev_RNAseq/scripts/GO_CC.R
~/Tvedora_dev_RNAseq/scripts/GO_MF.R

#GO analysis
#The GO analysis is run separately after the main script finishes.

~/Tvedora_dev_RNAseq
for i in $(ls | grep GO_0.05); do cd $i; source ~/Tvedora_dev_RNAseq/scripts/GO_analysis_05.sh; cd ..; done

#GO analysis results
#The main output is a table called Fisher.txt with the following columns, significant in topGO (Fisher) test are kept.


### Step7. Identify shared sex-biased or non-biased genes.

##this venn diagram can do the comparison of 5 at the maximum.
~/Tvedora_dev_RNAseq/scripts/Venn_diagram.R


### Step8. A series of R scripts to produce the result figures in this projects.

1. script to generate stacked bars.
~/Tvedora_dev_RNAseq/scripts/Stackedbars_barchat.R

2. script to display XY'/XX ratio of sex-biased genes, non-biased genes, as well as total genes.
~/Tvedora_dev_RNAseq/scripts/SB_G43G46.R
~/Tvedora_dev_RNAseq/scripts/Expression_ratio.R

### Step9. Analysis of evolutionary rate of sex-biased and unbiased genes, with one-to-one ortholog with Xenopus tropicalis.

scripts are here:
1. general pipeline to run Prank and PAMEL analysis.
~/Tvedora_dev_RNAseq/scripts/Prank_and_PAMEL.sh

2. Run PAMEL script.
~/Tvedora_dev_RNAseq/scripts/for_loop_codeml.sh

3. script to generate boxplot figure of dn/ds values.
~/Tvedora_dev_RNAseq/scripts/Dn_ds.R
